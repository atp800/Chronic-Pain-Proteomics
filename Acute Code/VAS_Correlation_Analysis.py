import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
import math


'''
This script checks for correlations between protein levels and VAS (pain score)
It checks for individual participants, and then across all participants
It exports significant correlations to an Excel file
'''


FILE_PATH = 'Combined_Data.xlsx'
METADATA_COLS = ['Sample ID', 'Condition']
VAS_COL = 'VAS (Pain Score)'
P_THRESHOLD = 0.05
REALTIVE_RANGE_THRESHOLD = 0.2  # Minimum relative range (20%) for a protein to be considered


# Read Excel file
df = pd.read_excel(FILE_PATH, sheet_name='Sheet1')



# ------------------------------
# Cleaning                      - THIS IS PAULINE STUDY SPECIFIC
# ------------------------------ 
# Remove rows with no protein data
protein_cols = [c for c in df.columns if c not in METADATA_COLS and c != VAS_COL]
df_clean = df.dropna(subset=protein_cols, how='all')


# # Remove any protein columns with >=3 missing values                                  --- Now redundant, checking for all 4 datapoints instead in correlation method (if len(x) < 4 or x.nunique() < 2:)
# protein_cols = [c for c in df.columns if c not in METADATA_COLS and c != VAS_COL]
# print(f"Proteins before dropping missing data: {len(protein_cols)}")

# filtered_proteins = [c for c in protein_cols if df[c].isna().sum() < 3]
# df_clean = df[METADATA_COLS + [VAS_COL] + filtered_proteins]
# print(f"Proteins remaining after dropping missing data: {len(filtered_proteins)}")

# protein_cols = filtered_proteins
# ------------------------------



# ------------------------------
# Method for computing correlations
# ------------------------------ 
all_results = []
def compute_correlations(sub_df, label):
    results = []

    for protein in protein_cols:
        # Create a temp pair for the protein and particpant and drop missing values
        df_pair = sub_df[[protein, VAS_COL]].dropna()
        x, y = df_pair[protein], df_pair[VAS_COL]

        # Check for all 4 datapoints, at least 2 unique values, and no zeros
        if len(x) < 4 or x.nunique() < 2 or (x == 0).any():
            continue

        # Check relative range is over threshold (protein variation isn't noise)
        if (x.max() - x.min()) / x.min() > REALTIVE_RANGE_THRESHOLD:  # Skip proteins with low relative range
            continue


        r, p = pearsonr(x, y)

        # Check if correlation is significant
        #if p < P_THRESHOLD and abs(r) < 0.998:                  # Exclude perfect correlations (caused by missing data)  #CHANGE FROM EXCLUSION TO A BOOLEAN COLUMN
        Significant = False
        if p < P_THRESHOLD:
            Significant = True
                                                
        results.append({'Protein': protein,
                        'Correlation': r,
                        'p-value': p,
                        'Label': label,
                        'Significant': Significant})                    # HOW IS P VALUE CALCULATED FOR INDIVIDUAL PARTICIPANTS

    # Print results
    print(f"\n=== Significant correlations for {label} ===")
    if not results:
        print("No significant correlations.")
    else:
        for res in sorted(results, key=lambda x: x['p-value']):
            print(f"{res['Protein']}: r={res['Correlation']:.3f}, p={res['p-value']:.4f}")

    # Add to global results
    all_results.extend(results)



# ------------------------------
# Correlations per individual participant
# ------------------------------
for pid, sub_df in df_clean.groupby('Sample ID'):
    compute_correlations(sub_df, label=f"participant {pid}")


# ------------------------------
# Correlations across all participants
# ------------------------------
compute_correlations(df_clean, label="ALL PARTICIPANTS")



# ------------------------------
# Export reuslts to Excel
# ------------------------------
if all_results:
    df_results = pd.DataFrame(all_results)
    df_results.fillna('n/a', inplace=True)
    df_results['Significance_Frequency'] = df_results.groupby('Protein')['Significant'].transform('sum')                               # Get frequency each significant protein appears
    df_results = df_results.sort_values(by=['Significance_Frequency', 'Protein', 'Correlation'], ascending=[False, True, False])     # Sort by frequency, then protein name, then correlation
    df_results.to_excel('Significant_VAS_Correlations.xlsx', index=False)
    print("\nSaved significant correlations to Significant_Correlations.xlsx")
else:
    print("\nNo significant correlations found; nothing to save.")






# ------------------------------
# Plot top X strongest correlations
# ------------------------------
NUMBER_TO_PLOT = 16

if all_results:
    # Sort by absolute correlation value
    top_results = sorted(all_results, key=lambda x: abs(x['Correlation']), reverse=True)[:NUMBER_TO_PLOT]

    # Determine subplot grid size
    n_cols = 4
    n_rows = math.ceil(len(top_results) / n_cols)

    # Set figure size
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*5.5, n_rows*4.5))
    axes = axes.flatten()  # flatten in case of single row/column

    for ax, res in zip(axes, top_results):
        protein = res['Protein']
        label = res['Label']
        r = res['Correlation']
        p = res['p-value']

        # Subset for plotting
        if label.startswith("participant"):
            pid = label.split()[-1]
            sub_df = df_clean[df_clean['Sample ID'] == pid]
        else:
            sub_df = df_clean

        plot_data = sub_df[[protein, VAS_COL]].dropna()             # Don't plot missing values
        x = sub_df[protein]
        y = sub_df[VAS_COL]

        sns.scatterplot(x=x, y=y, ax=ax)
        sns.regplot(x=x, y=y, scatter=False, color='red', ax=ax)    # optional regression line
        ax.set_ylim(-1, 11)                                         # Restrict to VAS scale (1-10) with some padding
        ax.set_title(f"{label}\nr={r:.4f}, p={p:.3f}", fontsize=9)
        ax.set_xlabel(protein, fontsize=9)
        ax.set_ylabel("VAS", fontsize=9)

    # Hide unused axes
    for i in range(len(top_results), len(axes)):
        fig.delaxes(axes[i])

    # Add spacing between plots
    plt.subplots_adjust(hspace=0.5, wspace=0.3, top=0.92)
    plt.show()
else:
    print("No significant correlations to plot.")
