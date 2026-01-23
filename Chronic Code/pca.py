import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from statsmodels.multivariate.manova import MANOVA

# ==========================================
# 1. SETUP
# ==========================================
# File Paths
INPUT_FILE_PATH  = "Chronic Code/Cleaned_Data - Soluble.xlsx"                         # NEED TO ADD BACK IN GROUP AND TIME COLUMNS (could get time from sample ID)
SHEET_NAME       = "Sheet 1"
DATA = pd.read_excel(INPUT_FILE_PATH, sheet_name=SHEET_NAME)


df = pd.DataFrame(DATA)

# ==========================================
# 2. PREPROCESSING
# ==========================================

# Extract Time and Subject from Sample_ID
# (Assumes format "Subject-Time")
df['Time'] = df['Sample_ID'].apply(lambda x: x.split('-')[1])
df['Subject'] = df['Sample_ID'].apply(lambda x: x.split('-')[0])

df['Study'] = df['Sample_ID'].apply(lambda x: x.split('P')[0])

# Separate Metadata from Expression Data
metadata_cols = ['Study','Sample_ID', 'Group', 'Replicate', 'Time', 'Subject']
expression_cols = [c for c in df.columns if c not in metadata_cols]

# # Logify
# epsilon = 1e-9 
# df_log = df.copy()
# df_log[expression_cols] = np.log2(df[expression_cols] + epsilon)
#df = df_log

# Impute missing vals
shift = 1.8
width = 0.3
np.random.seed(938) # For reproducibility

imputed_data = df[expression_cols].copy()       # Create a working copy of just the expression data

for col in imputed_data.columns:
    # Skip if no missing values
    if imputed_data[col].isna().sum() == 0:
        continue
    
    valid_data = imputed_data[col].dropna()         # Get statistics of valid values
    mu = valid_data.mean()
    sigma = valid_data.std()
    
    if pd.isna(sigma) or sigma == 0:                # If std is NaN or zero, fill with zeros
        imputed_data[col] = imputed_data[col].fillna(0)
        continue
        
    impute_mean = mu - (shift * sigma)              # Calculate mean and std for imputation
    impute_std = width * sigma
    
    n_missing = imputed_data[col].isna().sum()      # Generate random numbers for the missing entries
    noise_values = np.random.normal(loc=impute_mean, scale=impute_std, size=n_missing)

    mask = imputed_data[col].isna()                 # Fill NaN spots
    imputed_data.loc[mask, col] = noise_values

df[expression_cols] = imputed_data

# Standardisation (Z-score)
scaler = StandardScaler()
scaled_data = scaler.fit_transform(df[expression_cols])


# ==========================================
# 2.5 Determine number of principle components
# ==========================================
# Run PCA without specifying the number of components to see how much variance each one explains
pca_full = PCA().fit(scaled_data)

# Calculate the cumulative explained variance
cumulative_variance = np.cumsum(pca_full.explained_variance_ratio_)

# Find the number of components to explain 95% of the variance
n_components_95_variance = np.where(cumulative_variance >= 0.95)[0][0] + 1
print(f"Number of components to explain 95% of variance: {n_components_95_variance}")


# ==========================================
# 3. PCA (Principal Component Analysis)
# ==========================================

# Run PCA
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(scaled_data)
pca_df = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])

# Combine PCA results with metadata for plotting
final_df = pd.concat([pca_df, df[metadata_cols].reset_index(drop=True)], axis=1)

# Calculate explained variance
exp_var = pca.explained_variance_ratio_
print(f"Explained Variance: PC1={exp_var[0]*100:.2f}%, PC2={exp_var[1]*100:.2f}%")

# ==========================================
# 4. VISUALISATION
# ==========================================

plt.figure(figsize=(16, 12))

# Plot 1: Coloured by Biological Factor (Group)
plt.subplot(2, 2, 1)
sns.scatterplot(
    data=final_df, x='PC1', y='PC2', 
    hue='Group', style='Time', s=100
)
plt.title('PCA: Biological (Group)')
plt.grid(True, alpha=0.3)

# Plot 2: Coloured by Study (Technical 1)
plt.subplot(2, 2, 2)
sns.scatterplot(
    data=final_df, x='PC1', y='PC2', 
    hue='Study', style='Time', 
    palette='Set1', s=100
)
plt.title('PCA: Technical Check (Study)')
plt.grid(True, alpha=0.3)

# Plot 3: Coloured by Replicate (Technical 2)
plt.subplot(2, 2, 3)
sns.scatterplot(
    data=final_df, x='PC1', y='PC2', 
    hue='Replicate', style='Time', 
    palette='Set2', s=100
)
plt.title('PCA: Technical Check (Replicate)')
plt.grid(True, alpha=0.3)

# Plot 4: Coloured by Time (Biological/Longitudinal)
plt.subplot(2, 2, 4)
sns.scatterplot(
    data=final_df, x='PC1', y='PC2', 
    hue='Time', style='Group', 
    palette='viridis', s=100
)
plt.title('PCA: Longitudinal Check (Time)')
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# ==========================================
# 5. CLUSTERMAP (Hierarchical Clustering)
# ==========================================

# # Calculate correlation matrix
# corr_matrix = df_log[expression_cols].T.corr(method='pearson')


# # 1. Group Colors
# group_colors = df_log['Group'].map({'Resp': 'blue', 'NonResp': 'red'})

# # 2. Replicate Colors (Manual mapping)
# replicate_colors = df_log['Replicate'].map({1: 'green', 2: 'orange'})

# # 3. Study Colors (Automated mapping using a palette)
# unique_studies = df_log['Study'].unique()
# study_palette = sns.color_palette("husl", len(unique_studies))
# study_lut = dict(zip(unique_studies, study_palette))
# study_colors = df_log['Study'].map(study_lut)

# # Combine into a dataFrame for the row/col colors
# row_colors_df = pd.DataFrame({
#     'Group': group_colors,
#     'Replicate': replicate_colors,
#     'Study': study_colors
# })

# # Plot
# sns.clustermap(
#     corr_matrix, 
#     row_colors=row_colors_df,
#     col_colors=row_colors_df,
#     cmap='viridis',
#     figsize=(10, 10),
#     annot=False # I turned off annotations (numbers) as they are messy on large plots
# )
# plt.title("Sample Correlation Matrix")
# plt.show()




# ==========================================
# 6. MONVOA TEST
# ==========================================
print("\n" + "="*40)
print("MANOVA RESULTS (PCs ~ Factor)")
print("Tests for significant separation in PCA")
print("="*40)

# Testing if the PC1/PC2 combination differs significantly within each comparison
# P threshold 0.05
# For group (rep/nonresp group)/time: low p-value = good (disparity between groups)
# For replicate/study: high p-value = good (no batch effect)

factors_to_test = ['Group', 'Study', 'Replicate', 'Time']
results = []

# Ensure categorical types for the model
final_df['Replicate'] = final_df['Replicate'].astype(str)
final_df['Time'] = final_df['Time'].astype(str)

for factor in factors_to_test:
    # Check if the factor exists and has at least 2 levels (can't test if only 1 group)
    if factor in final_df.columns and len(final_df[factor].unique()) > 1:
        
        # Create the formula
        formula = f'PC1 + PC2 ~ {factor}'
        
        # Run MANOVA
        maov = MANOVA.from_formula(formula, data=final_df)
        test_res = maov.mv_test()
        
        # Wilks' Lambda p-value
        p_val = test_res.results[factor]['stat'].loc["Wilks' lambda", 'Pr > F']
        
        results.append({
            'Factor': factor, 
            'P-Value': p_val, 
            'Significant?': 'Yes' if p_val < 0.05 else 'No'
        })
    else:
        results.append({
            'Factor': factor, 
            'P-Value': 'N/A (Singular/Missing)', 
            'Significant?': '-'
        })

# Create a dataFrame for display
stats_df = pd.DataFrame(results)
pd.options.display.float_format = '{:.4g}'.format

print(stats_df)
print("\nInterpretation:")
print("- Group/Time: Significant (Yes) = Good (group disparity)")
print("- Study/Replicate:     Significant (Yes) = Bad (batch effect)")
print("="*40)