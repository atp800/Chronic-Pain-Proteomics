import pandas as pd


'''
- This script preprocess the proteomics data for Pauline's pain study
- It ingests the original Excel sheets and outputs a new Excel file combining all the data
- If the file (COMBINED_Data.xlsx) already exists, it will be overwritten
- It is specific to that data, and isn't efficient enoguh for large datasets
- To run, just put the source data in the same folder as the script and run the file. The output will appear in the same folder
'''


FILE_PATH = 'Acute_Raw-Abundance_VAS.xlsx'
METADATA_COLS = ['Sample ID', 'Condition']


# Read Excel files
df_vas = pd.read_excel(FILE_PATH, sheet_name='VAS')
df_soluble = pd.read_excel(FILE_PATH, sheet_name='Fraction_s')
df_SPE = pd.read_excel(FILE_PATH, sheet_name='Fraction_SPE')
df_insoluble = pd.read_excel(FILE_PATH, sheet_name='Fraction_i')


# Remove suffixes from sample IDs
df_insoluble['Sample ID'] = df_insoluble['Sample ID'].str.replace('_i$', '', regex=True)
df_SPE['Sample ID'] = df_SPE['Sample ID'].str.replace('_b$', '', regex=True)


# Drop Group column
for df in [df_soluble, df_SPE, df_insoluble]:
    df.drop(columns=['Group'], inplace=True)


# Drop replicate column and average values across replicates                        # FIll 0 as N/A - if N/A for one replicate only take existing value (don't average)
def collapse_replicates(df):
    df = df.copy()

    if 'Replicate' in df.columns:
        df = df.drop(columns=['Replicate'])

    id_cols = METADATA_COLS                                                         # Identify metadata columns

    df_avg = (                                                                      # Group + average ONLY protein columns
        df.groupby(id_cols, as_index=False)
          .mean(numeric_only=True)
    )
    return df_avg

df_soluble_avg = collapse_replicates(df_soluble)
df_SPE_avg = collapse_replicates(df_SPE)
df_insoluble_avg = collapse_replicates(df_insoluble)


# Identify metadata columns and get full set of protein columns
id_cols = METADATA_COLS

cols_sol = set(df_soluble_avg.columns) - set(id_cols)
cols_spe = set(df_SPE_avg.columns) - set(id_cols)
cols_i   = set(df_insoluble_avg.columns) - set(id_cols)


# Find duplicated protein columns and add sufixes
all_cols = list(cols_sol | cols_spe | cols_i)
duplicates = {
    c for c in all_cols 
    if sum([c in cols_sol, c in cols_spe, c in cols_i]) > 1
}

def add_suffix_only_if_duplicate(df, suffix, duplicates):                           # method for adding suffixes
    df = df.copy()
    df.rename(columns={c: f"{c}_{suffix}" for c in df.columns if c in duplicates}, inplace=True)
    return df

df_soluble_avg = add_suffix_only_if_duplicate(df_soluble_avg, "s", duplicates)      # applying suffixes
df_SPE_avg = add_suffix_only_if_duplicate(df_SPE_avg, "spe", duplicates)
df_insoluble_avg = add_suffix_only_if_duplicate(df_insoluble_avg, "i", duplicates)


# Merge datframes on metadata columns
dfs = [df_vas, df_soluble_avg, df_SPE_avg, df_insoluble_avg]

df_combined = dfs[0]
for df in dfs[1:]:
    df_combined = df_combined.merge(df, on=id_cols, how='outer')


# Remove time suffix from sample IDs
df_combined['Sample ID'] = df_combined['Sample ID'].str.replace(r'-T\d+$', '', regex=True)


# Sort protein columsn in alphabetical order
vas_col = ['VAS (Pain Score)']
protein_cols = sorted([c for c in df_combined.columns if c not in METADATA_COLS and c not in vas_col])
df_combined = df_combined[METADATA_COLS + vas_col + protein_cols]


# Export combined dataframe to Excel
df_combined.to_excel('Combined_Data.xlsx', index=False)
print("hooray!")