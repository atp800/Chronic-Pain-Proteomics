##########################################################
# IMPORTS 
##########################################################
import os
import sys
import warnings
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler



##########################################################
# PARAMTERS 
##########################################################
print("Scripted started")
INPUT_FILE_PATH  = "Chronic Code/AllTimePoints_PPT_Abundance_s.xlsx"
SHEET_NAME       = "Soluble"
ORIGINAL_DATA_DF = pd.read_excel(INPUT_FILE_PATH, sheet_name=SHEET_NAME)

ID_COLUMNS       = ["Sample", "Time", "%PPT (30N)", "Group", "Replicate"]  
PROTEIN_COLUMNS  = [col for col in ORIGINAL_DATA_DF.columns if col not in ID_COLUMNS]  
UNNEEDED_COLUMNS = None                                                     # use for unneeded columns

MISSINGNESS_THRESHOLD = 2                                                   # minimum number of biological samples of a proteins needed in a condition (e.g. day 14 responders) to retain it
ALREADY_LOG2 = False
ALREADY_NORMALISED = False
APPLY_OUTLIER_REMOVAL = True
OUTLIER_THRESHOLD = 3


##########################################################
# HANDLE MISSING VALUES                                                     # %PPT (30N) column: n/a values are missing, 0 values are true 0s
##########################################################                  # Protein columns: 0 values are missing (assumed MNAR)
print("Handling missing values...")

missing_val_df = ORIGINAL_DATA_DF.copy()                                    # create copy of original dataframe to handle missing values
for column in PROTEIN_COLUMNS:                                              # replace 0s with NaNs for protein columns
    missing_val_df[column] = missing_val_df[column].replace(0, np.nan)


# 1: check missing across replicates 
# - protein is considered present in sample if either replicate has non-missing value
sample_level_presence = (
    missing_val_df
    .groupby(["Sample", "Time", "Group"])[PROTEIN_COLUMNS]
    .apply(lambda g: g.notna().any())
    .reset_index()
)

# 2: check missing across each Group and each Time
presence_counts = (
    sample_level_presence
    .groupby(["Group", "Time"])[PROTEIN_COLUMNS]
    .sum()                                          # True→1, False→0
)

# 3. define planned comparisons
comparisons = [
    (("Resp", "T0"), ("Resp", "D14")),              # day 0 vs 14 responders
    (("Resp", "D14"), ("NonResp", "D14")),          # day 14 responders vs non-responders
]

# 4. for each comparison, check for at least two biological samples
# - union, so proteins are kept if present enough for either comparison
min_samples = MISSINGNESS_THRESHOLD

keep_mask = pd.Series(False, index=PROTEIN_COLUMNS)                                                     # create a boolean mask for proteins to keep - initialise as false for all

for c in comparisons:
    cond1, cond2 = c
    if cond1 not in presence_counts.index or cond2 not in presence_counts.index:                        # if the condition is not in the dataset, skip it to avoid errors
        continue
    mask = (presence_counts.loc[cond1] >= min_samples) & (presence_counts.loc[cond2] >= min_samples)    # mark proteins to keep for current comparison in 2nd mask
    keep_mask |= mask                                                                                   # union: update the proteins to keep in main mask with the 2nd one

# Final list of proteins
proteins_removed = keep_mask[~keep_mask].index.tolist()
proteins_to_keep = keep_mask[keep_mask].index.tolist()
filtered_df = missing_val_df[ID_COLUMNS + proteins_to_keep]

print(f"Kept proteins after missingness filter: {len(proteins_to_keep)}/{len(PROTEIN_COLUMNS)}")
print(f"Removed proteins: {len(proteins_removed)}")



##########################################################
# TRANSFORMATIONS
########################################################## 
print("Applying transformations...")

if not ALREADY_NORMALISED:          # Median-centered normalisation

    normalised_df = filtered_df.copy()
    protein_data = normalised_df[proteins_to_keep]
    sample_medians = protein_data.median(axis=1)
    grand_median = sample_medians.median()
    norm_factors = grand_median / sample_medians
    normalized_protein_data = protein_data.mul(norm_factors, axis=0)
    normalised_df[proteins_to_keep] = normalized_protein_data

    filtered_df = normalised_df

if not ALREADY_LOG2:

    transformed_df = filtered_df.copy()
    for column in proteins_to_keep:
        transformed_df[column] = np.log2(transformed_df[column])
    
    filtered_df = transformed_df


##########################################################
# REMOVE EXTREME OUTLIERS
##########################################################
if APPLY_OUTLIER_REMOVAL:
    print(f"Checking for outliers...")
    
    calc_data = filtered_df[proteins_to_keep].copy()

    if not ALREADY_LOG2:
        # Suppress warnings for log(0) if any remain, though we replaced 0s with NaNs earlier
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            calc_data = np.log2(calc_data)

    z_scores = (calc_data - calc_data.mean()) / calc_data.std()

    is_outlier = np.abs(z_scores) > OUTLIER_THRESHOLD


    total_outliers = is_outlier.sum().sum()
    filtered_df.loc[:, proteins_to_keep] = filtered_df.loc[:, proteins_to_keep].mask(is_outlier)

    print(f"Outlier removal complete.")
    print(f"Total data points detected as outliers and set to NaN: {total_outliers}")

    total_nans = filtered_df[proteins_to_keep].isna().sum().sum()
    print(f"Total missing values (NaNs) in dataset now: {total_nans}")




