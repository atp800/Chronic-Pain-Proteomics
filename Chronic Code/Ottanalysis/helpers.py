import numpy as np
import pandas as pd
import sys
from word2number import w2n
import re
import warnings


##################################################################
# HELPER FUCNTION: IMPUTATE MISSING VALUES
##################################################################
# Probabilistic minimum imputation (down-shifted normal distribution)
def impute_missing_vals(df, shift=1.8, width=0.3):
    """
    Imputes missing values using down-shifted normal distribution    
    Assumes missing values are MNAR (abundance beyond limit of detection)
    
    Parameters:
    - shift: How many standard deviations to shift the distribution left (default 1.8)
    - width: The width of the new noise distribution relative to original std (default 0.3)
    """

    # Set a random seed for reproducibility and create df copy to avoid modifying original
    np.random.seed(938)
    data = df.copy()
    
    # Iterate through each protein column
    for col in data.columns:
        # Skip if no missing values
        if data[col].isna().sum() == 0:
            continue
            
        # Get statistics of valid values
        valid_data = data[col].dropna()
        mu = valid_data.mean()
        sigma = valid_data.std()
        
        # If sigma is 0 (all values same) or NaN (1 value), fallback to 0 or min
        if pd.isna(sigma) or sigma == 0:
            data[col] = data[col].fillna(0)
            continue
            
        # Calculate parameters for the noise distribution
        # Shift the mean down by 1.8 std devs
        impute_mean = mu - (shift * sigma)
        # We make the spread narrower (0.3 of original)
        impute_std = width * sigma
        
        # Generate random numbers for the missing entries
        n_missing = data[col].isna().sum()
        noise_values = np.random.normal(loc=impute_mean, scale=impute_std, size=n_missing)
        
        # Fill NaN spots with these values
        # Use a mask to fill only the NaNs
        mask = data[col].isna()
        data.loc[mask, col] = noise_values
        
    return data



##################################################################
# HELPER FUNCTION: GROUP-WISE MISSING VALUE FILTERING
##################################################################
# FIlters missing values but retains if only missing in one group
def groupwise_missing_filter(df_proteins, series_groups, threshold=0.7):
    """
    Keeps proteins that are present in at least 'threshold' (e.g., 0.7 = 70%)
    of samples in AT LEAST ONE experimental group
    
    Preserves on/off biomarkers (e.g. present in disease but not in control)
    """


    print(f"Applying Group-Based Missingness Filter (Threshold: {threshold:.0%})")
    initial_count = df_proteins.shape[1]
    
    unique_groups = series_groups.unique()
    proteins_to_keep = []

    for protein in df_proteins.columns:
        keep_this_protein = False
        for group in unique_groups:
            # Get values for this specific group
            group_vals = df_proteins.loc[series_groups == group, protein]
            
            # Calculate valid ratio (count() excludes NaNs)
            if len(group_vals) > 0:
                valid_ratio = group_vals.count() / len(group_vals)
                if valid_ratio >= threshold:
                    keep_this_protein = True
                    break # Found a group where it's valid, so we keep it

        if keep_this_protein:
            proteins_to_keep.append(protein)

    filtered_df = df_proteins[proteins_to_keep]
    dropped_count = initial_count - filtered_df.shape[1]
    print(f"Dropped {dropped_count} proteins. Remaining: {filtered_df.shape[1]}")
    
    return filtered_df




##################################################################
# HELPER FUNCTION: DERIVE SUBJECT ID
##################################################################
def derive_subject_ids(sample_id_list, delimiter):
    """
    Extracts the subject identifier from a sample ID.
    Assumes the format "SubjectID-OtherInfo" (where '-' is the delimiter specified by the ID_DELIMITER parameter).
    """
    subject_ids = [str(s).split(delimiter)[0] for s in sample_id_list]

    # Check for missing delimiters by comapring input and output per item
    failures = [str(orig) for orig, new in zip(sample_id_list, subject_ids) if str(orig) == new]
    failure_count = len(failures)

    if failure_count > 0:
        print("\n" + "!"*65)
        print(f"WARNING: Subject ID derivation issue!!!!!")
        print(f"The delimiter '{delimiter}' was NOT found in {failure_count} out of {len(sample_id_list)} samples.")
        print(f"These samples retained their full ID and cannot be paired correctly.")
        
        # Show the first 3 examples of bad IDs so the user can debug
        if failure_count == len(sample_id_list):
            print("CRITICAL: Delimiter not found in ANY sample. Check 'ID_DELIMITER' setting.")
        else:
            print(f"Example IDs missing delimiter: {failures[:3]}")
        print("!"*65 + "\n")

    return subject_ids




##################################################################
# CALCULATE DEALTAS FOR DELTA MODE
##################################################################
def calculate_deltas(df, id_col, time_col, condition_col, protein_cols, id_delimiter, baseline_time, compare_time):
    '''
    If selected in GUI, ENTIRE dataset is converted to difference scores based on selected DELTA_COLUMN
    All subsequent tests (Limma, LogReg etc.) will run on the change values
    '''
    
    print("\n" + "="*65)
    print(f"TRANSFORMING DATA: Calculating Deltas ({compare_time} - {baseline_time})")
    print("="*65)

    subjects = derive_subject_ids(df[id_col].tolist(), delimiter=id_delimiter)
    
    df_work = df.copy()
    df_work[protein_cols] = impute_missing_vals(df_work[protein_cols])
    df_work['Temp_Subject_ID'] = subjects

    try:
        if df_work.duplicated(subset=['Temp_Subject_ID', time_col]).any():                          
            print("Error: Duplicate samples found (Same Subject + Same Timepoint). Cannot pivot.")
            sys.exit(1)

        # Pivot: Index=Subject, Columns=Time, Values=Proteins
        df_pivot = df_work.pivot(index='Temp_Subject_ID', columns=time_col, values=protein_cols) 
        
        vals_end = df_pivot.xs(compare_time, axis=1, level=time_col)
        vals_start = df_pivot.xs(baseline_time, axis=1, level=time_col)
        df_delta = vals_end - vals_start
        
        n_original = len(df_pivot)
        n_final = len(df_delta)
        print(f"Subjects with complete pairs: {n_final} (Dropped {n_original - n_final})")

        # Restore group data
        if condition_col:
            subj_group_map = df_work.drop_duplicates('Temp_Subject_ID').set_index('Temp_Subject_ID')[condition_col]
            df_delta[condition_col] = df_delta.index.map(subj_group_map)
            
        df_delta[id_col] = df_delta.index
        df_delta.reset_index(drop=True, inplace=True)
        
        return df_delta

    except Exception as e:
        print(f"Error during delta calculation: {e}")
        sys.exit(1)




##################################################################
# CONVERT TIMEPOINTS TO NUMBERS (FOR HINGE MODEL)
##################################################################
def convert_timepoint_to_number(original_timepoint):
    '''
    Converts a value in time column to a numerical floating point value
    E.g. five->5, D14->14, T0->0, day 5->5 etc.
    '''
    
    original_timepoint = original_timepoint.strip()         # Remove leading/trailing whitespace
    digits = re.search(r'-?\d+\.?\d*', original_timepoint)    # Look for digits
    numeric_value = None

    if digits:
        numeric_value = float(digits.group(0))
    else:
        try:
            numeric_value = w2n.word_to_num(original_timepoint)   # Convert number words (e.g. five)
        except ValueError:
            numeric_value = None

    if numeric_value is None:
        warnings.warn(f"Could not convert time label '{original_timepoint}' to a number")
    else:
        numeric_value = float(numeric_value)

    return numeric_value