import os
import sys
import pandas as pd
from config import CONFIG
from gui import run_gui

from helpers import calculate_deltas
from limma import run_limma
from logistic_regression import run_logistic_regression


# =========================================================
# 1. SETUP & CONFIGURATION
# =========================================================
if CONFIG.get("LAUNCH_GUI", True):
    run_gui()

# # Check for file
# if not CONFIG.get("INPUT_FILE_PATH") or not os.path.exists(CONFIG["INPUT_FILE_PATH"]):
#     sys.exit("Error: input file required")

# Setting print statements
MODE = CONFIG.get("ANALYSIS_MODE")
print("\n" + "="*60)
print(f"ANALYSIS MODE: {MODE}")
print("="*60)

print(f"Baseline: {CONFIG['BASELINE_VAL']}")
print(f"Comparisons to make: {CONFIG['COMPARE_VALS']}")
if CONFIG['LOOP_FILTER_VALS']:
    print(f"Will loop independently for: {CONFIG['LOOP_FILTER_VALS']}")
else:
    print("Will run on entire dataset (No subgroup loops selected).")

# Load data into dataframe
print("\nLoading Data...")
df_main = pd.read_excel(CONFIG["INPUT_FILE_PATH"], sheet_name=CONFIG["SHEET_NAME"])

# Determine protein columns (excludes id/condition/time columns and any other unneeded columns)
exclude_cols = {CONFIG["ID_COLUMN"], CONFIG["CONDITION_COLUMN"], CONFIG["TIME_COLUMN"]} | set(CONFIG["UNNEEDED_COLUMNS"])
exclude_cols = {c for c in exclude_cols if c is not None} # removes any none values
protein_cols =[c for c in df_main.columns if c not in exclude_cols]

# =========================================================
# 2. ANALYSES CALL FUNCTION
# =========================================================
def run_statistical_tests(df_subset, group_col, baseline_name, compare_name, output_dir, is_paired, is_delta):
    """
    This function expects a dataframe with EXACTLY 2 conditions or groups
    Insert Limma, Logistic Regression, etc., here.
    """
    print(f"\nRunning comparison: {baseline_name} vs {compare_name}")
    print(f"\nRows in subset: {len(df_subset)}")
    print(f"\nPaired Limma: {is_paired} | Is Delta Mode: {is_delta}")
    
    # Check output folder exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Delta transformation if used
    if is_delta:
        df_subset = calculate_deltas(
            df=df_subset, 
            id_col=CONFIG["ID_COLUMN"], 
            time_col=CONFIG["TIME_COLUMN"], 
            condition_col=CONFIG["CONDITION_COLUMN"], 
            protein_cols=protein_cols, 
            id_delimiter=CONFIG["ID_DELIMITER"], 
            baseline_time=baseline_name, 
            compare_time=compare_name
        )
        group_col = CONFIG["CONDITION_COLUMN"] # Once delta is calculated, Limma groups by Condition (Responder vs NonResp)
        is_paired = False # Deltas merge paired rows into single values

    # Run tests
    if CONFIG.get("RUN_LIMMA", False):
        run_limma(df_subset, group_col, protein_cols, baseline_name, compare_name, output_dir, is_paired, is_delta, CONFIG)
        
    if CONFIG.get("RUN_LOGISTIC_REGRESSION", False):
        run_logistic_regression(df_subset, group_col, protein_cols, baseline_name, compare_name, output_dir)

# =========================================================
# Loops for running analyses
# =========================================================

# If user left the "Loop" listbox empty, add a dummy value to let the loop run once
outer_loops = CONFIG["LOOP_FILTER_VALS"] if CONFIG["LOOP_FILTER_VALS"] else [None]
base_val = CONFIG["BASELINE_VAL"]
comp_vals = CONFIG["COMPARE_VALS"]

for current_filter in outer_loops:
    
    # --- SUBSET DATA BY OUTER LOOP (e.g. Filter to just "D0" or just "Responders") ---
    df_filtered = df_main.copy()
    output_subfolder = CONFIG["OUTPUT_FOLDER"]
    
    if current_filter is not None:
        print(f"\n--- Starting Analysis block for: {current_filter} ---")
        output_subfolder = os.path.join(CONFIG["OUTPUT_FOLDER"], str(current_filter))
        
        # Determine what column we are filtering on based on the Mode
        if MODE == "GROUP_COMPARISON":
            df_filtered = df_filtered[df_filtered[CONFIG["TIME_COLUMN"]].astype(str) == str(current_filter)]
        elif MODE in ["LONGITUDINAL", "DELTA"]:
            df_filtered = df_filtered[df_filtered[CONFIG["CONDITION_COLUMN"]].astype(str) == str(current_filter)]


    # --- INNER LOOP: Iterate through the targeted comparisons ---
    for target_val in comp_vals:
        
        # Setup the 2-condition dataframe
        if MODE == "GROUP_COMPARISON":
            # Keep only Baseline Group and Target Group
            col_to_filter = CONFIG["CONDITION_COLUMN"]
            df_2_groups = df_filtered[df_filtered[col_to_filter].astype(str).isin([str(base_val), str(target_val)])]
            
            run_statistical_tests(
                df_subset=df_2_groups, 
                group_col=CONFIG["CONDITION_COLUMN"], 
                baseline_name=base_val, 
                compare_name=target_val, 
                output_dir=output_subfolder,
                is_paired=False, 
                is_delta=False
            )

        elif MODE == "LONGITUDINAL":
            # Keep only Baseline Time and Target Time
            col_to_filter = CONFIG["TIME_COLUMN"]
            df_2_times = df_filtered[df_filtered[col_to_filter].astype(str).isin([str(base_val), str(target_val)])]
            
            run_statistical_tests(
                df_subset=df_2_times, 
                group_col=CONFIG["TIME_COLUMN"], 
                baseline_name=base_val, 
                compare_name=target_val, 
                output_dir=output_subfolder,
                is_paired=True,     # Longitudinal assumes paired tracking!
                is_delta=False
            )

        elif MODE == "DELTA":
            # Keep only Baseline Time and Target Time
            col_to_filter = CONFIG["TIME_COLUMN"]
            df_2_times = df_filtered[df_filtered[col_to_filter].astype(str).isin([str(base_val), str(target_val)])]
            
            run_statistical_tests(
                df_subset=df_2_times, 
                group_col=CONFIG["CONDITION_COLUMN"], # Stats group by Condition (Responder)
                baseline_name=base_val, 
                compare_name=target_val, 
                output_dir=output_subfolder,
                is_paired=False,    # Deltas collapse pairs into independent change values
                is_delta=True
            )

print("\n-----------------------------------")