import os
import sys
import pandas as pd
from config import CONFIG
from gui import run_gui

from helpers import calculate_deltas
from limma import run_limma
from logistic_regression import run_logistic_regression
from hinge_model import run_hinge_model


# =========================================================
# 1. SETUP & CONFIGURATION
# =========================================================
if CONFIG.get("LAUNCH_GUI", True):
    run_gui()

# # Check for file
# if not CONFIG.get("INPUT_FILE_PATH") or not os.path.exists(CONFIG["INPUT_FILE_PATH"]):
#     sys.exit("Error: input file required")

# Settings print statements
MODE = CONFIG.get("ANALYSIS_MODE")
print("\n" + "="*60)
print(f"ANALYSIS MODE: {MODE}")
print("="*60)

if MODE == "DELTA":
    print(f"Delta Calculation: Change from baseline {CONFIG.get('DELTA_BASELINE')} to {CONFIG.get('LOOP_VALS')}")

print(f"Baseline: {CONFIG['BASELINE_VAL']}")
print(f"Comparisons to make: {CONFIG['COMPARE_VALS']}")

if CONFIG['LOOP_VALS'] and CONFIG['LOOP_VALS'] != [None]:
    print(f"Will loop independently for: {CONFIG['LOOP_VALS']}")
else:
    print("Will run on entire dataset (no separate loops)")


# Load data into dataframe
print("\nLoading data...")
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
    """
    print(f"Running comparison: {baseline_name} vs {compare_name}")
    print(f"Rows in subset: {len(df_subset)}")
    print(f"Paired Limma: {is_paired} | Is Delta Mode: {is_delta}")
    
    # Check output folder exists
    os.makedirs(output_dir, exist_ok=True)

    # Run tests
    if CONFIG.get("RUN_LIMMA", False):
        run_limma(df_subset, group_col, protein_cols, baseline_name, compare_name, output_dir, is_paired, is_delta, CONFIG)
        
    if CONFIG.get("RUN_LOGISTIC_REGRESSION", False):
        run_logistic_regression(df_subset, group_col, protein_cols, baseline_name, compare_name, output_dir)

# =========================================================
# 3. Loops for running limma/logistic analyses
# =========================================================

# If user left the loop vals selection empty, add a dummy value to let the loop run once
outer_loops = CONFIG["LOOP_VALS"] if CONFIG["LOOP_VALS"] else [None]
base_val = CONFIG["BASELINE_VAL"]
comp_vals = CONFIG["COMPARE_VALS"]

for current_loop_val in outer_loops:
    
    # --- SUBSET DATA BY OUTER LOOP (e.g. Filter to just "D0" or just "Responders") ---
    df_filtered = df_main.copy()
    mode_folder = os.path.join(CONFIG["OUTPUT_FOLDER"], MODE)       # set folder for results based on analysis mode
    print("\n------------------------------------------------------------------")
    print("------------------------------------------------------------------")
    print("------------------------------------------------------------------")
    
    if current_loop_val is not None:
        print(f"\nStarting analysis for: {current_loop_val}")
        output_subfolder = os.path.join(mode_folder, str(current_loop_val))
        
        # Determine what column we are filtering on based on the Mode
        if MODE == "GROUP_COMPARISON":
            df_filtered = df_filtered[df_filtered[CONFIG["TIME_COLUMN"]].astype(str) == str(current_loop_val)]
        elif MODE == "LONGITUDINAL":
            df_filtered = df_filtered[df_filtered[CONFIG["CONDITION_COLUMN"]].astype(str) == str(current_loop_val)]
        elif MODE == "DELTA":
            # Calculate deltas if in delta mode
            print(f"Transforming data into deltas ({CONFIG['DELTA_BASELINE']} to {current_loop_val})")
            df_filtered = calculate_deltas(
                df=df_filtered, 
                id_col=CONFIG["ID_COLUMN"], 
                time_col=CONFIG["TIME_COLUMN"], 
                condition_col=CONFIG["CONDITION_COLUMN"], 
                protein_cols=protein_cols, 
                id_delimiter=CONFIG.get("ID_DELIMITER", "_"), 
                baseline_time=CONFIG["DELTA_BASELINE"], 
                compare_time=current_loop_val
            )
            delta_folder_name = f"{CONFIG['DELTA_BASELINE']}_to_{current_loop_val}"
            output_subfolder = os.path.join(mode_folder, delta_folder_name)
    
    else:
        # if current_loop_val is None (no timepoints/groups specified to iterate through),put results in "ALL" folder
        print("\nStarting analysis across entire dataset (no loops)")
        output_subfolder = os.path.join(mode_folder, "ALL")


    ##### INNER LOOP: Iterate through the targeted comparisons #########
    for target_val in comp_vals:
        comparison_folder = os.path.join(output_subfolder, f"{base_val}_vs_{target_val}")    # set subfolder for results based on comparison

        # Setup the 2-condition dataframe
        if MODE == "GROUP_COMPARISON":
            # Keep only baseline group and target group
            col_to_filter = CONFIG["CONDITION_COLUMN"]
            df_2_groups = df_filtered[df_filtered[col_to_filter].astype(str).isin([str(base_val), str(target_val)])]
            
            run_statistical_tests(
                df_subset=df_2_groups, 
                group_col=CONFIG["CONDITION_COLUMN"], 
                baseline_name=base_val, 
                compare_name=target_val, 
                output_dir=comparison_folder,
                is_paired=False, 
                is_delta=False
            )

        elif MODE == "LONGITUDINAL":
            # Keep only baseline time and target time
            col_to_filter = CONFIG["TIME_COLUMN"]
            df_2_times = df_filtered[df_filtered[col_to_filter].astype(str).isin([str(base_val), str(target_val)])]
            
            run_statistical_tests(
                df_subset=df_2_times, 
                group_col=CONFIG["TIME_COLUMN"], 
                baseline_name=base_val, 
                compare_name=target_val, 
                output_dir=comparison_folder,
                is_paired=True,     # Longitudinal assumes paired limma
                is_delta=False
            )

        elif MODE == "DELTA":
            # same as group coparison, but using globally calculated delta dataframe
            col_to_filter = CONFIG["CONDITION_COLUMN"]
            df_2_groups = df_filtered[df_filtered[col_to_filter].astype(str).isin([str(base_val), str(target_val)])]
            
            run_statistical_tests(
                df_subset=df_2_groups, 
                group_col=CONFIG["CONDITION_COLUMN"], 
                baseline_name=base_val, 
                compare_name=target_val, 
                output_dir=comparison_folder,
                is_paired=False,    # Deltas collapse pairs into independent change values
                is_delta=True
            )

print("\n-----------------------------------")


# =========================================================
# 3. Loops for running limma/logistic analyses
# =========================================================

if CONFIG.get("RUN_HINGE", False):
    
    # Import the hinge model function at the top of your script:
    # from hinge_model import run_hinge_model
    
    print("\n" + "="*60)
    print("ANALYSIS MODE: HINGE MODEL")
    print("="*60)

    # If the user left the groups listbox empty, run on the whole dataset as one group
    hinge_groups = CONFIG["HINGE_GROUPS_TO_RUN"] if CONFIG["HINGE_GROUPS_TO_RUN"] else [None]

    for group in hinge_groups:
        
        df_hinge_subset = df_main.copy()
        output_subfolder = os.path.join(CONFIG["OUTPUT_FOLDER"], "Hinge_Model")

        if group is not None:
            print(f"\n--- Running hinge model for group: {group} ---")
            df_hinge_subset = df_main[df_main[CONFIG["CONDITION_COLUMN"]].astype(str) == str(group)]
            output_subfolder = os.path.join(output_subfolder, str(group))
        else:
            print("\n--- Running hinge model on whole dataset (not split by group) ---")
            output_subfolder = os.path.join(output_subfolder, "ALL")

        os.makedirs(output_subfolder, exist_ok=True)

        
        run_hinge_model(
            df=df_hinge_subset,
            protein_cols=protein_cols,
            time_col=CONFIG["TIME_COLUMN"],
            patient_id_col=CONFIG["ID_COLUMN"],
            candidate_peaks=CONFIG["HINGE_PEAK_CANDIDATES"],
            num_plots=CONFIG["HINGE_NUM_PLOTS"],
            output_dir=output_subfolder
        )

print("\n-----------------------------------")
print("j'ai fini")