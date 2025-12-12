'''
This script intakes cleaned proteomics data for analysis.

It expects columns containing data to analyse (e.g. protein columns), 
at least one condition/group column, and a unique ID column.

If subject-wise analysis is needed (e.g. paired limma, comparing data 
from the same subject at different timepoints), the prefix of the ID column 
should be a subject identifier, separated by a delimiter (e.g. '-', '_', ' ', '.', '--').
For example "Subj01-Time0", "Subj01-Time1", for subject 1 at timepoints 0 and 1 respectively. 

Replicates and pre-processing should already be handled.
'''



##########################################################
# IMPORTS 
##########################################################
# System
import os
import sys
import glob
import datetime
import json
import hashlib
import warnings
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

# Data Handling
import numpy as np
import pandas as pd
import scipy.stats as stats
import openpyxl
from openpyxl.utils import get_column_letter

# Logistic Regression
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import roc_curve, auc, confusion_matrix, accuracy_score, f1_score
from sklearn.base import clone

# Visualisation
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


##################################################################
# SET LOCATION OF R INSTALLATION
##################################################################

# OPTION A: Manually to your R installation path here
r_home_path = None                                                      # Example: r_home_path = r"C:\Program Files\R\R-4.3.1"

# OPTION B: Auto-detect standard Windows R location if not set above
if not r_home_path:
    try:
        possible_paths = glob.glob(r"C:\Program Files\R\R-*")
        if possible_paths:
            # Pick the highest version number
            possible_paths.sort(reverse=True)
            r_home_path = possible_paths[0]
            print(f"Auto-detected R installation: {r_home_path}")
    except Exception:
        pass

if r_home_path and os.path.exists(r_home_path):
    # Set R_HOME environment variable
    os.environ['R_HOME'] = r_home_path
    
    # Set R_USER to current user to avoid permission errors
    os.environ['R_USER'] = os.path.expanduser("~")
    
    # Add R architecture-specific bin to PATH (Required for DLL loading)
    # Checks for x64 (64-bit) or i386 (32-bit)
    r_bin_x64 = os.path.join(r_home_path, 'bin', 'x64')
    r_bin_root = os.path.join(r_home_path, 'bin')
    
    if os.path.exists(r_bin_x64):
        os.environ['PATH'] = r_bin_x64 + ";" + os.environ['PATH']
    elif os.path.exists(r_bin_root):
        os.environ['PATH'] = r_bin_root + ";" + os.environ['PATH']
    else:
        print("Warning: Could not find R/bin folder. rpy2 might fail.")
else:
    print("Warning: Could not determine R_HOME. If rpy2 fails, set 'r_home_path' manually in the script.")


# R Imports
import rpy2
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.packages import importr


# ##################################################################
# # PARAMETERS 
# ##################################################################
LAUNCH_GUI = True

# File Paths
INPUT_FILE_PATH  = "Chronic Code/Cleaned_Data.xlsx"                         # NEED TO ADD BACK IN GROUP AND TIME COLUMNS (could get time from sample ID)
SHEET_NAME       = "Correlation_Matrix"
ORIGINAL_DATA_DF = pd.read_excel(INPUT_FILE_PATH, sheet_name=SHEET_NAME)
OUTPUT_FILE_PATH = "Chronic Code/Analysis_Output"

# Column Names
ID_COLUMN       = "Sample_ID" 
CONDITION_COLUMN = "Group"
DELTA_COLUMN = "Time"           # Column to calculate delta on if REPLACE_VALS_WITH_DELTAS is True
UNNEEDED_COLUMNS = ["Time"]    # CHECK SELECTED COLUMNS ARE ACTUALLY BEING REMOVED FROM ANALYSIS                                             # use for unneeded columns
PROTEIN_COLUMNS  = [col for col in ORIGINAL_DATA_DF.columns if (col not in ID_COLUMN and col not in UNNEEDED_COLUMNS and col != CONDITION_COLUMN and col != DELTA_COLUMN)]  

# Pre-Processing Settings
ALREADY_LOG2 = True             # NOT IMPLEMENTED
ALREADY_NORMALISED = True       # NOT IMPLEMENTED

# Tests to Run
RUN_PCA = False                 # NOT IMPLEMENTED
RUN_LIMMA = True
RUN_SPEARMANS = False           # NOT IMPLEMENTED
RUN_LASSO = False               # NOT IMPLEMENTED
RUN_LOGISTIC_REGRESSION = False

# Other Settings
ID_DELIMITER = "-"              # NOT IN GUI    # Delimiter to extract subject ID from sample ID for paired or subject-based analysis
P_THRESHOLD = 0.05
LOG_FC_THRESHOLD = 0         # 0.58 standard (1.5x), 0.38 relaxed (1.3x), 0.26 very relaxed (1.2x), 0 any change (p-value only)
LIMMA_IS_PAIRED = True          # NOT IN GUI    # Set to True for paired limma analysis (requires subject IDs in ID column)
                                # Recommend and give the option to switch if pairs exist/don't exist in condition column when limma is run
REPLACE_VALS_WITH_DELTAS = True         # If true makes and additional slector visible to choose which column to calculate delta for, Aand runs tests on delta values - turns limma into an interaction analysis
                                # e.g. get the difference from d0 to d14 for both groups and compare those differences (difference of differences)
                                # ALSO OVERWRITES LIMMA_IS_PAIRED: sets to false - pairs are now combined into single delta value


# # ADD A WAY TO CONTROL WHICH GROUP TO  USE FOR TIME COMPARISON
# ORIGINAL_DATA_DF = ORIGINAL_DATA_DF[ORIGINAL_DATA_DF['Group'] == 'Resp']
# ORIGINAL_DATA_DF.reset_index(drop=True, inplace=True)

# # ADD A WAY TO CONTROL WHICH TIME TO USE FOR GROUP COMPARISON
ORIGINAL_DATA_DF = ORIGINAL_DATA_DF[ORIGINAL_DATA_DF['Time'] == 'D14']
ORIGINAL_DATA_DF.reset_index(drop=True, inplace=True)



# add settings for groupwise missing value filtering threshold etc. 
# - make sure toggles added to gui
# - split toggles by tests they apply to??

# add ID_DELIMETER setting to gui
# add LIMMA_IS_PAIRED setting to gui

# Check ID coliumn is unique. 
# If not warn user and append condition with a hyphen.
# Check again is unqiue and if not warn again append a number

#####
# ASK CHATGPT TO CHECK FOR SETTINGS/MAGIC NUMBERS IN CODE AND ADD ALL AS PARAMETERS HERE AND IN GUI
#####


##################################################################
# GUI CONFIGURATION
##################################################################
def run_gui_selector():
    # Use global keywords to update the actual variables at the top of the script
    global INPUT_FILE_PATH, SHEET_NAME, OUTPUT_FILE_PATH
    global ID_COLUMN, CONDITION_COLUMN, UNNEEDED_COLUMNS, DELTA_COLUMN
    global ALREADY_LOG2, ALREADY_NORMALISED, REPLACE_VALS_WITH_DELTAS
    global RUN_PCA, RUN_LIMMA, RUN_SPEARMANS, RUN_LASSO, RUN_LOGISTIC_REGRESSION

    root = tk.Tk()
    root.title("Proteomics Analysis")
    root.geometry("650x850")
    icon_img = tk.PhotoImage(file="logo.png")   # otters are important and this should be considered when running the script
    root.iconphoto(True, icon_img)

    # --- Variables Linked to GUI Widgets ---
    v_input_path = tk.StringVar(value=INPUT_FILE_PATH)
    v_sheet = tk.StringVar(value=SHEET_NAME)
    v_output_path = tk.StringVar(value=OUTPUT_FILE_PATH)
    v_condition = tk.StringVar(value=CONDITION_COLUMN)
    v_id_col = tk.StringVar(value=ID_COLUMN)
    v_delta_col = tk.StringVar(value=DELTA_COLUMN)
    
    v_log2 = tk.BooleanVar(value=ALREADY_LOG2)
    v_norm = tk.BooleanVar(value=ALREADY_NORMALISED)
    v_pca = tk.BooleanVar(value=RUN_PCA)
    v_limma = tk.BooleanVar(value=RUN_LIMMA)
    v_spearman = tk.BooleanVar(value=RUN_SPEARMANS)
    v_lasso = tk.BooleanVar(value=RUN_LASSO)
    v_logistic = tk.BooleanVar(value=RUN_LOGISTIC_REGRESSION)
    v_calc_deltas = tk.BooleanVar(value=REPLACE_VALS_WITH_DELTAS)

    # Ensure this defaults to a standard boolean (False if None) to avoid third state (line through box)
    initial_delta = bool(REPLACE_VALS_WITH_DELTAS)
    v_calc_deltas = tk.BooleanVar(value=initial_delta)

    # Toggle function to active delta column selector
    def toggle_delta_selector(*args):
        if v_calc_deltas.get():
            cb_delta.config(state="readonly") # Enable dropdown
        else:
            cb_delta.config(state="disabled") # Disable dropdown
    # Attach the listener to variable to run every time it changes
    v_calc_deltas.trace_add("write", toggle_delta_selector)

    # --- Helper: Update Column Lists based on Sheet ---
    def load_columns(event=None):
        path = v_input_path.get()
        sheet = v_sheet.get()
        if os.path.exists(path) and sheet:
            try:
                # Read header only
                df = pd.read_excel(path, sheet_name=sheet, nrows=0)
                cols = list(df.columns)
                
                # Put columns from excel sheet into selectors
                cb_id['values'] = cols 
                cb_condition['values'] = cols
                cb_delta['values'] = cols 
                
                lb_unneeded.delete(0, tk.END)
                
                for i, col in enumerate(cols):
                    lb_unneeded.insert(tk.END, col)
                    
                    # Pre-select defaults if they match
                    if UNNEEDED_COLUMNS and col in UNNEEDED_COLUMNS:
                        lb_unneeded.selection_set(i)
                
                # Set Defaults
                if ID_COLUMN in cols: v_id_col.set(ID_COLUMN)
                if CONDITION_COLUMN in cols: v_condition.set(CONDITION_COLUMN)
                if DELTA_COLUMN in cols: v_delta_col.set(DELTA_COLUMN)

            except Exception as e:
                print(f"Error loading columns: {e}")

    def update_sheet_list(event=None):
        path = v_input_path.get()
        if os.path.exists(path):
            try:
                xl = pd.ExcelFile(path)
                cb_sheet['values'] = xl.sheet_names
                if SHEET_NAME in xl.sheet_names:
                    v_sheet.set(SHEET_NAME)
                else:
                    v_sheet.set(xl.sheet_names[0])
                load_columns() # Load columns for the default sheet
            except:
                pass

    def browse_in():
        f = filedialog.askopenfilename(filetypes=[("Excel", "*.xlsx *.xls")])
        if f: 
            v_input_path.set(f)
            update_sheet_list()

    def browse_out():
        d = filedialog.askdirectory()
        if d: v_output_path.set(d)

    def on_submit():
        # Update global variables
        global INPUT_FILE_PATH, SHEET_NAME, OUTPUT_FILE_PATH, CONDITION_COLUMN
        global ID_COLUMN, UNNEEDED_COLUMNS, DELTA_COLUMN
        global ALREADY_LOG2, ALREADY_NORMALISED, REPLACE_VALS_WITH_DELTAS
        global RUN_PCA, RUN_LIMMA, RUN_SPEARMANS, RUN_LASSO

        INPUT_FILE_PATH = v_input_path.get()
        SHEET_NAME = v_sheet.get()
        OUTPUT_FILE_PATH = v_output_path.get()
        
        # Get dropdown values
        CONDITION_COLUMN = v_condition.get()
        ID_COLUMN = v_id_col.get()
        DELTA_COLUMN = v_delta_col.get()
        
        # Get listbox selections
        all_options = lb_unneeded.get(0, tk.END)
        UNNEEDED_COLUMNS = [all_options[i] for i in lb_unneeded.curselection()]

        # Get toggle values
        ALREADY_LOG2 = v_log2.get()
        ALREADY_NORMALISED = v_norm.get()
        REPLACE_VALS_WITH_DELTAS = v_calc_deltas.get()

        RUN_PCA = v_pca.get()
        RUN_LIMMA = v_limma.get()
        RUN_SPEARMANS = v_spearman.get()
        RUN_LASSO = v_lasso.get()
        RUN_LOGISTIC_REGRESSION = v_logistic.get()

        root.destroy()

    # --- Layout ---
    pad = 5
    
    # 1. Files Frame
    lf_files = ttk.LabelFrame(root, text="Files & Data Source", padding=pad)
    lf_files.pack(fill="x", padx=pad, pady=pad)
    
    ttk.Label(lf_files, text="Input File:").grid(row=0, column=0, sticky="w")
    ttk.Entry(lf_files, textvariable=v_input_path, width=50).grid(row=0, column=1)
    ttk.Button(lf_files, text="Browse", command=browse_in).grid(row=0, column=2)

    ttk.Label(lf_files, text="Sheet:").grid(row=1, column=0, sticky="w")
    cb_sheet = ttk.Combobox(lf_files, textvariable=v_sheet, width=47)
    cb_sheet.grid(row=1, column=1)
    cb_sheet.bind("<<ComboboxSelected>>", load_columns)

    ttk.Label(lf_files, text="Output Folder:").grid(row=2, column=0, sticky="w")
    ttk.Entry(lf_files, textvariable=v_output_path, width=50).grid(row=2, column=1)
    ttk.Button(lf_files, text="Browse", command=browse_out).grid(row=2, column=2)


    # 2. Columns Frame
    lf_cols = ttk.LabelFrame(root, text="Column Definitions", padding=pad)
    lf_cols.pack(fill="both", expand=True, padx=pad, pady=pad)
    
    # Configure grid: Col 0 for Dropdowns, Col 1 for Listbox
    lf_cols.columnconfigure(0, weight=1)
    lf_cols.columnconfigure(1, weight=1)

    # LEFT SIDE    
    # ID Column (single-select dropdown)
    ttk.Label(lf_cols, text="Sample ID Column:").grid(row=0, column=0, sticky="w")
    cb_id = ttk.Combobox(lf_cols, textvariable=v_id_col, width=30)
    cb_id.grid(row=1, column=0, sticky="w", padx=pad, pady=(0, 15)) # Add spacing below

    # Condition Column (single-select dropdown)
    ttk.Label(lf_cols, text="Group/Condition Column:").grid(row=2, column=0, sticky="w")
    cb_condition = ttk.Combobox(lf_cols, textvariable=v_condition, width=30)
    cb_condition.grid(row=3, column=0, sticky="w", padx=pad)

    # Delta Column (single-select dropdown)
    ttk.Label(lf_cols, text="Delta Column (for limma interaction analysis):").grid(row=4, column=0, sticky="w")
    cb_delta = ttk.Combobox(lf_cols, textvariable=v_delta_col, width=30, state="disabled") # Start disabled
    cb_delta.grid(row=5, column=0, sticky="w", padx=pad, pady=(0, 5))

    # RIGHT SIDE    
    # Unneeded Columns (multi-select)
    ttk.Label(lf_cols, text="Unneeded/Metadata Columns (can select multiple):").grid(row=0, column=1, sticky="w")
    
    # Wrapper Frame for Listbox + Scrollbar
    frame_lb_un = ttk.Frame(lf_cols)
    # Span multiple rows so it stands tall next to the two dropdowns on the left
    frame_lb_un.grid(row=1, column=1, rowspan=4, sticky="nsew", padx=pad)
    
    sb_un = ttk.Scrollbar(frame_lb_un)
    sb_un.pack(side="right", fill="y")
    
    lb_unneeded = tk.Listbox(frame_lb_un, selectmode="multiple", height=8, exportselection=False, yscrollcommand=sb_un.set)
    lb_unneeded.pack(side="left", fill="both", expand=True)
    
    sb_un.config(command=lb_unneeded.yview)


    # 3. Settings Frame
    lf_sets = ttk.LabelFrame(root, text="Analysis Settings", padding=pad)
    lf_sets.pack(fill="x", padx=pad, pady=pad)

    ttk.Checkbutton(lf_sets, text="Data Already Log2 Transformed", variable=v_log2).grid(row=0, column=0, sticky="w")
    ttk.Checkbutton(lf_sets, text="Data Already Normalised", variable=v_norm).grid(row=0, column=1, sticky="w")
    ttk.Checkbutton(lf_sets, text="Replace Values With Deltas (e.g. Interaction Analysis)", variable=v_calc_deltas).grid(row=1, column=0, columnspan=2, sticky="w", pady=5)

    ttk.Separator(lf_sets, orient='horizontal').grid(row=2, column=0, columnspan=2, sticky="ew", pady=5)

    ttk.Checkbutton(lf_sets, text="Run PCA", variable=v_pca).grid(row=3, column=0, sticky="w")
    ttk.Checkbutton(lf_sets, text="Run Limma", variable=v_limma).grid(row=3, column=1, sticky="w")
    ttk.Checkbutton(lf_sets, text="Run Spearmans", variable=v_spearman).grid(row=4, column=0, sticky="w")
    ttk.Checkbutton(lf_sets, text="Run LASSO", variable=v_lasso).grid(row=4, column=1, sticky="w")
    ttk.Checkbutton(lf_sets, text="Run Logistic Regression", variable=v_logistic).grid(row=5, column=0, sticky="w")

    # Run Button
    ttk.Button(root, text="RUN ANALYSIS", command=on_submit).pack(pady=10, ipady=5, fill="x", padx=20)

    # CLose Button
    def on_closing():
        root.destroy()
        sys.exit("Analysis cancelled")
    root.protocol("WM_DELETE_WINDOW", on_closing)

    # Initial Load
    if os.path.exists(INPUT_FILE_PATH):
        update_sheet_list()
    toggle_delta_selector() # Run once to set initial state
    
    root.mainloop()

# Launch GUI if requested
if LAUNCH_GUI:
    run_gui_selector()
    print("Parameters updated via GUI")
else:
    print("Using hardcoded default parameters")

# Ensure Unneeded is a list if it wasn't set
if UNNEEDED_COLUMNS is None: UNNEEDED_COLUMNS = []



##################################################################
# START ANALYSIS
##################################################################

print(f"Loading: {INPUT_FILE_PATH} | Sheet: {SHEET_NAME}")

print("-- Parameters --")
print(f" Output Folder: {OUTPUT_FILE_PATH}")
print(f" ID Column: {ID_COLUMN}")
print(f" Condition Column: {CONDITION_COLUMN}")
print(f" Unneeded Columns: {UNNEEDED_COLUMNS}")

print("Running tests:")
if RUN_PCA: print(" - PCA")
if RUN_LIMMA: print(" - Limma")
if RUN_SPEARMANS: print(" - Spearmans")
if RUN_LASSO: print(" - LASSO")
if RUN_LOGISTIC_REGRESSION: print(" - LOGSITC REGRESSION")

try:
    ORIGINAL_DATA_DF = pd.read_excel(INPUT_FILE_PATH, sheet_name=SHEET_NAME)
    
    # Calculate protein columns based on others selected
    excluded_cols = {ID_COLUMN, CONDITION_COLUMN, DELTA_COLUMN} | set(UNNEEDED_COLUMNS)
    PROTEIN_COLUMNS  = [col for col in ORIGINAL_DATA_DF.columns if col not in excluded_cols]
    
    print(f"Condition: {CONDITION_COLUMN}")
    print(f"Proteins Identified: {len(PROTEIN_COLUMNS)}")

except Exception as e:
    print(f"Error loading data: {e}")
    sys.exit(1)

# Ensure output directory exists
if not os.path.exists(OUTPUT_FILE_PATH):
    os.makedirs(OUTPUT_FILE_PATH)



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
    # Create a copy to avoid modifying original
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
    of samples in AT LEAST ONE experimental group.
    
    Preserves on/off biomarkers (e.g., present in disease, absent in control).
    """
    print(f"Applying Group-Based Missingness Filter (Threshold: {threshold:.0%})...")
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
def derive_subject_ids(sample_id_list, delimiter=ID_DELIMITER):
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
# SWITCH TO DELTA MODE
##################################################################
'''
If selected in GUI, ENTIRE dataset is converted to difference scores based on selected DELTA_COLUMN
All subsequent tests (Limma, LogReg etc.) will run on the change values
'''

if REPLACE_VALS_WITH_DELTAS:
    print("\n" + "="*65)
    print(f"TRANSFORMING DATA: Calculating Deltas based on '{DELTA_COLUMN}'")
    print("="*65)


    # Check selected delta column exists
    if DELTA_COLUMN not in ORIGINAL_DATA_DF.columns:
        print(f"ERROR: Selected Delta Column '{DELTA_COLUMN}' not found in data.")
        sys.exit(1)

    # Get subject IDs to caluclate change for
    print(f"Deriving subjects using delimiter '{ID_DELIMITER}'...")
    subjects = derive_subject_ids(ORIGINAL_DATA_DF[ID_COLUMN].tolist(), delimiter=ID_DELIMITER)
    
    df_work = ORIGINAL_DATA_DF.copy()
    df_work = impute_missing_vals(df_work)
    df_work['Temp_Subject_ID'] = subjects

    # Identify timepoints - NEED MORE ROBUST METHOD HERE - NEED TO HANDLE IF THERE ARE MORE THAN TWO TIME POINTS OR THEY GET SORTE INTOT HE WRONG ORDER
    # Sort them to assume Time 2 - Time 1 (e.g., D14 - D0)
    timepoints = sorted(df_work[DELTA_COLUMN].unique())
    
    if len(timepoints) != 2:
        print(f"CRITICAL ERROR: Delta calculation requires exactly 2 timepoints.")
        print(f"Found {len(timepoints)} in column '{DELTA_COLUMN}': {timepoints}")
        sys.exit(1)

    t_start = timepoints[0]
    t_end   = timepoints[1]
    print(f"Calculating Change: {t_end} (Final) - {t_start} (Initial)")

    # Calculate delta
    try:
        # Check for duplicates
        if df_work.duplicated(subset=['Temp_Subject_ID', DELTA_COLUMN]).any():                          
            print("Error: Duplicate samples found (Same Subject + Same Timepoint). Cannot pivot.")
            sys.exit(1)

        df_pivot = df_work.pivot(index='Temp_Subject_ID', columns=DELTA_COLUMN, values=PROTEIN_COLUMNS) # Pivot: Index=Subject, Columns=Time, Values=Proteins
        vals_end = df_pivot.xs(t_end, axis=1, level=DELTA_COLUMN)
        vals_start = df_pivot.xs(t_start, axis=1, level=DELTA_COLUMN)
        df_delta = vals_end - vals_start
        
        n_original = len(df_pivot)
        n_final = len(df_delta)
        
        print(f"Subjects with complete pairs: {n_final} (Dropped {n_original - n_final})")

        # Restore group data (condition column) and ID column - assumes group doesn't change over time
        subj_group_map = df_work.drop_duplicates('Temp_Subject_ID').set_index('Temp_Subject_ID')[CONDITION_COLUMN]
        df_delta[CONDITION_COLUMN] = df_delta.index.map(subj_group_map)
        df_delta[ID_COLUMN] = df_delta.index
        df_delta.reset_index(drop=True, inplace=True)
        
        # Replace the main dataframe with delta dataframe
        ORIGINAL_DATA_DF = df_delta
        
        # Ensure limma runs as independent, since pairs have been combiend into single delta value
        LIMMA_IS_PAIRED = False 
        
        print("SUCCESS: dataset values replaced with delta values")
        print("All selected tests (Limma, Logistic Regression, etc.) will run on these deltas")
        print("-" * 65 + "\n")

    except Exception as e:
        print(f"Error during delta calculation: {e}")
        sys.exit(1)


##################################################################
# TEST: LIMMA DIFFERENTIAL EXPRESSION
##################################################################
# Limma as paired vs independent toggle in parameters section at top of script
if RUN_LIMMA:
    print("\n------------------------------------------------------------------")
    print("Running Limma...")
    if LIMMA_IS_PAIRED:
        print("Mode: PAIRED Analysis (Grouping by Subject ID)") # COMPARES ALPHABETICALLY FIRST TWO GROUPS - will need to add a gui or input option for selecting if there are more
    else:
        print("Mode: INDEPENDENT Analysis") # CHECK IDS ARE UNIQUE
    print("------------------------------------------------------------------")

    # R-Python DataFrame Conversion    
    pandas2ri.activate()
    limma = importr('limma')
    base = importr('base')
    stats = importr('stats')

    # 1. PREPARE DATA
    X_limma = ORIGINAL_DATA_DF[PROTEIN_COLUMNS]
    print("Imputing missing values...")                         # APPLY CUSTOM MISSING VALUE FILTERING first???
    if X_limma.isna().sum().sum() > 0:
        X_limma = impute_missing_vals(X_limma, shift=1.8, width=0.3)

    # Need to transpose: Rows = Proteins, Columns = Samples
    r_expression_matrix = X_limma.T

    # Prepare Groups
    clean_groups = [str(g).replace(" ", "_").replace("-", "_") for g in ORIGINAL_DATA_DF[CONDITION_COLUMN]]
    r_groups = ro.StrVector(clean_groups)
    
    # Prepare IDs
    raw_sample_ids = ORIGINAL_DATA_DF[ID_COLUMN].tolist()
     
    if LIMMA_IS_PAIRED:
        derived_subjects = derive_subject_ids(raw_sample_ids, delimiter=ID_DELIMITER)
        r_subjects = ro.StrVector(derived_subjects)
    else:
        r_subjects = ro.StrVector(raw_sample_ids)


    # 2. CREATE DESIGN MATRIX
    ro.globalenv['expression_data'] = r_expression_matrix
    ro.globalenv['groups'] = r_groups
    ro.globalenv['subjects'] = r_subjects
    ro.r('groups <- as.factor(groups)')                                             # force R to treat groups as a factor

    if LIMMA_IS_PAIRED:
        print("Creating Paired Design Matrix (~ 0 + groups + subjects)...")
        
        if len(set(derived_subjects)) == len(derived_subjects):                     # Validation: Check if pairing is actually possible
             print("CRITICAL WARNING: Derived Subject IDs are all unique.")
             print("Limma cannot pair samples if every ID is different.")
             ro.r('design <- model.matrix(~ 0 + groups)')
        else:
            ro.r('subjects <- as.factor(subjects)')                                 # force R to treat subjects as a factor
            ro.r('design <- model.matrix(~ 0 + groups + subjects)')
    else:
        print("Creating Independent Design Matrix...")
        ro.r('design <- model.matrix(~ 0 + groups)')
        
    # Get column names directly from the R object 'design'
    col_names = list(ro.r('colnames(design)'))
    unique_groups_clean = sorted(list(set(clean_groups)))
    new_col_names = []

    # Clean the column names
    for c in col_names:
        for g in unique_groups_clean:
            # Check if group name is in the column name (and ignore subject columns)
            if g in c and "subjects" not in c:
                new_col_names.append(g)
                break
        else:
            new_col_names.append(c)
            
    # Assign new names back to the R object 'design'
    ro.globalenv['new_cols'] = ro.StrVector(new_col_names)
    ro.r('colnames(design) <- new_cols')
    design_matrix = ro.globalenv['design']                      # Pull the design matrix object for use in lmFit


    # 3. FIT LINEAR MODEL & CONTRASTS
    print("Fitting Linear Model...")
    ro.r('expression_data <- as.matrix(expression_data)')
    ro.r('fit <- lmFit(expression_data, design)')
    
    # Define Contrast (e.g., "Resp - NonResp")
    if len(unique_groups_clean) >= 2:
        compare_string = f"{unique_groups_clean[1]} - {unique_groups_clean[0]}"
        print(f"Defining Contrast: {compare_string}")
        
        # Run makeContrasts, staying in R
        ro.r(f'contrast.matrix <- makeContrasts("{compare_string}", levels=design)')
        ro.r('fit2 <- contrasts.fit(fit, contrast.matrix)')
        ro.r('fit2 <- eBayes(fit2)')
        
        # Extract results back into python
        r_results = ro.r('topTable(fit2, number=Inf, sort.by="P")')        # number=Inf gets all proteins
        
        # Convert back to Pandas
        with localconverter(ro.default_converter + pandas2ri.converter):
            df_results = ro.conversion.rpy2py(r_results)
        
        # Add the Protein Names (which are the index of the R result)
        df_results.index.name = "Protein"
        df_results.reset_index(inplace=True)
        
        # Save
        limma_out_path = os.path.join(OUTPUT_FILE_PATH, "Limma_Results.csv")
        df_results.to_csv(limma_out_path, index=False)
        print(f"Saved Limma results to: {limma_out_path}")
        
        # Identify Significant Proteins
        sig_proteins = df_results[(df_results['adj.P.Val'] < P_THRESHOLD) & (abs(df_results['logFC']) > LOG_FC_THRESHOLD)]
        print(f"Significant Proteins (adj.P < {P_THRESHOLD} & |logFC| > {LOG_FC_THRESHOLD}): {len(sig_proteins)}")

        # 4. VOLCANO PLOT
        pdf_path = os.path.join(OUTPUT_FILE_PATH, "Limma_Volcano_Plot.pdf")
        with PdfPages(pdf_path) as pdf:
            plt.figure(figsize=(10, 8))
            
            # Scatter plot
            sns.scatterplot(data=df_results, x='logFC', y='adj.P.Val', alpha=0.5, edgecolor=None)
            
            # Transform y-axis to -log10
            plt.yscale('log')
            plt.gca().invert_yaxis() # Small p-values at top
            
            # Add thresholds
            plt.axhline(0.05, color='r', linestyle='--', alpha=0.5, label='p=0.05')
            plt.axvline(0.58, color='b', linestyle='--', alpha=0.5, label='1.5x Fold Change')
            plt.axvline(-0.58, color='b', linestyle='--', alpha=0.5)
            
            plt.title(f"Volcano Plot: {compare_string}")
            plt.xlabel("Log2 Fold Change")
            plt.ylabel("Adjusted P-Value (Log Scale)")
            plt.legend()
            
            pdf.savefig()
            plt.close()
        print(f"Saved Volcano Plot to: {pdf_path}")

    else:
        print("Error: Need at least 2 groups to perform Limma comparison.")







##################################################################
# TEST: LOGISTIC REGRESSION
##################################################################

if RUN_LOGISTIC_REGRESSION:
    print("\n------------------------------------------------------------------")
    print("Running Logistic Regression (Avec L1/L2 Regularisation)...")
    print("------------------------------------------------------------------")
    
    # 1. PREPARE DATA
    X = ORIGINAL_DATA_DF[PROTEIN_COLUMNS]
    y = ORIGINAL_DATA_DF[CONDITION_COLUMN]

    # More aggressive group-based missing protein filtering
    print(f"Original dimensionality: {X.shape[1]} proteins")
    X = groupwise_missing_filter(X, y, threshold=0.7)               # Keep proteins with <30% missing in at least one group
    print(f"Proteins after dropping >30% missing: {X.shape[1]}")

    # Remove low-variance features
    X_temp_for_var = X.fillna(0)                                    # Temporarily fill NaNs for variance calculation
    selector = VarianceThreshold(threshold=0.1) 
    try:
        selector.fit(X_temp_for_var)
        cols_kept = X.columns[selector.get_support()]
        X = X[cols_kept]
        print(f"Proteins after dropping low variance: {X.shape[1]}")
    except ValueError:
        print("Warning: Variance filter failed (variance too low globally?), skipping.")

    # Check for NaNs and impute
    if X.isna().sum().sum() > 0:
        print("Missing values detected: imputing with down-shifted normal distribution")
        X = impute_missing_vals(X)

    # Encode Target
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)
    classes = le.classes_
    print(f"Classes: {classes}")
    
    if len(classes) != 2:
        print(f"Warning: {len(classes)} groups detected. This script is designed for 2 groups.")

    # 2. SCALING
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # 3. TRAIN MODEL & EVALUATE (UNBIASED)
    print("Training Logistic Regression (~10 mins to run)...")
    
    # Define the model structure (Nested CV: The model itself does CV to find C)
    clf = LogisticRegressionCV(
        Cs=10, 
        cv=5, 
        penalty='elasticnet', 
        solver='saga', 
        l1_ratios=[0.1, 0.5, 0.7, 0.9, 0.95, 1],
        random_state=42, 
        max_iter=10000, 
        scoring='accuracy',
        verbose=0,                  # increase for more print statements while running
        n_jobs=-1
    )

    print("Performing leave-one-out corss-validation...")
    clf.fit(X_scaled, y_encoded)
    cv_outer = LeaveOneOut()            # Simulates a test set for every single sample.
    
    # Lists to store results
    y_pred_unbiased = []
    y_probs_unbiased = []
    
    total_samples = len(X_scaled)
    
    # Loop through every sample
    for i, (train_index, test_index) in enumerate(cv_outer.split(X_scaled)):
        
        # --- PROGRESS BAR ---
        percent = (i + 1) / total_samples
        bar_length = 40
        filled_length = int(bar_length * percent)
        bar = '█' * filled_length + '-' * (bar_length - filled_length)
        sys.stdout.write(f'\rProgress: |{bar}| {int(percent * 100)}% ({i+1}/{total_samples})')      # \r returns cursor to start of line, allowing overwrite
        sys.stdout.flush()
        # --------------------------

        # Split Data
        X_train, X_test = X_scaled[train_index], X_scaled[test_index]
        y_train, y_test = y_encoded[train_index], y_encoded[test_index]
        
        # Clone the model (start fresh every time)
        model_fold = clone(clf) 
        model_fold.fit(X_train, y_train)
        
        # Predict
        prediction = model_fold.predict(X_test)[0]
        y_pred_unbiased.append(prediction)
        
        if len(classes) == 2:
            prob = model_fold.predict_proba(X_test)[0][1]
            y_probs_unbiased.append(prob)

    print("\nValidation Complete.\n")
    
    # Convert lists to numpy arrays for the metrics calculations
    y_pred_unbiased = np.array(y_pred_unbiased)
    if len(classes) == 2:
        y_probs_unbiased = np.array(y_probs_unbiased)

    # 4. CALCULATE METRICS (On Unbiased Predictions)
    # -----------------------------------------------------------
    acc = accuracy_score(y_encoded, y_pred_unbiased)
    
    if len(classes) == 2:
        f1 = f1_score(y_encoded, y_pred_unbiased)
    else:
        f1 = f1_score(y_encoded, y_pred_unbiased, average='weighted')
        
    cm = confusion_matrix(y_encoded, y_pred_unbiased)

    print(f"\n--- Model Performance (Nested CV Estimate) ---")
    print(f"Accuracy: {acc:.4f}")
    print(f"F1 Score: {f1:.4f}")
    print(f"Confusion Matrix (Rows=Actual, Cols=Predicted):\n{cm}")
    print("----------------------------------------------\n")

    # 5. EXTRACT FEATURES (From the model trained on ALL data)
    # -----------------------------------------------------------
    print(f"Best Regularization Strength (C): {clf.C_[0]}")
    print(f"Best L1 Ratio: {clf.l1_ratio_[0]}")
    
    if len(classes) == 2:
        coefs = clf.coef_.flatten()
    else:
        coefs = np.max(np.abs(clf.coef_), axis=0)

    mask = coefs != 0
    current_protein_names = X.columns
    selected_proteins = np.array(current_protein_names)[mask]
    selected_coefs = coefs[mask]
    
    print(f"Selected {len(selected_proteins)} proteins.")

    logistic_df = pd.DataFrame({
        'Protein': selected_proteins,
        'Coefficient': selected_coefs,
        'Abs_Coefficient': np.abs(selected_coefs)
    }).sort_values(by='Abs_Coefficient', ascending=False)
    
    logistic_csv_path = os.path.join(OUTPUT_FILE_PATH, "Logistic_Selected_Features.csv")
    logistic_df.to_csv(logistic_csv_path, index=False)
    print(f"Saved features to: {logistic_csv_path}")


    # 6. VISUALISATION
    # -----------------------------------------------------------
    pdf_path = os.path.join(OUTPUT_FILE_PATH, "Logistic_Performance_Plots.pdf")
    with PdfPages(pdf_path) as pdf:
        
        # Plot A: ROC Curve (Using Unbiased Probabilities)
        if len(classes) == 2:
            fpr, tpr, _ = roc_curve(y_encoded, y_probs_unbiased)
            roc_auc = auc(fpr, tpr)
            
            plt.figure(figsize=(8, 6))
            plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (CV AUC = {roc_auc:.2f})')
            plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title(f'ROC Curve (Nested CV): {classes[1]} vs {classes[0]}')
            plt.legend(loc="lower right")
            plt.grid(True, alpha=0.3)
            pdf.savefig()
            plt.close()

        # Plot B: Top Coefficients (Fixed the Seaborn Warning)
        if not logistic_df.empty:
            top_n = 20
            plot_df = logistic_df.head(top_n)
            
            plt.figure(figsize=(10, 8))
            # FIX: Added hue and legend=False
            sns.barplot(data=plot_df, x='Coefficient', y='Protein', hue='Protein', palette='vlag', legend=False)
            plt.title(f'Top {top_n} Logistic Coefficients (Model fitted on all data)')
            plt.xlabel('Coefficient Value (Log Odds)')
            plt.axvline(0, color='k', linewidth=1)
            plt.grid(axis='x', alpha=0.3)
            plt.tight_layout()
            pdf.savefig()
            plt.close()
        
    print(f"Saved plots to: {pdf_path}")