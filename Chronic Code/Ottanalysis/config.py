'''
Configuration and settings
'''

CONFIG = {
    "LAUNCH_GUI": True,
    
    # File Paths
    "INPUT_FILE_PATH": "C:/Users/apana/OneDrive/Documents/Research/Uni of Western Sydney/Pain Project/Chronic Code/Cleaned Data - 23-3-26/Insoluble.xlsx",
    "SHEET_NAME": "Sheet 1",
    "OUTPUT_FOLDER": "Chronic Code/Ottanalysis/Analysis_Output/Insoluble",
    
    # Column Names
    "ID_COLUMN": "Sample_ID",
    "CONDITION_COLUMN": "Group",
    "TIME_COLUMN": "Time",
    "UNNEEDED_COLUMNS": [],
    # "FILTER_COLUMN": None,
    # "FILTER_VALUE": [],
    "PROTEIN_COLUMNS": [],
    
    # Pre-Processing Settings
    "ALREADY_LOG2": True,
    "ALREADY_NORMALISED": True,

    # Configuration Settings
    "ANALYSIS_MODE": "GROUP_COMPARISON",
    "BASELINE_VAL": None,
    "COMPARE_VALS": [],
    "LOOP_VALS": [],
    "DELTA_BASELINE": None,
    "HINGE_PEAK_CANDIDATES": [],
    "HINGE_GROUPS_TO_RUN": [],
    "HINGE_NUM_PLOTS": 20,
    # "DELTA_COMPARISON": None,
    # "REPLACE_VALS_WITH_DELTAS": True,
    # "LIMMA_IS_PAIRED": False,
    # "MULTI_FILTER_TOGGLE": False,

    # Other Settings
    "ID_DELIMITER": "-",
    "P_THRESHOLD": 0.05,
    "LOG_FC_THRESHOLD": 0.58,
    
    # Tests to Run
    "RUN_LIMMA": True,
    "RUN_HINGE": False,
    "RUN_PCA": False,
    "RUN_LOGISTIC_REGRESSION": False,
    "RUN_SPEARMANS": False
}