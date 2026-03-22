'''
Configuration and settings
'''

CONFIG = {
    "LAUNCH_GUI": True,
    
    # File Paths
    "INPUT_FILE_PATH": "Chronic Code/Cleaned_Data - Insoluble.xlsx",
    "SHEET_NAME": "Sheet 1",
    "OUTPUT_FOLDER": "Chronic Code/Analysis_Output",
    
    # Column Names
    "ID_COLUMN": "Sample_ID",
    "CONDITION_COLUMN": "Group",
    "DELTA_COLUMN": "Time",
    "UNNEEDED_COLUMNS": [],
    "FILTER_COLUMN": None,
    "FILTER_VALUE": [],
    "PROTEIN_COLUMNS": [],
    
    # Pre-Processing Settings
    "ALREADY_LOG2": True,
    "ALREADY_NORMALISED": True,

    # Other Settings
    "REPLACE_VALS_WITH_DELTAS": True,
    "LIMMA_IS_PAIRED": False,
    "MULTI_FILTER_TOGGLE": False,
    "ID_DELIMITER": "-",
    "P_THRESHOLD": 0.05,
    "LOG_FC_THRESHOLD": 0.58,
    
    # Tests to Run
    "RUN_PCA": False,
    "RUN_LIMMA": True,
    "RUN_SPEARMANS": False,
    "RUN_LASSO": False,
    "RUN_LOGISTIC_REGRESSION": False
}