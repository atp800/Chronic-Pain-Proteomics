import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import pandas as pd
import os
import sys

# Import the global CONFIG dictionary
from config import CONFIG

class ProteomicsGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Proteomics Analysis Wizard")
        self.root.geometry("750x850")
        
        # ----- Variables -----
        self.v_input = tk.StringVar(value=CONFIG.get("INPUT_FILE_PATH", ""))
        self.v_sheet = tk.StringVar(value=CONFIG.get("SHEET_NAME", ""))
        self.v_output = tk.StringVar(value=CONFIG.get("OUTPUT_FOLDER", ""))
        
        self.v_id_col = tk.StringVar(value=CONFIG.get("ID_COLUMN", ""))
        self.v_group_col = tk.StringVar(value=CONFIG.get("CONDITION_COLUMN", ""))
        self.v_time_col = tk.StringVar(value=CONFIG.get("TIME_COLUMN", "None"))
        
        # Mode specific Variables
        self.v_group_base = tk.StringVar()
        self.v_long_base_time = tk.StringVar()
        self.v_delta_base_time = tk.StringVar()
        self.v_delta_comp_time = tk.StringVar()
        self.v_delta_base_grp = tk.StringVar()
        
        # Test Variables
        self.v_run_limma = tk.BooleanVar(value=CONFIG.get("RUN_LIMMA", True))
        self.v_run_logreg = tk.BooleanVar(value=CONFIG.get("RUN_LOGISTIC_REGRESSION", False))

        # ----- Build UI -----
        self.build_data_section()
        self.build_mode_tabs()
        self.build_tests_section()
        
        # Load initial columns if file exists
        if os.path.exists(self.v_input.get()):
            self.update_sheet_list()

    # ==========================================
    # UI BUILDERS
    # ==========================================
    def build_data_section(self):
        lf = ttk.LabelFrame(self.root, text="1. Data Definitions", padding=10)
        lf.pack(fill="x", padx=10, pady=5)
        
        # File Selectors
        ttk.Label(lf, text="Input File:").grid(row=0, column=0, sticky="w")
        ttk.Entry(lf, textvariable=self.v_input, width=50).grid(row=0, column=1, padx=5)
        ttk.Button(lf, text="Browse", command=self.browse_in).grid(row=0, column=2)

        ttk.Label(lf, text="Sheet:").grid(row=1, column=0, sticky="w", pady=5)
        self.cb_sheet = ttk.Combobox(lf, textvariable=self.v_sheet, width=47, state="readonly")
        self.cb_sheet.grid(row=1, column=1, padx=5, sticky="w")
        self.cb_sheet.bind("<<ComboboxSelected>>", self.load_columns)

        ttk.Label(lf, text="Output Folder:").grid(row=2, column=0, sticky="w")
        ttk.Entry(lf, textvariable=self.v_output, width=50).grid(row=2, column=1, padx=5)
        ttk.Button(lf, text="Browse", command=self.browse_out).grid(row=2, column=2)
        
        ttk.Separator(lf, orient='horizontal').grid(row=3, column=0, columnspan=3, sticky="ew", pady=10)

        # Column Definitions
        col_frame = ttk.Frame(lf)
        col_frame.grid(row=4, column=0, columnspan=3, sticky="ew")
        
        ttk.Label(col_frame, text="ID Column:").grid(row=0, column=0, sticky="w")
        self.cb_id = ttk.Combobox(col_frame, textvariable=self.v_id_col, state="readonly")
        self.cb_id.grid(row=0, column=1, padx=5, sticky="w")

        ttk.Label(col_frame, text="Group Column:").grid(row=1, column=0, sticky="w", pady=5)
        self.cb_group = ttk.Combobox(col_frame, textvariable=self.v_group_col, state="readonly")
        self.cb_group.grid(row=1, column=1, padx=5, sticky="w")
        self.cb_group.bind("<<ComboboxSelected>>", self.update_group_values)

        ttk.Label(col_frame, text="Time Column (used for delta):").grid(row=2, column=0, sticky="w")
        self.cb_time = ttk.Combobox(col_frame, textvariable=self.v_time_col, state="readonly")
        self.cb_time.grid(row=2, column=1, padx=5, sticky="w")
        self.cb_time.bind("<<ComboboxSelected>>", self.update_time_values)

        ttk.Label(col_frame, text="Columns to Ignore:").grid(row=0, column=2, sticky="nw", padx=(20, 5))
        self.lb_ignore = tk.Listbox(col_frame, selectmode="extended", height=4, exportselection=False)
        self.lb_ignore.grid(row=0, column=3, rowspan=3, sticky="ew")

    def build_mode_tabs(self):
        lf = ttk.LabelFrame(self.root, text="2. Analysis Mode Setup", padding=10)
        lf.pack(fill="both", expand=True, padx=10, pady=5)
        
        ttk.Label(lf, text="Select the tab for the analysis you want to run", 
                  font=("Arial", 9, "italic")).pack(anchor="w", pady=(0, 5))

        self.notebook = ttk.Notebook(lf)
        self.notebook.pack(fill="both", expand=True)

        # --- TAB 1: Group Comparison ---
        self.tab_group = ttk.Frame(self.notebook, padding=10)
        self.notebook.add(self.tab_group, text="Group Comparison")
        
        ttk.Label(self.tab_group, text="Baseline Group:").grid(row=0, column=0, sticky="w", pady=5)
        self.cb_grp_base = ttk.Combobox(self.tab_group, textvariable=self.v_group_base, state="readonly")
        self.cb_grp_base.grid(row=0, column=1, sticky="w", padx=5)

        ttk.Label(self.tab_group, text="Compare Against \n(Select multiple):").grid(row=1, column=0, sticky="nw", pady=5)
        self.lb_grp_comp = tk.Listbox(self.tab_group, selectmode="extended", height=4, exportselection=False)
        self.lb_grp_comp.grid(row=1, column=1, sticky="ew", padx=5, pady=5)

        ttk.Label(self.tab_group, text="Run for Timepoint(s) \n(Runs separately for each \nselection, or leave blank\nto run on all combined):").grid(row=2, column=0, sticky="nw", pady=5)
        self.lb_grp_times = tk.Listbox(self.tab_group, selectmode="extended", height=4, exportselection=False)
        self.lb_grp_times.grid(row=2, column=1, sticky="ew", padx=5, pady=5)

        # --- TAB 2: Longitudinal ---
        self.tab_long = ttk.Frame(self.notebook, padding=10)
        self.notebook.add(self.tab_long, text="Longitudinal Comparison")
        
        ttk.Label(self.tab_long, text="Baseline Time:").grid(row=0, column=0, sticky="w", pady=5)
        self.cb_long_base = ttk.Combobox(self.tab_long, textvariable=self.v_long_base_time, state="readonly")
        self.cb_long_base.grid(row=0, column=1, sticky="w", padx=5)

        ttk.Label(self.tab_long, text="Compare Against \n(Select multiple):").grid(row=1, column=0, sticky="nw", pady=5)
        self.lb_long_comp = tk.Listbox(self.tab_long, selectmode="extended", height=4, exportselection=False)
        self.lb_long_comp.grid(row=1, column=1, sticky="ew", padx=5, pady=5)

        ttk.Label(self.tab_long, text="Run for Group(s) \n(Runs separately for each \nselection, or leave blank\nto run on all combined):").grid(row=2, column=0, sticky="nw", pady=5)
        self.lb_long_grps = tk.Listbox(self.tab_long, selectmode="extended", height=4, exportselection=False)
        self.lb_long_grps.grid(row=2, column=1, sticky="ew", padx=5, pady=5)

        # --- TAB 3: Delta Comparison ---
        self.tab_delta = ttk.Frame(self.notebook, padding=10)
        self.notebook.add(self.tab_delta, text="Delta Comparison")
        
        ttk.Label(self.tab_delta, text="Delta Baseline Time:").grid(row=0, column=0, sticky="w", pady=5)
        self.cb_delta_base_time = ttk.Combobox(self.tab_delta, textvariable=self.v_delta_base_time, state="readonly")
        self.cb_delta_base_time.grid(row=0, column=1, sticky="w", padx=5)

        ttk.Label(self.tab_delta, text="Calculate Delta to:").grid(row=1, column=0, sticky="nw", pady=5)
        self.cb_delta_comp_time = ttk.Combobox(self.tab_delta, textvariable=self.v_delta_comp_time, state="readonly")
        self.cb_delta_comp_time.grid(row=1, column=1, sticky="w", padx=5)

        ttk.Label(self.tab_delta, text="Baseline Group:").grid(row=3, column=0, sticky="w", pady=2)
        self.cb_delta_base_grp = ttk.Combobox(self.tab_delta, textvariable=self.v_delta_base_grp, state="readonly")
        self.cb_delta_base_grp.grid(row=3, column=1, sticky="w", padx=5)

        ttk.Label(self.tab_delta, text="Compare Against \n(Select multiple):").grid(row=4, column=0, sticky="nw", pady=2)
        self.lb_delta_comp_grp = tk.Listbox(self.tab_delta, selectmode="extended", height=4, exportselection=False)
        self.lb_delta_comp_grp.grid(row=4, column=1, sticky="ew", padx=5, pady=2)

    def build_tests_section(self):
        lf = ttk.LabelFrame(self.root, text="3. Tests to Run", padding=10)
        lf.pack(fill="x", padx=10, pady=5)
        
        ttk.Checkbutton(lf, text="Limma", variable=self.v_run_limma).pack(side="left", padx=10)
        ttk.Checkbutton(lf, text="Logistic Regression", variable=self.v_run_logreg).pack(side="left", padx=10)
        
        ttk.Button(self.root, text="RUN ANALYSIS", command=self.on_submit).pack(pady=15, ipady=5, fill="x", padx=20)

    # ==========================================
    # LOGIC & UPDATES
    # ==========================================
    def browse_in(self):
        f = filedialog.askopenfilename(filetypes=[("Excel", "*.xlsx *.xls")])
        if f: 
            self.v_input.set(f)
            self.update_sheet_list()

    def browse_out(self):
        d = filedialog.askdirectory()
        if d: self.v_output.set(d)

    def update_sheet_list(self):
        path = self.v_input.get()
        if os.path.exists(path):
            try:
                xl = pd.ExcelFile(path)
                self.cb_sheet['values'] = xl.sheet_names
                if self.v_sheet.get() not in xl.sheet_names:
                    self.v_sheet.set(xl.sheet_names[0])
                self.load_columns()
            except Exception as e:
                print(f"Error reading sheets: {e}")

    def load_columns(self, event=None):
        path = self.v_input.get()
        sheet = self.v_sheet.get()
        if os.path.exists(path) and sheet:
            try:
                df = pd.read_excel(path, sheet_name=sheet, nrows=0)
                cols = list(df.columns)
                
                self.cb_id['values'] = cols
                self.cb_group['values'] = ["None"] + cols
                self.cb_time['values'] = ["None"] + cols
                
                self.lb_ignore.delete(0, tk.END)
                for col in cols:
                    self.lb_ignore.insert(tk.END, col)
                    
                self.update_group_values()
                self.update_time_values()
            except Exception as e:
                print(f"Error loading columns: {e}")

    def update_group_values(self, event=None):
        path = self.v_input.get()
        sheet = self.v_sheet.get()
        grp_col = self.v_group_col.get()
        
        self.cb_grp_base.set("")
        self.cb_delta_base_grp.set("")
        self.lb_grp_comp.delete(0, tk.END)
        self.lb_long_grps.delete(0, tk.END)
        self.lb_delta_comp_grp.delete(0, tk.END)
        
        if grp_col and grp_col != "None" and os.path.exists(path):
            try:
                df = pd.read_excel(path, sheet_name=sheet, usecols=[grp_col])
                uniques = sorted([str(x) for x in df[grp_col].dropna().unique()])
                
                self.cb_grp_base['values'] = uniques
                self.cb_delta_base_grp['values'] = uniques
                for u in uniques:
                    self.lb_grp_comp.insert(tk.END, u)
                    self.lb_long_grps.insert(tk.END, u)
                    self.lb_delta_comp_grp.insert(tk.END, u)
            except Exception:
                print("Error loading group values")

    def update_time_values(self, event=None):
        path = self.v_input.get()
        sheet = self.v_sheet.get()
        time_col = self.v_time_col.get()
        
        self.cb_long_base.set("")
        self.cb_delta_base_time.set("")
        self.cb_delta_comp_time.set("")
        self.lb_grp_times.delete(0, tk.END)
        self.lb_long_comp.delete(0, tk.END)
        
        # Enable/Disable tabs based on time column presence
        if not time_col or time_col == "None":
            self.notebook.tab(1, state="disabled")
            self.notebook.tab(2, state="disabled")
            self.lb_grp_times.config(state="disabled")
            return
        else:
            self.notebook.tab(1, state="normal")
            self.notebook.tab(2, state="normal")
            self.lb_grp_times.config(state="normal")
            
        if os.path.exists(path):
            try:
                df = pd.read_excel(path, sheet_name=sheet, usecols=[time_col])
                uniques = sorted([str(x) for x in df[time_col].dropna().unique()])
                
                self.cb_long_base['values'] = uniques
                self.cb_delta_base_time['values'] = uniques
                self.cb_delta_comp_time['values'] = uniques
                for u in uniques:
                    self.lb_grp_times.insert(tk.END, u)
                    self.lb_long_comp.insert(tk.END, u)
            except Exception:
                print("Error loading time values")

    def get_listbox_vals(self, listbox):
        return[listbox.get(i) for i in listbox.curselection()]

    def on_submit(self):
        # 1. Update Standard Global Configurations
        CONFIG["INPUT_FILE_PATH"] = self.v_input.get()
        CONFIG["SHEET_NAME"] = self.v_sheet.get()
        CONFIG["OUTPUT_FOLDER"] = self.v_output.get()
        
        CONFIG["ID_COLUMN"] = self.v_id_col.get()
        CONFIG["CONDITION_COLUMN"] = self.v_group_col.get() if self.v_group_col.get() != "None" else None
        CONFIG["TIME_COLUMN"] = self.v_time_col.get() if self.v_time_col.get() != "None" else None
        CONFIG["UNNEEDED_COLUMNS"] = self.get_listbox_vals(self.lb_ignore)
        
        CONFIG["RUN_LIMMA"] = self.v_run_limma.get()
        CONFIG["RUN_LOGISTIC_REGRESSION"] = self.v_run_logreg.get()

        # 2. Extract Mode-Specific Configuration
        selected_tab = self.notebook.index("current")
        
        if selected_tab == 0:
            CONFIG["ANALYSIS_MODE"] = "GROUP_COMPARISON"
            CONFIG["BASELINE_VAL"] = self.v_group_base.get()
            CONFIG["COMPARE_VALS"] = self.get_listbox_vals(self.lb_grp_comp)
            CONFIG["LOOP_VALS"] = self.get_listbox_vals(self.lb_grp_times) # Timepoints to loop through
            
        elif selected_tab == 1:
            CONFIG["ANALYSIS_MODE"] = "LONGITUDINAL"
            CONFIG["BASELINE_VAL"] = self.v_long_base_time.get()
            CONFIG["COMPARE_VALS"] = self.get_listbox_vals(self.lb_long_comp)
            CONFIG["LOOP_VALS"] = self.get_listbox_vals(self.lb_long_grps) # Groups to loop through
            
        elif selected_tab == 2:
            CONFIG["ANALYSIS_MODE"] = "DELTA"
            CONFIG["DELTA_BASELINE"] = self.v_delta_base_time.get()
            CONFIG["DELTA_COMPARISON"] = self.v_delta_comp_time.get()
            CONFIG["BASELINE_VAL"] = self.v_delta_base_grp.get()
            CONFIG["COMPARE_VALS"] = self.get_listbox_vals(self.lb_delta_comp_grp)
            CONFIG["LOOP_VALS"] = [None]    #self.get_listbox_vals(self.lb_delta_grps) # Groups to loop through

            # Validation for delta mode
            if not CONFIG["DELTA_BASELINE"] or not CONFIG["DELTA_COMPARISON"]:
                messagebox.showerror("Error", "Please select a Delta Baseline and Comparison time.")
                return


        # Validate minimum inputs
        if not CONFIG["BASELINE_VAL"] or not CONFIG["COMPARE_VALS"]:
            messagebox.showerror("Error", "Please select a baseline value and at least one comparison value.")
            return

        self.root.destroy()

def run_gui():
    root = tk.Tk()
    try:
        icon_img = tk.PhotoImage(file="logo.png") 
        root.iconphoto(True, icon_img)
    except:
        pass # Gracefully fail if logo missing
    app = ProteomicsGUI(root)
    root.mainloop()


if __name__ == "__main__":
    run_gui()