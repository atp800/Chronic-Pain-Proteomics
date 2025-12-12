# Pain Proteomics Analysis

**Project Summary**

- **Purpose:** Identify proteins associated with pain either by direct correlation with self-reported pain scores (VAS) or by detecting rise/fall patterns over time after a pain-inducing injection.

## Scripts

- **`Preprocessing.py`:** Ingests raw Excel sheets and produces `Combined_Data.xlsx`. Run this first.
- **`VAS_Correlation_Analysis.py`:** Computes Pearson correlations between protein levels and VAS pain scores for individual participants and across all participants. Exports `Significant_VAS_Correlations.xlsx` and shows scatter/regression plots for top hits.
- **`Time_Analysis.py`:** Fits quadratic curves across timepoints to detect bell/inverted-bell patterns (rise then fall or fall then rise). Exports `Time_Bell_Curves.xlsx` and plots top results.

## Workflow

1. Ensure raw source Excel file is next to the scripts (default file names are set in each script).
2. Run preprocessing to generate `Combined_Data.xlsx`:

```powershell
python Preprocessing.py
```

3. Run the VAS correlation analysis:

```powershell
python VAS_Correlation_Analysis.py
```

4. Run the time/bell-curve analysis:

```powershell
python Time_Analysis.py
```

## Dependencies

Python packages: `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`, `statsmodels`, `openpyxl` (for Excel I/O). Install with:

```powershell
pip install pandas numpy scipy matplotlib seaborn statsmodels openpyxl
```

## Input / Output

- Input: Script-specific Excel files (defaults in each script). `Preprocessing.py` expects `Acute_Raw-Abundance_VAS.xlsx` and writes `Combined_Data.xlsx`.
- Outputs:
  - `Combined_Data.xlsx` (from preprocessing)
  - `Significant_VAS_Correlations.xlsx` (from VAS correlation analysis)
  - `Time_Bell_Curves.xlsx` (from time analysis)
  - Plots are shown interactively via matplotlib/seaborn when analyses run.

## Data Handling & Filtering Decisions

### `Preprocessing.py`

**Input:** Raw Excel file `Acute_Raw-Abundance_VAS.xlsx` with separate sheets for VAS scores, soluble fraction, SPE fraction, and insoluble fraction.

**Processing Steps:**
1. **Sample ID cleanup:** Removes fraction suffixes (_i, _b) from sample IDs.
2. **Drop Group column:** Removes the Group column from protein fraction sheets.
3. **Collapse replicates:** Averages values across replicates within each Sample ID + Condition group.
4. **Merge fractions:** Combines all fractions (soluble, SPE, insoluble) with VAS data on Sample ID and Condition.
5. **Deduplicate proteins:** Proteins found in multiple fractions receive suffixes (_s, _spe, _i).
6. **Remove time suffix:** Strips time indicators (e.g., -T10) from Sample IDs to consolidate samples.
7. **Sort columns:** Arranges output with metadata, VAS, then proteins alphabetically.

**Output:** `Combined_Data.xlsx` with one row per Sample ID × Condition combination and all protein measurements.

---

### `VAS_Correlation_Analysis.py`

**Input:** `Combined_Data.xlsx`

**Filtering & Thresholds:**
- **P-value threshold:** p < 0.05 for significance.
- **Minimum datapoints:** ≥ 4 valid (non-missing) datapoints per protein–participant pair.
- **Unique values:** ≥ 2 unique protein values (excludes constant measurements).
- **Zero exclusion:** Skips any protein measurement with a zero value.
- **Relative range check:** Protein must show ≥ 20% relative range: `(max – min) / min ≥ 0.2`. This filters out proteins with low variation (noise).
- **Perfect correlation exclusion:** Correlations with |r| ≥ 0.998 are excluded (typically artifacts from missing data patterns).
- **Only significant correlations are outputted:** Some proteins may have significant correlations with individual participants, but not overall (or vice versa), and only the results with significant correlations will be included in the output file.

**Analysis:**
- Computes Pearson correlation between each protein and VAS for:
- Each participant individually.
- All participants combined.
- Significant results are aggregated and ranked by frequency (how often each protein appears across all correlations).

**Output:** 
- `Significant_VAS_Correlations.xlsx`: All significant correlations sorted by frequency then correlation strength.
- Interactive plots: Scatter + regression line for top 16 correlations.

---

### `Time_Analysis.py`

**Input:** `Combined_Data.xlsx`

**Filtering & Thresholds:**
- **Protein columns:** Exclude proteins with ≥ 3 missing values across all rows.
- **Time mapping:** Converts Condition (T0, T10, T20, T30) to numeric time in minutes.
- **Complete timepoints:** Bell curve detection requires all four timepoints (T0, T10, T20, T30) present.
- **P-value threshold:** p < 0.05 for quadratic term significance. --REMOVED REQUIRMENT TO FOCUS ONLY ON BELL STRENGTH. Can be reintroduced by adding back 'if a is not None and pval < 0.05' to the 'Quadratics and bell curves' sections

**Analysis:**
- Fits a quadratic model to each protein across timepoints for:
  - Each participant individually.
  - All participants combined.
- **Bell curve detection:** Identifies patterns where:
  - **Normal bell:** protein rises T0→T10/T20, then falls by T30.
  - **Inverted bell:** protein falls T0→T10/T20, then rises by T30.
- **Bell strength:** Measures magnitude as `(peak – end_average) / end_average`, normalised to endpoints.

**Output:**
- `Time_Bell_Curves.xlsx`: Significant quadratic fits with bell curve patterns and strength scores.
- Interactive plots: Line plots for top 16 strongest bell curves.

---

## Notes

- Run `Preprocessing.py` first to ensure analyses use the combined dataset.
- Edit the `FILE_PATH` constant at the top of each script if your input filenames differ.
- Adjust key thresholds by editing top-level variables:
  - `P_THRESHOLD = 0.05` (significance level)
  - `REALTIVE_RANGE_THRESHOLD = 0.2` (VAS correlation, minimum 20% protein variation)
  - `NUMBER_TO_PLOT = 16` (number of top results to visualize)
- All scripts are study-specific and would need to be adjusted to apply to new datasets.
- Due to missing data, 575 proteins are excluded in the time analysis and a similar number in the correlation analysis.
