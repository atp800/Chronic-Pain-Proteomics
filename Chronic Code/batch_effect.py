import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer

# ==========================================
# 1. SETUP
# ==========================================
# File Paths
INPUT_FILE_PATH  = "Chronic Code/Data_To_Clean.xlsx"                         # NEED TO ADD BACK IN GROUP AND TIME COLUMNS (could get time from sample ID)
SHEET_NAME       = "Soluble - Group"
DATA = pd.read_excel(INPUT_FILE_PATH, sheet_name=SHEET_NAME)


df = pd.DataFrame(DATA)

# ==========================================
# 2. PREPROCESSING
# ==========================================

# Extract Timepoint and Subject from Sample_ID
# (Assumes format "Subject-Timepoint")
df['Timepoint'] = df['Sample_ID'].apply(lambda x: x.split('-')[1])
df['Subject'] = df['Sample_ID'].apply(lambda x: x.split('-')[0])

df['Study'] = df['Sample_ID'].apply(lambda x: x.split('P')[0])

# Separate Metadata from Expression Data
metadata_cols = ['Study','Sample_ID', 'Condition', 'Replicate', 'Timepoint', 'Subject']
expression_cols = [c for c in df.columns if c not in metadata_cols]

# Logify
epsilon = 1e-9 
df_log = df.copy()
df_log[expression_cols] = np.log2(df[expression_cols] + epsilon)

# Impute missing vals
imputer = SimpleImputer(strategy='mean') 
df_log[expression_cols] = imputer.fit_transform(df_log[expression_cols])

# Standardisation (Z-score)
scaler = StandardScaler()
scaled_data = scaler.fit_transform(df_log[expression_cols])

# ==========================================
# 3. PCA (Principal Component Analysis)
# ==========================================

# Run PCA
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(scaled_data)
pca_df = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])

# Combine PCA results with metadata for plotting
final_df = pd.concat([pca_df, df_log[metadata_cols].reset_index(drop=True)], axis=1)

# Calculate explained variance
exp_var = pca.explained_variance_ratio_
print(f"Explained Variance: PC1={exp_var[0]*100:.2f}%, PC2={exp_var[1]*100:.2f}%")

# ==========================================
# 4. VISUALISATION
# ==========================================

# Increased figure size to accommodate 3 plots side-by-side (18x6 instead of 12x5)
plt.figure(figsize=(18, 6))

# Plot 1: Coloured by Biological Factor (Condition)
plt.subplot(1, 3, 1)
sns.scatterplot(
    data=final_df, x='PC1', y='PC2', 
    hue='Condition', style='Timepoint', s=100
)
plt.title('PCA: Biological (Condition)')
plt.grid(True, alpha=0.3)

# Plot 2: Coloured by Study (technical 1)
plt.subplot(1, 3, 2)
sns.scatterplot(
    data=final_df, x='PC1', y='PC2', 
    hue='Study', style='Timepoint', 
    palette='Set1', s=100
)
plt.title('PCA: Batch Effect (Study)')
plt.grid(True, alpha=0.3)

# Plot 3: Coloured by Replicate (technical 2 2)
plt.subplot(1, 3, 3)
sns.scatterplot(
    data=final_df, x='PC1', y='PC2', 
    hue='Replicate', style='Timepoint', 
    palette='Set2', s=100
)
plt.title('PCA: Batch Effect (Replicate)')
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# ==========================================
# 5. CLUSTERMAP (Hierarchical Clustering)
# ==========================================

# # Calculate correlation matrix
# corr_matrix = df_log[expression_cols].T.corr(method='pearson')


# # 1. Condition Colors
# condition_colors = df_log['Condition'].map({'Resp': 'blue', 'NonResp': 'red'})

# # 2. Replicate Colors (Manual mapping)
# replicate_colors = df_log['Replicate'].map({1: 'green', 2: 'orange'})

# # 3. Study Colors (Automated mapping using a palette)
# unique_studies = df_log['Study'].unique()
# study_palette = sns.color_palette("husl", len(unique_studies))
# study_lut = dict(zip(unique_studies, study_palette))
# study_colors = df_log['Study'].map(study_lut)

# # Combine into a dataFrame for the row/col colors
# row_colors_df = pd.DataFrame({
#     'Condition': condition_colors,
#     'Replicate': replicate_colors,
#     'Study': study_colors
# })

# # Plot
# sns.clustermap(
#     corr_matrix, 
#     row_colors=row_colors_df,
#     col_colors=row_colors_df,
#     cmap='viridis',
#     figsize=(10, 10),
#     annot=False # I turned off annotations (numbers) as they are messy on large plots
# )
# plt.title("Sample Correlation Matrix")
# plt.show()