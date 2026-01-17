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
# 4. VISUALIsATION
# ==========================================

plt.figure(figsize=(12, 5))

# Plot 1: Coloired by Biological Factor (Condition)
plt.subplot(1, 2, 1)
sns.scatterplot(
    data=final_df, x='PC1', y='PC2', 
    hue='Condition', style='Timepoint', s=100
)
plt.title('PCA: Biological Factors (Condition/Time)')
plt.grid(True, alpha=0.3)

# # Plot 2: Coloured by Replicate: separation = batch effect
# plt.subplot(1, 2, 2)
# sns.scatterplot(
#     data=final_df, x='PC1', y='PC2', 
#     hue='Replicate', style='Timepoint', 
#     palette='Set1', s=100
# )
# plt.title('PCA: Technical Factors (Replicate Check)')
# plt.grid(True, alpha=0.3)

# plt.tight_layout()
# plt.show()

# Plot 2: Coloured by Replicate: separation = batch effect
plt.subplot(1, 2, 2)
sns.scatterplot(
    data=final_df, x='PC1', y='PC2', 
    hue='Study', style='Timepoint', 
    palette='Set1', s=100
)
plt.title('PCA: Technical Factors (Replicate Check)')
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# ==========================================
# 5. CLUSTERMAP (Hierarchical Clustering)
# ==========================================

# Calculate correlation matrix (transposed so correlate samples, not proteins)
corr_matrix = df_log[expression_cols].T.corr(method='pearson')

# Map metadata to colours for  sidebar
condition_colors = df_log['Condition'].map({'Resp': 'blue', 'NonResp': 'red'})
replicate_colors = df_log['Replicate'].map({1: 'green', 2: 'orange'})

# Plot
sns.clustermap(
    corr_matrix, 
    row_colors=[condition_colors, replicate_colors],
    col_colors=[condition_colors, replicate_colors],
    cmap='viridis',
    figsize=(8, 8),
    annot=True
)
plt.title("Sample Correlation Matrix")
plt.show()