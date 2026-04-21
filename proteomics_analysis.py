# ===============================
# Proteomics Analysis Script
# Heatmap + Volcano Plot
# ===============================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------
# 1. Load Data
# -------------------------------
file_path = "Juhi_Control_Vs_Treated_Proteins.xlsx"
df = pd.read_excel(file_path)

# -------------------------------
# 2. Compute metrics
# -------------------------------
df['log2FC'] = np.log2(df['Abundance Ratio: (Treated) / (Control)'])
df['neglog10p'] = -np.log10(df['Abundance Ratio Adj. P-Value: (Treated) / (Control)'])

# Remove invalid values
df.replace([np.inf, -np.inf], np.nan, inplace=True)
df.dropna(subset=['log2FC', 'neglog10p'], inplace=True)

# -------------------------------
# 3. Define protein subset (heatmap proteins)
# -------------------------------
proteins = [
    "PPP2CA","ARG2","STK25","EZR","RPL38","CLIC4","CASP3","BABAM1",
    "PALLD","KIF2A","CHEK2","SUGT1","ACTG1","CLTCL1","VPS4B","CDC23","ADAMTS9"
]

subset = df[df['Gene Symbol'].isin(proteins)].copy()

# -------------------------------
# 4. Add negative control
# -------------------------------
neg_control = df[
    (df['log2FC'].between(-1, 1)) &
    (df['neglog10p'] < 1.3)
].head(1)

subset = pd.concat([subset, neg_control])

# -------------------------------
# 5. Volcano Plot
# -------------------------------

plt.figure(figsize=(6,5))

# Background (all proteins)
plt.scatter(df['log2FC'], df['neglog10p'], alpha=0.2)

# Color coding
colors = []
for _, r in subset.iterrows():
    if r['log2FC'] > 1:
        colors.append('red')
    elif r['log2FC'] < -1:
        colors.append('blue')
    else:
        colors.append('gray')

# Plot subset
plt.scatter(subset['log2FC'], subset['neglog10p'], c=colors)

# Key proteins (repelled labels)
key_proteins = ["CASP3", "EZR", "ACTG1"]

offsets = {
    "CASP3": (0.3, 0.6),
    "EZR": (-0.4, 0.6),
    "ACTG1": (0.3, -0.6)
}

for _, row in subset.iterrows():
    name = row['Gene Symbol']
    if name in key_proteins:
        dx, dy = offsets[name]
        plt.annotate(
            name,
            (row['log2FC'], row['neglog10p']),
            xytext=(row['log2FC'] + dx, row['neglog10p'] + dy),
            arrowprops=dict(arrowstyle='-'),
            fontsize=8
        )

# Threshold lines
plt.axvline(x=1, linestyle='--')
plt.axvline(x=-1, linestyle='--')
plt.axhline(y=1.3, linestyle='--')

plt.xlabel("log2 Fold Change")
plt.ylabel("-log10 p-value")
plt.title("Volcano Plot (Limb-development proteins)")

plt.tight_layout()
plt.savefig("volcano_plot.png", dpi=300)
plt.close()

# -------------------------------
# 6. Heatmap (Z-score normalized)
# -------------------------------

# Extract abundance data
heat = subset[[
    'Abundances (Grouped): Control',
    'Abundances (Grouped): Treated'
]].astype(float)

# Clean values
heat.replace([np.inf, -np.inf], np.nan, inplace=True)
heat.fillna(heat.mean(), inplace=True)

# Row-wise z-score
heat_z = (heat - heat.mean(axis=1).values.reshape(-1,1)) / heat.std(axis=1).values.reshape(-1,1)
heat_z = np.nan_to_num(heat_z)

labels = subset['Gene Symbol'].fillna(subset['Accession'])

# Plot heatmap
plt.figure(figsize=(5,10))
plt.imshow(heat_z, aspect='auto', cmap='bwr', vmin=-2, vmax=2)

plt.yticks(range(len(labels)), labels, fontsize=6)
plt.xticks([0,1], ["Control", "Treated"])

cbar = plt.colorbar()
cbar.set_label("Row-wise z-score (mean-centered)")

plt.title("Heatmap (Limb-development proteins)")

plt.tight_layout()
plt.savefig("heatmap.png", dpi=300)
plt.close()

# -------------------------------
# DONE
# -------------------------------
print("Analysis complete. Files saved:")
print(" - volcano_plot.png")
print(" - heatmap.png")