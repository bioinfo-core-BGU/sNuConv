import sys
import pandas as pd
import anndata

# Variables
Correlations = sys.argv[1]
pseudo_h5 = sys.argv[2]
outDir = sys.argv[3]

# Import data
pseudo_h5ad = anndata.read_h5ad(pseudo_h5)
df_corr = pd.read_table(Correlations)
df_corr.index = df_corr['Gene']

# Normalize pseudo counts
pseudo_counts = pseudo_h5ad.to_df().transpose()
pseudo_counts = (pseudo_counts / pseudo_counts.sum(axis=0))*1000000
pseudo_counts = pseudo_counts[pseudo_counts.index.isin(df_corr.index)]

# Correct Counts
pseudo_corrected_counts = pseudo_counts.copy(deep=True)
for gene in df_corr.index:
	pseudo_corrected_counts.loc[gene,] = pseudo_counts.loc[gene] * df_corr.loc[gene,'Slope'] + df_corr.loc[gene,'Intercept']
pseudo_corrected_counts[pseudo_corrected_counts<0] = 0
pseudo_corrected_counts = pseudo_corrected_counts.transpose()

# Save results
corrected_h5ad = anndata.AnnData(X=pseudo_corrected_counts, obs=pseudo_h5ad.obs, uns=pseudo_h5ad.uns)
corrected_h5ad.write_h5ad(outDir + 'training_data_corrected.h5ad')