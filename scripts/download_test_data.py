#!/usr/bin/env python3
"""Download test dataset for pipeline testing."""

import scanpy as sc
from pathlib import Path

print("Downloading PBMC 3k dataset...")
adata = sc.datasets.pbmc3k()

print(f"Downloaded {adata.n_obs} cells × {adata.n_vars} genes")

# Basic preprocessing
print("Preprocessing...")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print(f"After filtering: {adata.n_obs} cells × {adata.n_vars} genes")

# Add ensembl_id - use the actual Ensembl IDs from gene_ids column
# Tahoe model uses Ensembl gene IDs (ENSG format)
adata.var['ensembl_id'] = adata.var['gene_ids']

# Add some fake cell type annotations for testing
# Just use simple binning for quick testing
print("Adding cell type annotations...")
import numpy as np
n_cells = adata.n_obs
cell_types = ['T cells', 'B cells', 'Monocytes', 'NK cells', 'Other']
adata.obs['cell_type'] = np.random.choice(cell_types, size=n_cells)

print(f"Cell types: {adata.obs['cell_type'].unique()}")

# Save
output_path = Path("data/pbmc3k_test.h5ad")
output_path.parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(output_path)

print(f"✓ Saved to {output_path}")
print(f"\nDataset info:")
print(f"  Cells: {adata.n_obs}")
print(f"  Genes: {adata.n_vars}")
print(f"  Cell types: {len(adata.obs['cell_type'].unique())}")
