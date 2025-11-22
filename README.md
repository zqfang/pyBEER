# BEER Python Implementation

**BEER** (Batch EffEct Remover) - Python implementation for removing batch effects from single-cell data using mutual nearest neighbors and PC selection.

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This is a Python port of the BEER algorithm originally implemented in R. BEER removes batch effects in single-cell RNA-seq data by:

1. **Detecting batch-specific principal components** using mutual nearest neighbors (MN pairs)
2. **Selecting biologically relevant PCs** while excluding batch-effect-dominated PCs
3. **Correcting expression values** using ComBat (optional)

**Original Publication:**
> Zhang, F., Wu, Y., & Tian, W. (2019). A novel approach to remove the batch effect of single-cell data. *Cell Discovery*, 5, 46. https://doi.org/10.1038/s41421-019-0114-x

## Key Features
- **Detect Batch Effect quantitatively**, instead of visual inspection
- **Modern Python implementation** using NumPy, SciPy, AnnData, and Scanpy
- **Functional API** matching the original R interface
- **Easy integration** with Scanpy workflows
- **ComBat integration** for expression correction (optional)

## Installation

### Requirements

- Python >= 3.8
- numpy >= 1.20.0
- scipy >= 1.7.0
- pandas >= 1.3.0
- scikit-learn >= 1.0.0
- anndata >= 0.8.0
- scanpy >= 1.9.0
- statsmodels >= 0.13.0

### Install Dependencies

```bash
pip install numpy scipy pandas scikit-learn anndata scanpy statsmodels
```

### Optional Dependencies

For ComBat correction:
```bash
pip install combat-python
```

## Quick Start

### Basic Usage

```python
import scanpy as sc
from beer import beer, get_use, select_use

# Load your data
adata = sc.read_h5ad('your_data.h5ad')

# Run BEER
result = beer(
    adata,
    batch_key='batch',
    gnum=30,          # groups per batch
    pcnum=50,         # number of PCs
    gn=2000,          # variable genes per batch
    combat=False,     # set True for ComBat correction
    seed=123,
    n_rounds=1
)

# Access results
selected_pcs = result.select  # 0-based indices
corrected_adata = result.adata

# Compute UMAP with selected PCs
sc.pp.neighbors(result.adata, n_pcs=len(selected_pcs), use_rep='X_pca')
sc.tl.umap(result.adata)

# Visualize
sc.pl.umap(result.adata, color='batch')

# View correlations
print(f"Spearman correlations: {result.cor}")
print(f"Pearson correlations: {result.lcor}")
```

### Custom PC Selection

```python
# Strict threshold selection
strict_pcs = get_use(result, cut_r=0.8, cut_l=0.8)

# With rank ratios
custom_pcs = select_use(result, cut_r=0.7, cut_l=0.7, rr=0.6, rl=0.6)
```

### Multiple Batches

```python
# Handle multiple batches (>2)
result = beer(
    adata,
    batch_key='batch',
    gnum=30,
    pcnum=50,
    n_rounds=2,  # Increase rounds for better cross-batch pairing
    combat=False
)
```

### Remove Specific Genes

```python
# Remove cell cycle genes from analysis
cell_cycle_genes = ['MKI67', 'TOP2A', 'PCNA']
result = beer(
    adata,
    batch_key='batch',
    rmg=cell_cycle_genes
)
```

## Result Object: `BEERResult`

```python
@dataclass
class BEERResult:
    adata: AnnData          # Processed AnnData with PCA/UMAP
    vp: np.ndarray          # Valid MN pairs (n_pairs x 2)
    cor: np.ndarray         # Spearman correlations per PC
    pv: np.ndarray          # P-values for Spearman
    fdr: np.ndarray         # FDR-adjusted p-values
    lcor: np.ndarray        # Pearson correlations per PC
    lpv: np.ndarray         # P-values for Pearson
    lc1: np.ndarray         # Linear model intercept p-values
    lc2: np.ndarray         # Linear model slope p-values
    lfdr: np.ndarray        # FDR-adjusted Pearson p-values
    select: np.ndarray      # Selected PC indices (0-based)

    # Parameters used
    round: int
    combat: bool
    gnum: int
    gn: int
    pcnum: int
    seed: int
```

## API Reference

### Main Function

```python
beer(
    adata,                    # AnnData object (cells x genes)
    batch_key='batch',        # Key in adata.obs for batch labels
    gnum=30,                  # Number of groups per batch
    pcnum=50,                 # Number of PCs to compute
    gn=2000,                  # HVGs per batch
    combat=False,             # Apply ComBat correction
    seed=123,                 # Random seed
    n_components=2,           # UMAP dimensions
    n_rounds=1,               # MN detection rounds
    rmg=None                  # Genes to remove
) -> BEERResult
```

### Helper Functions

```python
# Get PCs above correlation cutoffs
get_use(result, cut_r=0.7, cut_l=0.7) -> np.ndarray

# Select PCs with rank and correlation criteria
select_use(result, cut_r=0.7, cut_l=0.7, rr=0.5, rl=0.5) -> np.ndarray
```

## Common Issues

| Issue | Symptom | Solution |
|-------|---------|----------|
| Too few PCs selected | < 10 PCs selected | Increase `gnum` or lower thresholds |
| Over-correction | Biological signal lost | Decrease `n_rounds`, use `rmg` |
| Under-correction | Batches still separated | Increase `n_rounds` |
| Memory limitations | Out of memory | Reduce `pcnum` or `gn` |

## Workflow Integration

### Scanpy Pipeline

```python
import scanpy as sc
from beer import beer

# Standard Scanpy preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Run BEER (will normalize internally)
result = beer(adata, batch_key='batch', gnum=30, pcnum=50)

# Continue with Scanpy workflow using selected PCs
sc.pp.neighbors(result.adata, n_pcs=len(result.select), use_rep='X_pca')
sc.tl.umap(result.adata)
sc.tl.leiden(result.adata)

# Visualize
sc.pl.umap(result.adata, color=['batch', 'leiden'])
```

### Save/Load Results

```python
# Save corrected data
result.adata.write_h5ad('beer_corrected.h5ad')

# Load for downstream analysis
adata_corrected = sc.read_h5ad('beer_corrected.h5ad')
```

## Examples

See [example_usage.py](example_usage.py) for comprehensive examples:

1. **Basic Usage**: Simple batch correction with simulated data
2. **Real Data**: Analysis with PBMC3k dataset
3. **Multiple Batches**: Handling >2 batches
4. **Custom PC Selection**: Manual threshold adjustment
5. **Gene Removal**: Excluding specific genes
6. **Detailed Results**: Exploring all output fields

Run all examples:
```bash
python example_usage.py
```

## Comparison with R Version

| Parameter | R | Python |
|-----------|---|--------|
| Data orientation | genes x cells | cells x genes (AnnData) |
| PC indices | 1-based | 0-based |
| Result object | list | BEERResult dataclass |
| `cor` | Spearman correlations | `result.cor` |
| `lcor` | Pearson correlations | `result.lcor` |
| `select` | Selected PCs | `result.select` |
| `vp` | Valid pairs | `result.vp` |

## Resources

- **Original R Implementation**: https://github.com/jumphone/BEER
- **Paper**: https://doi.org/10.1038/s41421-019-0114-x
- **Scanpy**: https://scanpy.readthedocs.io/
- **AnnData**: https://anndata.readthedocs.io/

## Acknowledgments

- **Original Author**: Feng Zhang
- **Original R Implementation**: https://github.com/jumphone/BEER
- **Python Port**: Based on BEER v0.1.9

## Changelog

### v0.2.1 (November 2024)
- Simplified functional API matching R version
- Uses `beer()` function instead of class
- Result fields match R: `cor`, `lcor`, `vp`, `select`
- Added `get_use()` and `select_use()` helper functions
- Fixed matrix orientation handling
- Improved documentation

### v0.2.0 (October 2024)
- Initial Python implementation
- Full compatibility with Scanpy/AnnData ecosystem
