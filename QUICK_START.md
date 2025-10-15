# BEER Python - Quick Start Guide

## 5-Minute Quick Start

### Option 1: Your data is already log-normalized (most common)

```python
import anndata as ad
from beer import BEER

# Load your preprocessed data
adata = ad.read_h5ad('your_data.h5ad')

# Run BEER (data is already normalized)
beer = BEER(n_pcs=50, n_groups=30)
result = beer.fit_transform(adata, batch_key='batch', normalize=False)

# Done! Access results
corrected_adata = result.adata
selected_pcs = result.selected_pcs
```

### Option 2: Your data is raw counts

```python
import anndata as ad
from beer import BEER

# Load raw counts
adata = ad.read_h5ad('raw_counts.h5ad')

# Run BEER (will normalize automatically)
beer = BEER(n_pcs=50, n_groups=30)
result = beer.fit_transform(adata, batch_key='batch', normalize=True)

# Done!
corrected_adata = result.adata
```

## Quick Reference

### When to use `normalize=True` vs `normalize=False`

| Your Data Type | Use This |
|----------------|----------|
| ✅ Already log-normalized (Scanpy preprocessed) | `normalize=False` |
| ✅ Already log-normalized (from R/Seurat) | `normalize=False` |
| ✅ Raw UMI counts | `normalize=True` |
| ✅ FPKM/RPKM/TPM | `normalize=False` |
| ❌ Scaled (z-scores) | Not supported - use log-normalized instead |

### Basic Parameters

```python
BEER(
    n_pcs=50,              # Number of PCs to compute
    n_groups=30,           # Groups per batch (k-means)
    n_variable_genes=2000, # Variable genes to use
    n_rounds=1,            # MN detection rounds (1-3)
    use_combat=True,       # ComBat correction
    random_seed=123        # For reproducibility
)
```

### Common Use Cases

#### Case 1: Fix the shape mismatch error

If you see: `ValueError: operands could not be broadcast together`

**Solution:** Your data is already normalized. Use `normalize=False`:

```python
result = beer.fit_transform(adata, batch_key='batch', normalize=False)
```

#### Case 2: Integrate multiple batches

```python
beer = BEER(n_pcs=50, n_groups=30, n_rounds=2)  # Increase rounds
result = beer.fit_transform(adata, batch_key='batch', normalize=False)
```

#### Case 3: Remove specific genes (e.g., cell cycle)

```python
cell_cycle_genes = ['PCNA', 'MKI67', 'TOP2A', ...]
result = beer.fit_transform(
    adata,
    batch_key='batch',
    normalize=False,
    remove_genes=cell_cycle_genes
)
```

#### Case 4: BBKNN enhancement for better mixing

```python
from beer import BEER, apply_bbknn

beer = BEER(n_pcs=50, n_groups=30)
result = beer.fit_transform(adata, batch_key='batch', normalize=False)

# Apply BBKNN
umap_bbknn = apply_bbknn(result.adata, result.selected_pcs, batch_key='batch')
result.adata.obsm['X_umap_bbknn'] = umap_bbknn
```

## Installation

```bash
pip install -r requirements_python.txt
```

## Verify It's Working

```python
import numpy as np
import anndata as ad
from beer import BEER

# Create test data
X = np.random.randn(500, 1000)
batch = np.array(['A'] * 250 + ['B'] * 250)
adata = ad.AnnData(X)
adata.obs['batch'] = batch

# Run BEER
beer = BEER(n_pcs=20, n_groups=10)
result = beer.fit_transform(adata, batch_key='batch', normalize=False)

print(f"✓ BEER works! Selected {len(result.selected_pcs)} PCs")
```

## Next Steps

- **Full documentation:** [README_PYTHON.md](README_PYTHON.md)
- **Pre-normalized data guide:** [PRE_NORMALIZED_DATA_USAGE.md](PRE_NORMALIZED_DATA_USAGE.md)
- **More examples:** [example_usage.py](example_usage.py)
- **Installation help:** [INSTALL_PYTHON.md](INSTALL_PYTHON.md)

## Common Errors & Quick Fixes

| Error | Fix |
|-------|-----|
| `ValueError: operands could not be broadcast` | Use `normalize=False` |
| `ValueError: Need at least 2 batches` | Check your batch_key is correct |
| `ImportError: combat-python not installed` | `pip install combat-python` or set `use_combat=False` |
| `No variable genes found` | Check your data format, use `normalize=False` if pre-normalized |

## Need Help?

1. Check if your data is log-normalized: `print(adata.X.max())` - should be < 20
2. Try `normalize=False` if data is already processed
3. See [PRE_NORMALIZED_DATA_USAGE.md](PRE_NORMALIZED_DATA_USAGE.md) for detailed guide
4. Open an issue on GitHub with your error message

---

**That's it!** BEER is ready to use. 🍺
