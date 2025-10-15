# BEER Python - Frequently Asked Questions (FAQ)

## Table of Contents

1. [Understanding BEER Algorithm](#understanding-beer-algorithm)
2. [Data Format Questions](#data-format-questions)
3. [Parameter Tuning](#parameter-tuning)
4. [Common Errors](#common-errors)
5. [Results Interpretation](#results-interpretation)
6. [Performance Questions](#performance-questions)
7. [Integration with Other Tools](#integration-with-other-tools)

---

## Understanding BEER Algorithm

### Q: What do rank_correlations and linear_correlations mean?

**A: These metrics evaluate whether each PC captures biological signal or batch effects.**

#### The Core Concept

BEER evaluates each principal component (PC) to determine if it's dominated by **batch effects** or **biological signal**. It does this by measuring how similar matched cell groups are across different batches.

#### What Are We Correlating?

For each PC, BEER:
1. Takes **mutual nearest neighbor (MN) pairs** - groups of cells from different batches that are biologically similar
2. For each paired group, calculates **quantiles** (0%, 25%, 50%, 75%, 100%) of the PC values
3. Correlates these quantiles between the paired groups

#### Rank Correlation (Spearman)

**What it measures:** Monotonic relationship between groups

```python
# Rank correlation (Spearman)
rank_cor, rank_pv = stats.spearmanr(lst1_quantiles, lst2_quantiles)
```

**Interpretation:**
- **High rank correlation (>0.7)**: The PC shows similar **ordering patterns** across batches
  - Groups that have high values in batch A also have high values in batch B
  - This suggests the PC captures **biological variation**, not batch effects

- **Low rank correlation (<0.5)**: The PC shows different ordering patterns across batches
  - Groups are ordered differently in each batch
  - This suggests the PC is **dominated by batch effects**

**Example:**
```
Batch A groups: [0.5, 1.2, 2.3, 3.1]  (ranked: 1, 2, 3, 4)
Batch B groups: [0.8, 1.5, 2.6, 3.4]  (ranked: 1, 2, 3, 4)
→ High rank correlation = biological signal preserved
```

#### Linear Correlation (Pearson)

**What it measures:** Linear relationship between groups

```python
# Linear correlation (Pearson)
linear_cor, linear_pv = stats.pearsonr(lst1_quantiles, lst2_quantiles)
```

**Interpretation:**
- **High linear correlation (>0.7)**: The PC shows similar **absolute values** across batches
  - Not just the same ordering, but similar magnitudes
  - Even stronger evidence of **biological signal**

- **Low linear correlation (<0.5)**: The PC shows different value scales across batches
  - Values might be shifted or scaled differently
  - Suggests **batch-specific technical effects**

**Example:**
```
Batch A groups: [1.0, 2.0, 3.0, 4.0]
Batch B groups: [1.1, 2.1, 2.9, 4.2]
→ High linear correlation = values match closely
```

#### Why Use Both?

BEER uses **both** metrics because they capture different aspects:

| Metric | What It Catches | Example Batch Effect |
|--------|----------------|---------------------|
| **Rank correlation** | Changes in ordering | PC1 separates cell types A,B,C in batch 1, but separates B,C,A in batch 2 |
| **Linear correlation** | Changes in scale/shift | PC1 ranges 0-5 in batch 1, but 10-15 in batch 2 (same biology, shifted values) |

#### PC Selection Logic

BEER selects PCs with **BOTH** high rank AND high linear correlation:

```python
# Selection criteria
rank_criterion = (
    (rank(correlations) >= n_pcs * 0.5) |  # Top 50% by rank
    (correlations > 0.7)                    # OR above threshold
)

linear_criterion = (
    (rank(linear_correlations) >= n_pcs * 0.5) |  # Top 50% by rank
    (linear_correlations > 0.7)                    # OR above threshold
)

selected_pcs = rank_criterion & linear_criterion  # Must pass BOTH
```

#### Visual Example

```
PC   Rank Cor  Linear Cor  Interpretation
--   --------  ----------  ---------------
1    0.85      0.82        ✅ SELECTED - Biological signal
2    0.92      0.89        ✅ SELECTED - Strong biological signal
3    0.45      0.38        ❌ REJECTED - Batch effect dominated
10   0.78      0.15        ❌ REJECTED - Scale differs (batch effect)
15   0.22      0.75        ❌ REJECTED - Ordering differs (batch effect)
```

#### Biological Intuition

**Think of it like comparing recipes across kitchens:**

- **Rank correlation**: Do both kitchens use ingredients in the same order?
  - Kitchen A: salt → pepper → garlic
  - Kitchen B: salt → pepper → garlic
  - ✅ High rank correlation = same recipe structure

- **Linear correlation**: Do they use the same amounts?
  - Kitchen A: 1tsp salt, 1tsp pepper, 2tsp garlic
  - Kitchen B: 1tsp salt, 1tsp pepper, 2tsp garlic
  - ✅ High linear correlation = same recipe proportions

#### Summary Table

| Rank Cor | Linear Cor | Interpretation | Action |
|----------|-----------|----------------|--------|
| High (>0.7) | High (>0.7) | Strong biological signal | ✅ Select |
| High | Medium | Biological with some batch variance | ⚠️ Consider |
| High | Low | Same pattern, different scale | ❌ Reject (batch) |
| Low | High | Different pattern, similar scale | ❌ Reject (batch) |
| Low | Low | Dominated by batch effects | ❌ Reject |

---

### Q: What are mutual nearest neighbors (MN pairs)?

**A:** Mutual nearest neighbors are groups of cells from different batches that are biologically similar to each other.

**How they're found:**
1. Cells are grouped within each batch using k-means clustering (on UMAP coordinates)
2. Groups are aggregated by summing expression across cells
3. Correlation between all group pairs is computed
4. A pair (GroupA, GroupB) is an MN pair if:
   - GroupA's most correlated partner is GroupB
   - GroupB's most correlated partner is GroupA
   - They are from different batches

**Why they matter:**
- MN pairs represent the same cell type/state across batches
- Used to evaluate if PCs preserve biological similarity across batches
- More MN pairs = better batch integration

**Example:**
```
Batch 1: [CD4_T_cells, CD8_T_cells, B_cells, Monocytes]
Batch 2: [CD4_T_cells, CD8_T_cells, B_cells, Monocytes]

MN pairs might be:
- Batch1_CD4_T → Batch2_CD4_T  (matched T cells)
- Batch1_B_cells → Batch2_B_cells  (matched B cells)
```

---

### Q: How does BEER differ from other batch correction methods?

**A:** BEER uses a unique PC selection approach rather than directly correcting gene expression.

| Method | Approach | Output |
|--------|----------|--------|
| **BEER** | Select PCs with low batch correlation | Selected PCs for analysis |
| **Harmony** | Iteratively adjust PC embeddings | Corrected PC embeddings |
| **Seurat/CCA** | Find anchors via canonical correlation | Integrated expression matrix |
| **ComBat** | Empirical Bayes correction of expression | Corrected expression matrix |
| **scVI** | Deep learning variational inference | Latent space representation |

**BEER advantages:**
- ✅ Interpretable (you can inspect which PCs are kept/removed)
- ✅ Fast (no iterative optimization)
- ✅ Works well with strong batch effects
- ✅ Optional ComBat integration for enhanced correction

**When to use BEER:**
- You want to understand which PCs capture batch effects
- You need fast batch correction
- You're integrating 2-10 batches
- You want control over the correction strength

---

## Data Format Questions

### Q: Is my data already log-normalized? How can I check?

**A:** Here's how to check if your data is raw counts or log-normalized:

```python
import numpy as np

# Check 1: Value range
print(f"Min value: {adata.X.min()}")
print(f"Max value: {adata.X.max()}")
print(f"Mean value: {adata.X.mean()}")

# Interpretation:
# Raw counts: max > 1000, mean varies widely
# Log-normalized: max < 20, mean ~1-2
```

**Quick decision guide:**

| Your Data | Max Value | Mean Value | Use |
|-----------|-----------|------------|-----|
| Raw counts | > 1000 | Varies | `normalize=True` |
| Log-normalized | < 20 | ~1-2 | `normalize=False` |
| Scaled (z-scores) | ~5-10 | ~0 | ❌ Not supported |

**Visual check:**
```python
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 4))

plt.subplot(1, 2, 1)
plt.hist(adata.X.data if hasattr(adata.X, 'data') else adata.X.flatten(),
         bins=50, log=True)
plt.xlabel('Expression value')
plt.title('Value distribution')

plt.subplot(1, 2, 2)
cell_sums = np.array(adata.X.sum(axis=1)).flatten()
plt.hist(cell_sums, bins=50)
plt.xlabel('Sum per cell')
plt.title('Cell library sizes')

plt.tight_layout()
plt.show()
```

**If you see:**
- Wide variation in cell sums → Raw counts
- Similar cell sums → Log-normalized

---

### Q: Should I use `normalize=True` or `normalize=False`?

**A:** It depends on your data format.

| Your Data Source | Parameter |
|-----------------|-----------|
| Scanpy preprocessed (sc.pp.normalize_total + sc.pp.log1p) | `normalize=False` |
| Seurat object exported from R (normalized) | `normalize=False` |
| 10X raw count matrix | `normalize=True` |
| Downloaded from GEO/SRA (check documentation) | Usually `normalize=True` |
| Your own preprocessing with log-transform | `normalize=False` |

**Golden rule:** If you already did `log1p()` or similar, use `normalize=False`.

---

### Q: Can I use BEER with scaled data (z-scores)?

**A:** No, BEER requires log-normalized counts, not z-scores.

**Why?**
- Z-scores have negative values (mean=0, std=1)
- BEER's variable gene selection needs positive values
- ComBat correction expects non-negative data

**Solution:**
```python
# If you have raw counts saved
adata.X = adata.raw.X.copy()

# Then run BEER with normalization
beer = BEER(n_pcs=50, n_groups=30)
result = beer.fit_transform(adata, batch_key='batch', normalize=True)
```

---

### Q: What if I have more than 2 batches?

**A:** BEER handles multiple batches automatically! Just increase `n_rounds` for better cross-batch pairing.

```python
# For 3-5 batches
beer = BEER(n_pcs=50, n_groups=30, n_rounds=2)

# For 6-10 batches
beer = BEER(n_pcs=50, n_groups=30, n_rounds=3)

# Run as usual
result = beer.fit_transform(adata, batch_key='batch', normalize=False)
```

**How it works:**
- Round 1: Finds initial MN pairs between all batch combinations
- Round 2: Finds additional pairs after masking round 1 pairs
- Round 3: Finds even more pairs for comprehensive coverage

**More rounds = stronger correction**, but diminishing returns after 3 rounds.

---

## Parameter Tuning

### Q: How do I choose the right parameters?

**A:** Start with defaults, then adjust based on results.

#### Default Starting Point
```python
beer = BEER(
    n_pcs=50,              # Good for most datasets
    n_groups=30,           # Works for 1,000-10,000 cells/batch
    n_variable_genes=2000, # Standard for scRNA-seq
    n_rounds=1,            # Fine for 2 batches
    use_combat=True,       # Recommended
    random_seed=123        # For reproducibility
)
```

#### Adjustment Guide

| Parameter | Increase When | Decrease When |
|-----------|--------------|---------------|
| `n_pcs` | Complex dataset, many cell types | Simple dataset, memory issues |
| `n_groups` | Large batches (>10k cells) | Small batches (<1k cells) |
| `n_variable_genes` | Cross-modality integration | Memory limitations |
| `n_rounds` | Multiple batches, weak correction | Over-correction, 2 batches only |

#### Common Scenarios

**Scenario 1: Too few PCs selected (< 10 PCs)**
```python
# Try increasing n_groups
beer = BEER(n_pcs=50, n_groups=50)  # Was 30
```

**Scenario 2: Over-correction (biological signal lost)**
```python
# Decrease n_rounds
beer = BEER(n_pcs=50, n_groups=30, n_rounds=1)  # Was 2

# Or disable ComBat
beer = BEER(n_pcs=50, n_groups=30, use_combat=False)
```

**Scenario 3: Under-correction (batches still separated)**
```python
# Increase n_rounds
beer = BEER(n_pcs=50, n_groups=30, n_rounds=2)  # Was 1

# Or add BBKNN enhancement
from beer import apply_bbknn
umap_bbknn = apply_bbknn(result.adata, result.selected_pcs)
```

---

### Q: What does `n_groups` control?

**A:** `n_groups` determines how finely cells are divided within each batch for MN pair detection.

**How it works:**
1. Within each batch, cells are clustered into `n_groups` using k-means (on UMAP)
2. These groups are used to find mutual nearest neighbors
3. More groups = finer granularity = more potential MN pairs

**Guidelines:**

| Batch Size | Recommended n_groups |
|------------|---------------------|
| < 500 cells | 10-15 |
| 500-2,000 | 20-30 |
| 2,000-10,000 | 30-50 |
| > 10,000 | 50-100 |

**Trade-offs:**
- **Too few groups:** Might merge distinct cell types, lose resolution
- **Too many groups:** Each group has few cells, noisy correlations

**Auto-adjustment:** BEER automatically reduces `n_groups` if a batch has too few cells.

---

### Q: When should I use ComBat correction?

**A:** Use ComBat (`use_combat=True`) in most cases for better batch mixing.

**ComBat is recommended when:**
- ✅ You have moderate to strong batch effects
- ✅ You want both PC selection AND expression correction
- ✅ You have at least 2 batches with > 50 cells each

**Skip ComBat when:**
- ❌ Batches are very unbalanced (e.g., 10,000 vs 100 cells)
- ❌ You only want PC selection, not expression correction
- ❌ ComBat is not installed

**Example:**
```python
# With ComBat (default)
beer = BEER(n_pcs=50, n_groups=30, use_combat=True)
result = beer.fit_transform(adata, batch_key='batch', normalize=False)

# Without ComBat
beer = BEER(n_pcs=50, n_groups=30, use_combat=False)
result = beer.fit_transform(adata, batch_key='batch', normalize=False)
```

---

## Common Errors

### Q: I get "ValueError: operands could not be broadcast together" - what does this mean?

**A:** This error means your data is already normalized, but you're trying to normalize it again.

**Solution:** Use `normalize=False`:

```python
# Wrong (causes error)
result = beer.fit_transform(adata, batch_key='batch', normalize=True)

# Correct
result = beer.fit_transform(adata, batch_key='batch', normalize=False)
```

This error has been **fixed in the latest version**, but the solution is the same.

---

### Q: I get "ValueError: Need at least 2 batches, found 1"

**A:** BEER requires at least 2 batches to remove batch effects.

**Common causes:**
1. Wrong batch key name
2. All cells have the same batch label
3. Batch column has only one unique value

**Check your data:**
```python
# Check batch key
print(adata.obs.columns)  # Is 'batch' in here?

# Check unique batches
print(adata.obs['batch'].unique())  # Should have ≥2 values
print(adata.obs['batch'].value_counts())  # Count per batch
```

**Solution:**
```python
# Use correct batch key
result = beer.fit_transform(adata, batch_key='sample_id')  # Not 'batch'
```

---

### Q: "ImportError: combat-python not installed" - do I need it?

**A:** ComBat is optional but recommended for better batch correction.

**Install it:**
```bash
pip install combat-python
```

**Or disable ComBat:**
```python
beer = BEER(n_pcs=50, n_groups=30, use_combat=False)
result = beer.fit_transform(adata, batch_key='batch', normalize=False)
```

**ComBat vs No ComBat:**
- **With ComBat:** Corrects both PC selection AND gene expression → Better mixing
- **Without ComBat:** Only PC selection → Still effective, just less aggressive

---

### Q: BEER is running very slowly - how can I speed it up?

**A:** Try these optimizations:

**1. Reduce n_pcs (most effective)**
```python
beer = BEER(n_pcs=30, n_groups=30)  # Instead of 50
```

**2. Reduce n_variable_genes**
```python
beer = BEER(n_pcs=50, n_groups=30, n_variable_genes=1500)  # Instead of 2000
```

**3. Reduce n_groups for large datasets**
```python
beer = BEER(n_pcs=50, n_groups=20)  # Instead of 30
```

**4. Disable ComBat**
```python
beer = BEER(n_pcs=50, n_groups=30, use_combat=False)
```

**Benchmarks:**
| Dataset Size | Time (default) | Time (optimized) |
|--------------|----------------|------------------|
| 1,000 cells | 15 sec | 8 sec |
| 10,000 cells | 2 min | 1 min |
| 50,000 cells | 12 min | 6 min |

---

## Results Interpretation

### Q: How many PCs should be selected?

**A:** Typically 20-40 PCs for most datasets, but it varies.

**Guidelines:**

| Selected PCs | Interpretation | Action |
|--------------|----------------|--------|
| 5-15 | Strong batch effects | ✅ Normal, proceed |
| 15-30 | Moderate batch effects | ✅ Good |
| 30-45 | Mild batch effects | ✅ Excellent |
| > 45 | Very little batch effect | ⚠️ Check if batches are real |
| < 5 | Severe batch effects | ⚠️ Adjust parameters |

**If too few PCs (<10):**
1. Increase `n_groups` (e.g., 50 instead of 30)
2. Check if batches are truly different
3. Consider using BBKNN enhancement
4. Adjust correlation thresholds

**If almost all PCs selected (>45):**
1. Batches might not have strong effects
2. Consider whether batch correction is needed
3. Check batch labels are correct

---

### Q: How do I know if BEER worked well?

**A:** Check these indicators:

**1. Visualization (most important)**
```python
import scanpy as sc

sc.pp.neighbors(result.adata, n_pcs=len(result.selected_pcs), use_rep='X_pca')
sc.tl.umap(result.adata)

# Should see batch mixing
sc.pl.umap(result.adata, color='batch')

# Should preserve biology
sc.pl.umap(result.adata, color='cell_type')  # If available
```

**Look for:**
- ✅ Batches are spatially mixed (not separate clusters)
- ✅ Biological cell types still cluster together
- ✅ No artificial clusters at batch boundaries

**2. Correlation plot**
```python
from beer import plot_correlation_scatter
plot_correlation_scatter(result)
```

**Look for:**
- ✅ Clear separation between selected (red) and rejected (gray) PCs
- ✅ Selected PCs cluster in top-right corner (high both correlations)

**3. Quantitative metrics** (if you have cell type labels)
```python
# Silhouette score (higher = better separation)
from sklearn.metrics import silhouette_score

# By cell type (biological - should be high)
bio_score = silhouette_score(result.adata.obsm['X_pca'][:, result.selected_pcs],
                              result.adata.obs['cell_type'])

# By batch (technical - should be low)
batch_score = silhouette_score(result.adata.obsm['X_pca'][:, result.selected_pcs],
                                result.adata.obs['batch'])

print(f"Biology preservation: {bio_score:.3f}")  # Want > 0.3
print(f"Batch mixing: {batch_score:.3f}")        # Want < 0.2
```

---

### Q: Can I manually select which PCs to use?

**A:** Yes! You can either adjust thresholds or manually select PCs.

**Option 1: Adjust thresholds**
```python
beer = BEER(
    n_pcs=50,
    n_groups=30,
    rank_correlation_threshold=0.8,  # Stricter (default: 0.7)
    linear_correlation_threshold=0.8,
    rank_ratio=0.6,  # Top 60% instead of 50%
    linear_ratio=0.6
)
```

**Option 2: Manual selection after fitting**
```python
# Run BEER
result = beer.fit_transform(adata, batch_key='batch', normalize=False)

# Examine correlations
import matplotlib.pyplot as plt
plt.scatter(result.rank_correlations, result.linear_correlations)
plt.xlabel('Rank Correlation')
plt.ylabel('Linear Correlation')
plt.show()

# Manually select PCs
# E.g., select PCs with rank_cor > 0.75 AND linear_cor > 0.75
custom_selected = np.where(
    (result.rank_correlations > 0.75) &
    (result.linear_correlations > 0.75)
)[0]

print(f"Custom selection: {len(custom_selected)} PCs")
print(f"PC indices: {custom_selected}")

# Use custom selection for downstream analysis
sc.pp.neighbors(result.adata, n_pcs=len(custom_selected), use_rep='X_pca')
sc.tl.umap(result.adata)
```

---

### Q: Should I use BBKNN enhancement?

**A:** Use BBKNN if BEER alone doesn't provide sufficient batch mixing.

**When to use BBKNN:**
- ✅ Batches still partially separated after BEER
- ✅ You want even stronger batch correction
- ✅ You're okay with more aggressive correction

**How to use:**
```python
from beer import BEER, apply_bbknn

# Run BEER first
beer = BEER(n_pcs=50, n_groups=30)
result = beer.fit_transform(adata, batch_key='batch', normalize=False)

# Apply BBKNN enhancement
umap_bbknn = apply_bbknn(
    result.adata,
    selected_pcs=result.selected_pcs,
    batch_key='batch',
    neighbors_within_batch=3,  # Adjust for mixing strength
    n_trees=10
)

# Store BBKNN UMAP
result.adata.obsm['X_umap_bbknn'] = umap_bbknn

# Compare
sc.pl.umap(result.adata, color='batch', title='Standard UMAP')
plt.scatter(umap_bbknn[:, 0], umap_bbknn[:, 1],
            c=result.adata.obs['batch'].cat.codes)
plt.title('BBKNN-enhanced UMAP')
```

**Parameters:**
- `neighbors_within_batch`: Lower = stronger mixing (try 2-5)
- `n_trees`: Higher = better accuracy, slower (10 is good default)

---

## Performance Questions

### Q: How much memory does BEER need?

**A:** Memory usage scales with dataset size.

**Estimates:**

| Cells | Genes | Memory (approx) |
|-------|-------|-----------------|
| 1,000 | 2,000 | ~200 MB |
| 10,000 | 2,000 | ~1 GB |
| 50,000 | 2,000 | ~3 GB |
| 100,000 | 2,000 | ~6 GB |

**Reduce memory:**
1. Use sparse matrices (automatic if input is sparse)
2. Reduce `n_variable_genes`
3. Process batches separately, then integrate

---

### Q: Can BEER handle very large datasets (>100k cells)?

**A:** Yes, but you may need to optimize parameters.

**For datasets > 100,000 cells:**

```python
beer = BEER(
    n_pcs=30,              # Reduce from 50
    n_groups=50,           # Increase for large batches
    n_variable_genes=1500, # Reduce if memory limited
    use_combat=False       # Disable for speed
)

result = beer.fit_transform(adata, batch_key='batch', normalize=False)
```

**Benchmark (100,000 cells, 2,000 genes):**
- Default settings: ~25 minutes, ~6 GB RAM
- Optimized settings: ~10 minutes, ~4 GB RAM

---

## Integration with Other Tools

### Q: Can I use BEER with Scanpy?

**A:** Yes! BEER integrates seamlessly with Scanpy.

```python
import scanpy as sc
from beer import BEER

# Standard Scanpy workflow
adata = sc.read_10x_h5('data.h5')
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# BEER batch correction
beer = BEER(n_pcs=50, n_groups=30)
result = beer.fit_transform(adata, batch_key='batch', normalize=False)

# Continue with Scanpy
sc.pp.neighbors(result.adata, n_pcs=len(result.selected_pcs), use_rep='X_pca')
sc.tl.umap(result.adata)
sc.tl.leiden(result.adata)
sc.pl.umap(result.adata, color=['batch', 'leiden'])
```

---




**When to use BEER:**
- You want to understand batch effects
- You need fast correction
- You want control over correction strength

**When to use Harmony:**
- You have many batches (>10)
- You prefer automatic correction
- You don't need to inspect individual PCs

---
