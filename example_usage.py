"""
Example Usage of BEER Python Implementation
============================================

This script demonstrates how to use the BEER Python package for
batch effect removal in single-cell data.

Author: Feng Zhang (original R), Python port examples
Date: November 2024
"""

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt

from beer import beer, get_use, select_use, BEERResult


# ===========================================================================
# Example 1: Basic Usage with Simulated Data
# ===========================================================================

def example1_basic_usage():
    """Example 1: Basic BEER usage with simulated data."""
    print("=" * 80)
    print("Example 1: Basic Usage")
    print("=" * 80)

    # Create simulated data
    np.random.seed(123)
    n_cells_per_batch = 500
    n_genes = 2000

    # Batch 1
    X1 = np.random.negative_binomial(5, 0.3, size=(n_cells_per_batch, n_genes))
    X1 = X1.astype(float)

    # Batch 2 (with batch effect: shifted distribution)
    X2 = np.random.negative_binomial(5, 0.3, size=(n_cells_per_batch, n_genes))
    X2 = X2 * 1.5  # Batch effect
    X2 = X2.astype(float)

    # Combine
    X = np.vstack([X1, X2])
    batch = np.array(['Batch1'] * n_cells_per_batch + ['Batch2'] * n_cells_per_batch)

    # Create AnnData object
    adata = ad.AnnData(X)
    adata.obs['batch'] = batch
    adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
    adata.obs_names = [f"Cell_{i}" for i in range(len(batch))]

    print(f"Created simulated data: {adata.n_obs} cells, {adata.n_vars} genes")

    # Run BEER
    result = beer(
        adata,
        batch_key='batch',
        gnum=30,
        pcnum=50,
        gn=1000,
        combat=False,
        seed=123,
        n_rounds=1
    )

    print(f"\nSelected {len(result.select)} PCs with low batch correlation")
    print(f"Selected PC indices (0-based): {result.select}")

    # Compute UMAP with selected PCs
    sc.pp.neighbors(result.adata, n_pcs=len(result.select), use_rep='X_pca')
    sc.tl.umap(result.adata)

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    sc.pl.umap(result.adata, color='batch', ax=axes[0], show=False, title='After BEER')
    sc.pl.umap(result.adata, color='group', ax=axes[1], show=False,
               title='BEER Groups', legend_loc='none')

    plt.tight_layout()
    plt.savefig('example1_beer_results.png', dpi=300, bbox_inches='tight')
    print("\nSaved visualization to example1_beer_results.png")

    return result


# ===========================================================================
# Example 2: Real Data Analysis
# ===========================================================================

def example2_real_data():
    """Example 2: BEER with real single-cell data."""
    print("=" * 80)
    print("Example 2: Real Data Analysis")
    print("=" * 80)

    # Download example dataset (PBMC3k from 10X Genomics)
    try:
        adata = sc.datasets.pbmc3k()
        print(f"Loaded PBMC3k dataset: {adata.n_obs} cells, {adata.n_vars} genes")

        # Create artificial batches for demonstration
        n_cells = adata.n_obs
        batch = np.array(['Batch1'] * (n_cells // 2) + ['Batch2'] * (n_cells - n_cells // 2))
        adata.obs['batch'] = batch

        # Basic QC
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)

        print(f"After QC: {adata.n_obs} cells, {adata.n_vars} genes")

        # Run BEER
        result = beer(
            adata,
            batch_key='batch',
            gnum=30,
            pcnum=50,
            gn=2000,
            combat=False,
            seed=123
        )

        print(f"\nSelected {len(result.select)} PCs")

        # Compute UMAP with selected PCs
        sc.pp.neighbors(result.adata, n_pcs=len(result.select), use_rep='X_pca')
        sc.tl.umap(result.adata)

        # Plot
        sc.pl.umap(result.adata, color='batch', title='After BEER Correction')
        plt.savefig('example2_pbmc3k_beer.png', dpi=300, bbox_inches='tight')
        print("Saved visualization to example2_pbmc3k_beer.png")

        return result

    except Exception as e:
        print(f"Could not load PBMC3k dataset: {e}")
        print("Please provide your own data path")
        return None


# ===========================================================================
# Example 3: Multiple Batches
# ===========================================================================

def example3_multiple_batches():
    """Example 3: Handling multiple batches (>2)."""
    print("=" * 80)
    print("Example 3: Multiple Batches")
    print("=" * 80)

    # Create 4 batches
    np.random.seed(123)
    n_cells_per_batch = 200
    n_genes = 1500

    batches_data = []
    batch_labels = []

    for i in range(4):
        # Each batch has slightly different distribution
        X_batch = np.random.negative_binomial(5 + i, 0.3, size=(n_cells_per_batch, n_genes))
        batches_data.append(X_batch)
        batch_labels.extend([f'Batch{i+1}'] * n_cells_per_batch)

    X = np.vstack(batches_data).astype(float)
    batch = np.array(batch_labels)

    adata = ad.AnnData(X)
    adata.obs['batch'] = batch

    print(f"Created data with {len(np.unique(batch))} batches")
    print(f"Batch distribution: {pd.Series(batch).value_counts().to_dict()}")

    # Run BEER
    result = beer(
        adata,
        batch_key='batch',
        gnum=20,
        pcnum=50,
        gn=1000,
        n_rounds=2,  # Use 2 rounds for multiple batches
        combat=False,
        seed=123
    )

    print(f"\nSelected {len(result.select)} PCs")
    print(f"Found {len(result.vp)} mutual nearest neighbor pairs")

    # Visualize
    sc.pp.neighbors(result.adata, n_pcs=len(result.select), use_rep='X_pca')
    sc.tl.umap(result.adata)

    sc.pl.umap(result.adata, color='batch', title='Multiple Batches after BEER')
    plt.savefig('example3_multiple_batches.png', dpi=300, bbox_inches='tight')
    print("Saved visualization to example3_multiple_batches.png")

    return result


# ===========================================================================
# Example 4: Custom PC Selection
# ===========================================================================

def example4_custom_pc_selection():
    """Example 4: Manual PC selection based on correlation thresholds."""
    print("=" * 80)
    print("Example 4: Custom PC Selection")
    print("=" * 80)

    # Create data
    np.random.seed(123)
    X = np.random.negative_binomial(5, 0.3, size=(400, 1200)).astype(float)
    batch = np.array(['Batch1'] * 200 + ['Batch2'] * 200)

    adata = ad.AnnData(X)
    adata.obs['batch'] = batch

    # Run BEER
    result = beer(
        adata,
        batch_key='batch',
        gnum=30,
        pcnum=50,
        gn=1000,
        combat=False,
        seed=123
    )

    print(f"Default selection: {len(result.select)} PCs")

    # Custom selection with strict threshold
    strict_pcs = get_use(result, cut_r=0.8, cut_l=0.8)
    print(f"Strict threshold (0.8): {len(strict_pcs)} PCs")

    # Custom selection with lenient threshold
    lenient_pcs = get_use(result, cut_r=0.5, cut_l=0.5)
    print(f"Lenient threshold (0.5): {len(lenient_pcs)} PCs")

    # Using select_use with rank ratios
    custom_pcs = select_use(result, cut_r=0.7, cut_l=0.7, rr=0.6, rl=0.6)
    print(f"Custom rank ratio (0.6): {len(custom_pcs)} PCs")

    # Visualize correlation scatter
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.scatter(result.cor, result.lcor, c='gray', alpha=0.5, label='All PCs')
    ax.scatter(result.cor[result.select], result.lcor[result.select],
               c='red', s=50, label=f'Selected (n={len(result.select)})')

    ax.axhline(y=0.7, color='blue', linestyle='--', label='Threshold')
    ax.axvline(x=0.7, color='blue', linestyle='--')

    ax.set_xlabel('Spearman Correlation (cor)')
    ax.set_ylabel('Pearson Correlation (lcor)')
    ax.set_title('PC Selection: Batch Correlation')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.savefig('example4_custom_selection.png', dpi=300, bbox_inches='tight')
    print("Saved correlation plot to example4_custom_selection.png")

    return result


# ===========================================================================
# Example 5: Using with Different Flavors of HVG Selection
# ===========================================================================

def example5_with_rmg():
    """Example 5: Remove specific genes from analysis."""
    print("=" * 80)
    print("Example 5: Remove Genes from Analysis")
    print("=" * 80)

    # Create data
    np.random.seed(123)
    n_cells = 400
    n_genes = 1000

    X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(float)
    batch = np.array(['Batch1'] * 200 + ['Batch2'] * 200)

    adata = ad.AnnData(X)
    adata.obs['batch'] = batch
    adata.var_names = [f"Gene_{i}" for i in range(n_genes)]

    # Define genes to remove (e.g., cell cycle genes)
    genes_to_remove = [f"Gene_{i}" for i in range(0, 50)]  # Remove first 50 genes

    print(f"Removing {len(genes_to_remove)} genes from analysis")

    # Run BEER with gene removal
    result = beer(
        adata,
        batch_key='batch',
        gnum=20,
        pcnum=30,
        gn=500,
        rmg=genes_to_remove,
        combat=False,
        seed=123
    )

    print(f"Selected {len(result.select)} PCs")

    # Visualize
    sc.pp.neighbors(result.adata, n_pcs=len(result.select), use_rep='X_pca')
    sc.tl.umap(result.adata)

    sc.pl.umap(result.adata, color='batch', title='BEER with Gene Removal')
    plt.savefig('example5_gene_removal.png', dpi=300, bbox_inches='tight')
    print("Saved visualization to example5_gene_removal.png")

    return result


# ===========================================================================
# Example 6: Accessing Detailed Results
# ===========================================================================

def example6_detailed_results():
    """Example 6: Exploring detailed BEER results."""
    print("=" * 80)
    print("Example 6: Detailed Results")
    print("=" * 80)

    # Create data
    np.random.seed(123)
    X = np.random.negative_binomial(5, 0.3, size=(300, 800)).astype(float)
    batch = np.array(['Batch1'] * 150 + ['Batch2'] * 150)

    adata = ad.AnnData(X)
    adata.obs['batch'] = batch

    # Run BEER
    result = beer(
        adata,
        batch_key='batch',
        gnum=20,
        pcnum=30,
        gn=500,
        seed=123
    )

    # Print detailed results
    print("\n--- BEER Results Summary ---")
    print(f"Number of cells: {result.adata.n_obs}")
    print(f"Number of genes: {result.adata.n_vars}")
    print(f"Number of PCs computed: {result.pcnum}")
    print(f"Number of PCs selected: {len(result.select)}")
    print(f"Number of MN pairs: {len(result.vp)}")

    print("\n--- PC Correlations ---")
    print(f"Spearman correlations (first 10): {result.cor[:10]}")
    print(f"Pearson correlations (first 10): {result.lcor[:10]}")

    print("\n--- P-values and FDR ---")
    print(f"Spearman p-values (first 5): {result.pv[:5]}")
    print(f"FDR-adjusted (first 5): {result.fdr[:5]}")

    print("\n--- Selected PC indices (0-based) ---")
    print(result.select)

    print("\n--- Parameters used ---")
    print(f"GNUM: {result.gnum}")
    print(f"GN: {result.gn}")
    print(f"PCNUM: {result.pcnum}")
    print(f"ROUND: {result.round}")
    print(f"SEED: {result.seed}")
    print(f"COMBAT: {result.combat}")

    # Access AnnData metadata
    print("\n--- AnnData metadata ---")
    print(f"Obs columns: {list(result.adata.obs.columns)}")
    print(f"Layers: {list(result.adata.layers.keys())}")
    print(f"Obsm keys: {list(result.adata.obsm.keys())}")

    return result


# ===========================================================================
# Main Function
# ===========================================================================

def main():
    """Run all examples."""
    print("\n" + "=" * 80)
    print("BEER Python Implementation - Example Usage")
    print("=" * 80 + "\n")

    examples = [
        ("Basic Usage", example1_basic_usage),
        ("Real Data", example2_real_data),
        ("Multiple Batches", example3_multiple_batches),
        ("Custom PC Selection", example4_custom_pc_selection),
        ("Gene Removal", example5_with_rmg),
        ("Detailed Results", example6_detailed_results),
    ]

    for name, func in examples:
        try:
            print(f"\n{'='*80}")
            print(f"Running: {name}")
            print(f"{'='*80}\n")
            func()
            print(f"\n{name} completed successfully\n")
        except Exception as e:
            print(f"\n{name} failed: {e}\n")
            import traceback
            traceback.print_exc()

    print("\n" + "=" * 80)
    print("All examples completed!")
    print("=" * 80 + "\n")


if __name__ == '__main__':
    main()
