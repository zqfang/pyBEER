"""
Example Usage of BEER Python Implementation
============================================

This script demonstrates how to use the BEER Python package for
batch effect removal in single-cell data.

Author: Feng Zhang (original R), Python port examples
Date: October 2025
"""

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

from beer import (
    BEER,
    BEERResult,
    apply_combat_to_pca,
    apply_bbknn,
    plot_correlation_scatter
)


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
    beer = BEER(
        n_pcs=50,
        n_groups=30,
        n_variable_genes=1000,
        n_rounds=1,
        use_combat=True,
        random_seed=123
    )

    result = beer.fit_transform(adata, batch_key='batch')

    print(f"\nSelected {len(result.selected_pcs)} PCs with low batch correlation")
    print(f"Selected PC indices: {result.selected_pcs}")

    # Visualize results
    plot_correlation_scatter(result, figsize=(10, 6))

    # Compute UMAP with selected PCs
    sc.pp.neighbors(result.adata, n_pcs=len(result.selected_pcs), use_rep='X_pca')
    sc.tl.umap(result.adata)

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    sc.pl.umap(result.adata, color='batch', ax=axes[0], show=False, title='After BEER')
    sc.pl.umap(result.adata, color='beer_group', ax=axes[1], show=False,
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
    # This assumes you have the data; replace with your own data path
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
        beer = BEER(
            n_pcs=50,
            n_groups=30,
            n_variable_genes=2000,
            n_rounds=1,
            use_combat=True,
            random_seed=123
        )

        result = beer.fit_transform(adata, batch_key='batch', normalize=True)

        print(f"\nSelected {len(result.selected_pcs)} PCs")

        # Compute UMAP with selected PCs
        sc.pp.neighbors(result.adata, n_pcs=len(result.selected_pcs), use_rep='X_pca')
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
# Example 3: BBKNN Enhancement
# ===========================================================================

def example3_bbknn_enhancement():
    """Example 3: Using BBKNN enhancement for better batch mixing."""
    print("=" * 80)
    print("Example 3: BBKNN Enhancement")
    print("=" * 80)

    # Create simulated data with strong batch effect
    np.random.seed(123)
    n_cells_per_batch = 300
    n_genes = 1500

    X1 = np.random.negative_binomial(5, 0.3, size=(n_cells_per_batch, n_genes))
    X2 = np.random.negative_binomial(8, 0.4, size=(n_cells_per_batch, n_genes))  # Strong batch effect

    X = np.vstack([X1, X2]).astype(float)
    batch = np.array(['Batch1'] * n_cells_per_batch + ['Batch2'] * n_cells_per_batch)

    adata = ad.AnnData(X)
    adata.obs['batch'] = batch

    # Run BEER
    beer = BEER(n_pcs=50, n_groups=20, n_variable_genes=1000, use_combat=True)
    result = beer.fit_transform(adata, batch_key='batch')

    print(f"Selected {len(result.selected_pcs)} PCs")

    # Apply BBKNN enhancement
    try:
        umap_bbknn = apply_bbknn(
            result.adata,
            selected_pcs=result.selected_pcs,
            batch_key='batch',
            neighbors_within_batch=3,
            n_trees=10,
            n_umap_components=2
        )

        result.adata.obsm['X_umap_bbknn'] = umap_bbknn

        # Compare standard vs BBKNN UMAP
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        sc.pl.umap(result.adata, color='batch', ax=axes[0], show=False,
                   title='Standard UMAP')

        # Plot BBKNN UMAP
        axes[1].scatter(
            umap_bbknn[:, 0],
            umap_bbknn[:, 1],
            c=[0 if b == 'Batch1' else 1 for b in result.adata.obs['batch']],
            cmap='Set1',
            s=5,
            alpha=0.7
        )
        axes[1].set_title('BBKNN-enhanced UMAP')
        axes[1].set_xlabel('UMAP 1')
        axes[1].set_ylabel('UMAP 2')

        plt.tight_layout()
        plt.savefig('example3_bbknn_comparison.png', dpi=300, bbox_inches='tight')
        print("Saved BBKNN comparison to example3_bbknn_comparison.png")

    except ImportError:
        print("BBKNN not installed. Install with: pip install bbknn")

    return result


# ===========================================================================
# Example 4: Multiple Batches
# ===========================================================================

def example4_multiple_batches():
    """Example 4: Handling multiple batches (>2)."""
    print("=" * 80)
    print("Example 4: Multiple Batches")
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
    beer = BEER(
        n_pcs=50,
        n_groups=20,
        n_variable_genes=1000,
        n_rounds=2,  # Use 2 rounds for multiple batches
        use_combat=True
    )

    result = beer.fit_transform(adata, batch_key='batch')

    print(f"\nSelected {len(result.selected_pcs)} PCs")
    print(f"Found {len(result.mutual_pairs)} mutual nearest neighbor pairs")

    # Visualize
    sc.pp.neighbors(result.adata, n_pcs=len(result.selected_pcs), use_rep='X_pca')
    sc.tl.umap(result.adata)

    sc.pl.umap(result.adata, color='batch', title='Multiple Batches after BEER')
    plt.savefig('example4_multiple_batches.png', dpi=300, bbox_inches='tight')
    print("Saved visualization to example4_multiple_batches.png")

    return result


# ===========================================================================
# Example 5: Refitting with Different Parameters
# ===========================================================================

def example5_refitting():
    """Example 5: Re-run BEER with adjusted parameters."""
    print("=" * 80)
    print("Example 5: Parameter Refitting")
    print("=" * 80)

    # Create data
    np.random.seed(123)
    X = np.random.negative_binomial(5, 0.3, size=(500, 1500)).astype(float)
    batch = np.array(['Batch1'] * 250 + ['Batch2'] * 250)

    adata = ad.AnnData(X)
    adata.obs['batch'] = batch

    # Initial run
    beer = BEER(n_pcs=30, n_groups=20, n_rounds=1)
    result1 = beer.fit_transform(adata, batch_key='batch')

    print(f"Initial run: Selected {len(result1.selected_pcs)} PCs")

    # Refit with different parameters
    result2 = beer.refit(n_groups=40, n_rounds=2)

    print(f"After refit: Selected {len(result2.selected_pcs)} PCs")

    # Compare
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    axes[0].scatter(result1.rank_correlations, result1.linear_correlations,
                    c='blue', alpha=0.5, label='Initial')
    axes[0].set_title('Initial Run')
    axes[0].set_xlabel('Rank Correlation')
    axes[0].set_ylabel('Linear Correlation')
    axes[0].grid(True, alpha=0.3)

    axes[1].scatter(result2.rank_correlations, result2.linear_correlations,
                    c='red', alpha=0.5, label='Refit')
    axes[1].set_title('After Refitting')
    axes[1].set_xlabel('Rank Correlation')
    axes[1].set_ylabel('Linear Correlation')
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('example5_refitting_comparison.png', dpi=300, bbox_inches='tight')
    print("Saved comparison to example5_refitting_comparison.png")

    return result1, result2


# ===========================================================================
# Example 6: Custom PC Selection
# ===========================================================================

def example6_custom_pc_selection():
    """Example 6: Manual PC selection based on correlation thresholds."""
    print("=" * 80)
    print("Example 6: Custom PC Selection")
    print("=" * 80)

    # Create data
    np.random.seed(123)
    X = np.random.negative_binomial(5, 0.3, size=(400, 1200)).astype(float)
    batch = np.array(['Batch1'] * 200 + ['Batch2'] * 200)

    adata = ad.AnnData(X)
    adata.obs['batch'] = batch

    # Run BEER with custom thresholds
    beer = BEER(
        n_pcs=50,
        n_groups=30,
        rank_correlation_threshold=0.8,  # Stricter threshold
        linear_correlation_threshold=0.8,
        use_combat=True
    )

    result = beer.fit_transform(adata, batch_key='batch')

    print(f"Selected {len(result.selected_pcs)} PCs with strict threshold (0.8)")

    # Compare with lenient threshold
    beer_lenient = BEER(
        n_pcs=50,
        n_groups=30,
        rank_correlation_threshold=0.5,
        linear_correlation_threshold=0.5
    )

    result_lenient = beer_lenient.fit_transform(adata, batch_key='batch')

    print(f"Selected {len(result_lenient.selected_pcs)} PCs with lenient threshold (0.5)")

    # Visualize difference
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.scatter(result.rank_correlations, result.linear_correlations,
               c='gray', alpha=0.3, label='All PCs')
    ax.scatter(result.rank_correlations[result.selected_pcs],
               result.linear_correlations[result.selected_pcs],
               c='red', s=50, label=f'Strict (n={len(result.selected_pcs)})')
    ax.scatter(result_lenient.rank_correlations[result_lenient.selected_pcs],
               result_lenient.linear_correlations[result_lenient.selected_pcs],
               c='blue', s=30, alpha=0.5, label=f'Lenient (n={len(result_lenient.selected_pcs)})')

    ax.set_xlabel('Rank Correlation')
    ax.set_ylabel('Linear Correlation')
    ax.set_title('PC Selection with Different Thresholds')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.savefig('example6_custom_selection.png', dpi=300, bbox_inches='tight')
    print("Saved comparison to example6_custom_selection.png")

    return result, result_lenient


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
        ("BBKNN Enhancement", example3_bbknn_enhancement),
        ("Multiple Batches", example4_multiple_batches),
        ("Parameter Refitting", example5_refitting),
        ("Custom PC Selection", example6_custom_pc_selection),
    ]

    for name, func in examples:
        try:
            print(f"\n{'='*80}")
            print(f"Running: {name}")
            print(f"{'='*80}\n")
            func()
            print(f"\n✓ {name} completed successfully\n")
        except Exception as e:
            print(f"\n✗ {name} failed: {e}\n")
            import traceback
            traceback.print_exc()

    print("\n" + "=" * 80)
    print("All examples completed!")
    print("=" * 80 + "\n")


if __name__ == '__main__':
    main()
