"""
BEER: Batch EffEct Remover for Single-Cell Data (Python Implementation)
===========================================================================

Version: 0.2.1
Author: Feng Zhang (original R), Python port
Date: November 2024
License: MIT

Description:
    BEER removes batch effects from single-cell RNA-seq data by detecting
    batch-specific principal components using mutual nearest neighbors
    and selecting biologically relevant PCs for downstream analysis.

Note on matrix orientation:
    - R (Seurat): genes x cells
    - Python (AnnData): cells x genes
    This implementation handles the conversion automatically.

Dependencies:
    - numpy >= 1.20.0
    - scipy >= 1.7.0
    - pandas >= 1.3.0
    - anndata >= 0.8.0
    - scanpy >= 1.9.0
    - scikit-learn >= 1.0.0
    - combat (optional, for ComBat correction)

Usage:
    from beer import beer
    result = beer(adata, batch_key='batch', gnum=30, pcnum=50)
    selected_pcs = result.select
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass
from datetime import datetime
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats
from scipy.stats import rankdata
from sklearn.cluster import KMeans
from statsmodels.stats.multitest import multipletests

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Constants
BEER_VERSION = "0.2.1"
CORMETHOD = "spearman"


@dataclass
class BEERResult:
    """Container for BEER analysis results."""

    adata: ad.AnnData
    vp: np.ndarray  # Valid pairs (MN pairs) - group names
    cor: np.ndarray  # Spearman correlations for each PC
    pv: np.ndarray  # P-values for Spearman correlations
    fdr: np.ndarray  # FDR-adjusted p-values
    lcor: np.ndarray  # Pearson correlations for each PC
    lpv: np.ndarray  # P-values for Pearson correlations
    lc1: np.ndarray  # Linear model intercept p-values
    lc2: np.ndarray  # Linear model slope p-values
    lfdr: np.ndarray  # FDR-adjusted Pearson p-values
    select: np.ndarray  # Selected PC indices (0-based)

    # Parameters
    round: int = 1
    combat: bool = False
    combat_exp: Optional[np.ndarray] = None
    rmg: Optional[list] = None
    gnum: int = 30
    gn: int = 2000
    pcnum: int = 50
    seed: int = 123
    n: int = 2
    app: str = "BEER"


def _to_dense(X) -> np.ndarray:
    """Convert sparse matrix to dense array."""
    if sparse.issparse(X):
        return X.toarray()
    return np.asarray(X)


def _generate_agg(
    data: np.ndarray,
    groups: np.ndarray,
) -> tuple[np.ndarray, list[str]]:
    """
    Aggregate expression by group using sum.

    Args:
        data: Expression matrix (cells x genes) - Python convention
        groups: Group label for each cell

    Returns:
        Aggregated matrix (n_groups x n_genes) and group names
    """
    unique_groups = list(dict.fromkeys(groups))  # Preserve order
    n_groups = len(unique_groups)
    n_genes = data.shape[1]

    agg = np.zeros((n_groups, n_genes))

    for i, group in enumerate(unique_groups):
        mask = groups == group
        n_cells = np.sum(mask)
        if n_cells > 1:
            agg[i, :] = np.sum(data[mask, :], axis=0)
        else:
            agg[i, :] = data[mask, :].flatten()

    return agg, unique_groups


def _get_group(
    embeddings: np.ndarray,
    batch_name: str,
    n_clusters: int,
    random_state: int = 123
) -> np.ndarray:
    """
    Cluster cells using k-means and create group labels.

    Args:
        embeddings: Cell embeddings (n_cells x n_dims)
        batch_name: Name of the batch
        n_clusters: Number of clusters
        random_state: Random seed

    Returns:
        Group labels for each cell (batch_cluster format)
    """
    logger.info(f"Get group for: {batch_name}")

    kmeans = KMeans(
        n_clusters=n_clusters,
        max_iter=100,
        random_state=random_state,
        n_init=10
    )
    clusters = kmeans.fit_predict(embeddings)

    # Create group labels: batch_cluster
    groups = np.array([f"{batch_name}_{c}" for c in clusters])

    logger.info(f"Group Number: {len(np.unique(groups))}")

    return groups


def _get_mn_pairs(
    cor_matrix: np.ndarray,
    row_names: list,
    col_names: list
) -> list[tuple[str, str]]:
    """
    Find mutual nearest neighbor pairs from correlation matrix.

    A pair (i, j) is MN if cor[i,j] is max in both row i and column j.
    """
    pairs = []

    for i in range(cor_matrix.shape[0]):
        for j in range(cor_matrix.shape[1]):
            cor_val = cor_matrix[i, j]
            if cor_val != -99999:
                # Check mutual maximum
                if (cor_val == np.max(cor_matrix[i, :]) and
                    cor_val == np.max(cor_matrix[:, j])):
                    pairs.append((row_names[i], col_names[j]))

    return pairs


def _get_vp_all(
    adata: ad.AnnData,
    n_rounds: int = 1
) -> np.ndarray:
    """
    Find all valid mutual nearest neighbor pairs across batches.

    Args:
        adata: AnnData with 'batch' and 'group' in obs,
               and log-normalized data in layers['log_normalized']
        n_rounds: Number of rounds for finding pairs

    Returns:
        Array of valid pairs (n_pairs x 2) containing group names
    """
    logger.info("Finding MN pairs...")

    # Get normalized data (cells x genes)
    data = _to_dense(adata.layers['log_normalized'])
    groups = adata.obs['group'].values

    # Aggregate by group (n_groups x n_genes)
    agg_matrix, group_names = _generate_agg(data, groups)

    # Compute correlation between groups
    # Spearman: rank the genes for each group, then correlate
    if CORMETHOD == "spearman":
        # Rank along genes (axis=1) for each group
        ranked = np.apply_along_axis(rankdata, 1, agg_matrix)
        cor_matrix = np.corrcoef(ranked)
    else:
        cor_matrix = np.corrcoef(agg_matrix)

    # Get unique batches
    unique_batches = list(dict.fromkeys(adata.obs['batch'].values))

    # Get batch for each group
    def get_batch(group_name):
        return group_name.rsplit('_', 1)[0]

    group_batch = np.array([get_batch(g) for g in group_names])

    all_pairs = []

    for round_idx in range(n_rounds):
        # Mask already found pairs
        if all_pairs:
            for p1, p2 in all_pairs:
                b1 = get_batch(p1)
                b2 = get_batch(p2)

                b1_idx = np.where(group_batch == b1)[0]
                b2_idx = np.where(group_batch == b2)[0]

                p1_idx = group_names.index(p1)
                p2_idx = group_names.index(p2)

                # Mask cross-batch correlations
                cor_matrix[p1_idx, b2_idx] = -99999
                cor_matrix[p2_idx, b1_idx] = -99999
                cor_matrix[b2_idx, p1_idx] = -99999
                cor_matrix[b1_idx, p2_idx] = -99999

        # Find pairs between all batch combinations
        for i in range(len(unique_batches) - 1):
            for j in range(i + 1, len(unique_batches)):
                b1 = unique_batches[i]
                b2 = unique_batches[j]

                b1_idx = np.where(group_batch == b1)[0]
                b2_idx = np.where(group_batch == b2)[0]

                # Extract submatrix
                sub_matrix = cor_matrix[np.ix_(b1_idx, b2_idx)]
                row_names = [group_names[k] for k in b1_idx]
                col_names = [group_names[k] for k in b2_idx]

                pairs = _get_mn_pairs(sub_matrix, row_names, col_names)
                all_pairs.extend(pairs)

        logger.info(f"ROUND: {round_idx + 1}")

    # Remove duplicates while preserving order
    seen = set()
    unique_pairs = []
    for pair in all_pairs:
        pair_key = tuple(sorted(pair))
        if pair_key not in seen:
            seen.add(pair_key)
            unique_pairs.append(pair)

    logger.info(f"Number of MN pairs: {len(unique_pairs)}")

    return np.array(unique_pairs) if unique_pairs else np.array([]).reshape(0, 2)


def _evaluate_pro_beer(
    pca_embeddings: np.ndarray,
    groups: np.ndarray,
    valid_pairs: np.ndarray
) -> dict:
    """
    Evaluate each PC for batch correlation using MN pair quantiles.

    Args:
        pca_embeddings: PCA embeddings (n_cells x n_pcs)
        groups: Group labels for each cell
        valid_pairs: Array of valid pairs (n_pairs x 2) with group names

    Returns:
        Dictionary with correlation statistics for each PC
    """
    logger.info("Evaluating PCs...")
    logger.info("Start")

    n_pcs = pca_embeddings.shape[1]

    all_cor = []
    all_pv = []
    all_lcor = []
    all_lpv = []
    all_lc1 = []
    all_lc2 = []

    for pc_idx in range(n_pcs):
        pc_values = pca_embeddings[:, pc_idx]

        lst1_quantile = []
        lst2_quantile = []

        for pair in valid_pairs:
            p1, p2 = pair[0], pair[1]

            # Get cell indices for each group
            idx1 = np.where(groups == p1)[0]
            idx2 = np.where(groups == p2)[0]

            # Get quantiles (0%, 25%, 50%, 75%, 100%)
            q1 = np.percentile(pc_values[idx1], [0, 25, 50, 75, 100])
            q2 = np.percentile(pc_values[idx2], [0, 25, 50, 75, 100])

            lst1_quantile.extend(q1)
            lst2_quantile.extend(q2)

        lst1_quantile = np.array(lst1_quantile)
        lst2_quantile = np.array(lst2_quantile)

        # Spearman correlation
        spearman_result = stats.spearmanr(lst1_quantile, lst2_quantile)
        all_cor.append(spearman_result.correlation)
        all_pv.append(spearman_result.pvalue)

        # Pearson correlation
        pearson_result = stats.pearsonr(lst1_quantile, lst2_quantile)
        all_lcor.append(pearson_result[0])
        all_lpv.append(pearson_result[1])

        # Linear regression for intercept/slope p-values
        slope, intercept, r_value, p_slope, std_err = stats.linregress(
            lst2_quantile, lst1_quantile
        )

        # Calculate intercept p-value
        n = len(lst1_quantile)
        if n > 2 and std_err > 0:
            ss_x = np.sum((lst2_quantile - np.mean(lst2_quantile))**2)
            if ss_x > 0:
                se_intercept = std_err * np.sqrt(
                    np.sum(lst2_quantile**2) / (n * ss_x)
                )
                if se_intercept > 0:
                    t_intercept = intercept / se_intercept
                    p_intercept = 2 * (1 - stats.t.cdf(abs(t_intercept), n - 2))
                else:
                    p_intercept = 1.0
            else:
                p_intercept = 1.0
        else:
            p_intercept = 1.0

        all_lc1.append(p_intercept)
        all_lc2.append(p_slope)

        logger.info(f"{pc_idx + 1}")

    # FDR correction
    _, fdr, _, _ = multipletests(all_pv, method='fdr_bh')
    _, lfdr, _, _ = multipletests(all_lpv, method='fdr_bh')

    logger.info("Finished!!!")

    return {
        'cor': np.array(all_cor),
        'pv': np.array(all_pv),
        'fdr': fdr,
        'lcor': np.array(all_lcor),
        'lpv': np.array(all_lpv),
        'lc1': np.array(all_lc1),
        'lc2': np.array(all_lc2),
        'lfdr': lfdr
    }


def beer(
    adata: ad.AnnData,
    batch_key: str = 'batch',
    gnum: int = 30,
    pcnum: int = 50,
    gn: int = 2000,
    combat: bool = False,
    seed: int = 123,
    n_components: int = 2,
    n_rounds: int = 1,
    rmg: Optional[list] = None,
) -> BEERResult:
    """
    BEER: Batch Effect Remover for single-cell data.

    This function implements the BEER algorithm for batch effect removal
    by identifying and selecting PCs that are not batch-specific.

    Args:
        adata: AnnData object (cells x genes)
        batch_key: Key in adata.obs containing batch labels
        gnum: Number of groups per batch for MNN finding
        pcnum: Number of PCs to compute
        gn: Number of highly variable genes per batch
        combat: Whether to apply ComBat batch correction (requires pycombat)
        seed: Random seed for reproducibility
        n_components: Number of UMAP components
        n_rounds: Number of rounds for finding MN pairs
        rmg: List of genes to remove from analysis

    Returns:
        BEERResult object containing:
        - adata: Processed AnnData with PCA/UMAP
        - vp: Valid MN pairs
        - cor/lcor: Spearman/Pearson correlations per PC
        - select: Indices of selected PCs (0-based)

    Example:
        >>> import scanpy as sc
        >>> adata = sc.read_h5ad('data.h5ad')
        >>> result = beer(adata, batch_key='batch', gnum=30, pcnum=50)
        >>> selected_pcs = result.select
        >>> # Use selected PCs for downstream analysis
        >>> sc.pp.neighbors(adata, n_pcs=len(selected_pcs), use_rep='X_pca')
    """
    np.random.seed(seed)

    logger.info("BEER start!")
    logger.info(datetime.now())

    # Validate input
    if batch_key not in adata.obs.columns:
        raise ValueError(f"Batch key '{batch_key}' not found in adata.obs")

    # Copy data
    adata = adata.copy()

    # Process batch labels (replace _ with .)
    batch = adata.obs[batch_key].astype(str).str.replace('_', '.', regex=False).values
    adata.obs['batch'] = batch
    unique_batches = list(dict.fromkeys(batch))

    logger.info(f"Group number (GNUM) is: {gnum}")
    logger.info(f"Variable gene number (GN) of each batch is: {gn}")
    logger.info(f"ROUND is: {n_rounds}")

    # Find highly variable genes per batch
    all_hvg = []
    for i, this_batch in enumerate(unique_batches):
        logger.info(f"{i + 1}")
        logger.info(this_batch)

        batch_mask = batch == this_batch
        batch_adata = adata[batch_mask].copy()

        # Normalize
        sc.pp.normalize_total(batch_adata, target_sum=1e4)
        sc.pp.log1p(batch_adata)

        # Find HVGs
        try:
            sc.pp.highly_variable_genes(
                batch_adata,
                n_top_genes=min(gn, batch_adata.n_vars),
                flavor='seurat_v3',
                span=0.3
            )
        except Exception:
            # Fallback to seurat flavor if seurat_v3 fails
            sc.pp.highly_variable_genes(
                batch_adata,
                n_top_genes=min(gn, batch_adata.n_vars),
                flavor='seurat'
            )

        hvg = batch_adata.var_names[batch_adata.var['highly_variable']].tolist()
        all_hvg.extend(hvg)

    # Get unique HVGs
    var_genes = list(dict.fromkeys(all_hvg))

    logger.info(f"Total variable gene number (GN) is: {len(var_genes)}")

    # Remove specified genes
    if rmg is not None:
        logger.info(f"Total removed gene number is: {len(rmg)}")
        var_genes = [g for g in var_genes if g not in rmg]
        logger.info(f"Total used gene number is: {len(var_genes)}")

    # Mark variable genes
    adata.var['highly_variable'] = adata.var_names.isin(var_genes)

    # Normalize full dataset
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Store log-normalized data
    adata.layers['log_normalized'] = adata.X.copy()

    combat_exp = None

    if combat:
        logger.info("Applying ComBat correction...")
        try:
            from combat.pycombat import pycombat

            # Get expression for variable genes (cells x genes)
            var_idx = adata.var_names.isin(var_genes)
            X_var = _to_dense(adata.X[:, var_idx])

            # pycombat expects genes x samples
            combat_data = pycombat(X_var.T, batch)

            # Clean up
            combat_data = np.clip(combat_data.T, 0, None)  # Back to cells x genes
            combat_data = np.nan_to_num(combat_data)

            combat_exp = combat_data

            # Create temp AnnData for scaling and PCA
            temp_adata = ad.AnnData(
                X=combat_data,
                obs=adata.obs.copy()
            )
            temp_adata.var_names = adata.var_names[var_idx]

            # Scale
            sc.pp.scale(temp_adata)

            # PCA
            logger.info("Calculating PCs...")
            sc.tl.pca(temp_adata, n_comps=pcnum, random_state=seed)

            # Copy results
            adata.obsm['X_pca'] = temp_adata.obsm['X_pca']
            adata.uns['pca'] = temp_adata.uns['pca']

        except ImportError:
            logger.warning("pycombat not available, skipping ComBat")
            combat = False

    if not combat:
        # Standard processing
        adata_hvg = adata[:, var_genes].copy()
        sc.pp.scale(adata_hvg)

        logger.info("Calculating PCs...")
        sc.tl.pca(adata_hvg, n_comps=pcnum, random_state=seed)

        adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
        adata.uns['pca'] = adata_hvg.uns['pca']

    # UMAP
    sc.pp.neighbors(adata, n_pcs=pcnum, random_state=seed)
    sc.tl.umap(adata, n_components=n_components, random_state=seed)

    # Cluster cells within each batch
    umap_embeddings = adata.obsm['X_umap']
    groups = np.empty(adata.n_obs, dtype=object)

    for this_batch in unique_batches:
        batch_mask = batch == this_batch
        batch_idx = np.where(batch_mask)[0]
        batch_umap = umap_embeddings[batch_mask]

        this_gnum = min(gnum, len(batch_idx) - 1)
        batch_groups = _get_group(batch_umap, this_batch, this_gnum, seed)
        groups[batch_idx] = batch_groups

    adata.obs['group'] = groups

    # Find MN pairs
    vp = _get_vp_all(adata, n_rounds)

    # Get PCA embeddings
    pca_embeddings = adata.obsm['X_pca']

    # Mark cells in valid pairs
    map_labels = np.full(adata.n_obs, 'NA', dtype=object)
    if len(vp) > 0:
        map_labels[np.isin(groups, vp[:, 0])] = 'V1'
        map_labels[np.isin(groups, vp[:, 1])] = 'V2'
    adata.obs['map'] = map_labels

    # Evaluate PCs
    if len(vp) == 0:
        logger.warning("No MN pairs found! Using all PCs.")
        n_pcs = pca_embeddings.shape[1]
        eval_results = {
            'cor': np.ones(n_pcs),
            'pv': np.zeros(n_pcs),
            'fdr': np.zeros(n_pcs),
            'lcor': np.ones(n_pcs),
            'lpv': np.zeros(n_pcs),
            'lc1': np.ones(n_pcs),
            'lc2': np.zeros(n_pcs),
            'lfdr': np.zeros(n_pcs)
        }
        selected = np.arange(n_pcs)
    else:
        eval_results = _evaluate_pro_beer(pca_embeddings, groups, vp)

        # Select PCs based on correlations
        n_pcs = len(eval_results['cor'])
        cor_ranks = rankdata(eval_results['cor'])
        lcor_ranks = rankdata(eval_results['lcor'])

        selected = np.where(
            ((cor_ranks >= n_pcs / 2) | (eval_results['cor'] > 0.7)) &
            ((lcor_ranks >= n_pcs / 2) | (eval_results['lcor'] > 0.7))
        )[0]

    logger.info("############################################################################")
    logger.info("BEER cheers !!! All main steps finished.")
    logger.info("############################################################################")
    logger.info(datetime.now())

    return BEERResult(
        adata=adata,
        vp=vp,
        cor=eval_results['cor'],
        pv=eval_results['pv'],
        fdr=eval_results['fdr'],
        lcor=eval_results['lcor'],
        lpv=eval_results['lpv'],
        lc1=eval_results['lc1'],
        lc2=eval_results['lc2'],
        lfdr=eval_results['lfdr'],
        select=selected,
        round=n_rounds,
        combat=combat,
        combat_exp=combat_exp,
        rmg=rmg,
        gnum=gnum,
        gn=gn,
        pcnum=pcnum,
        seed=seed,
        n=n_components,
        app='BEER'
    )


def get_use(result: BEERResult, cut_r: float = 0.7, cut_l: float = 0.7) -> np.ndarray:
    """
    Get PC indices based on correlation cutoffs.

    Args:
        result: BEERResult object
        cut_r: Cutoff for Spearman correlation
        cut_l: Cutoff for Pearson correlation

    Returns:
        Array of selected PC indices (0-based)
    """
    return np.where((result.cor > cut_r) & (result.lcor > cut_l))[0]


def select_use(
    result: BEERResult,
    cut_r: float = 0.7,
    cut_l: float = 0.7,
    rr: float = 0.5,
    rl: float = 0.5
) -> np.ndarray:
    """
    Select PC indices based on rank and correlation cutoffs.

    Args:
        result: BEERResult object
        cut_r: Cutoff for Spearman correlation
        cut_l: Cutoff for Pearson correlation
        rr: Rank ratio for Spearman
        rl: Rank ratio for Pearson

    Returns:
        Array of selected PC indices (0-based)
    """
    n_pcs = len(result.cor)
    cor_ranks = rankdata(result.cor)
    lcor_ranks = rankdata(result.lcor)

    return np.where(
        ((cor_ranks >= n_pcs * rr) | (result.cor > cut_r)) &
        ((lcor_ranks >= n_pcs * rl) | (result.lcor > cut_l))
    )[0]


# Alias for compatibility
MBEER = beer


if __name__ == "__main__":
    print(f"BEER v{BEER_VERSION} - Batch Effect Remover for single-cell data")
    print("")
    print("Usage:")
    print("  from beer import beer")
    print("  result = beer(adata, batch_key='batch', gnum=30, pcnum=50)")
    print("  selected_pcs = result.select")
