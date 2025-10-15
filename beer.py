"""
BEER: Batch EffEct Remover for Single-Cell Data (Python Implementation)
===========================================================================

Version: 0.2.0
Author: Feng Zhang (original R), Python port by Claude
Date: October 2025
License: MIT

Description:
    BEER removes batch effects from single-cell RNA-seq and ATAC-seq data
    by detecting batch-specific principal components using mutual nearest
    neighbors and selecting biologically relevant PCs for downstream analysis.

Reference:
    Zhang, F., Wu, Y., & Tian, W. (2019). A novel approach to remove the
    batch effect of single-cell data. Cell Discovery, 5, 46.
    https://doi.org/10.1038/s41421-019-0114-x

Dependencies:
    - numpy >= 1.20.0
    - scipy >= 1.7.0
    - pandas >= 1.3.0
    - anndata >= 0.8.0
    - scanpy >= 1.9.0
    - scikit-learn >= 1.0.0
    - bbknn >= 1.5.0 (optional)
    - combat-python >= 0.3.0 (for ComBat correction)

Usage:
    from beer_python import BEER

    # Create BEER instance
    beer = BEER(n_pcs=50, n_groups=30, n_variable_genes=2000, random_seed=123)

    # Fit and transform data
    result = beer.fit_transform(adata, batch_key='batch')

    # Access results
    selected_pcs = result['selected_pcs']
    corrected_adata = result['adata']
"""

import warnings
from typing import Optional, Union, List, Tuple, Dict, Any
from dataclasses import dataclass
import logging

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.spatial.distance import cdist
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import anndata as ad

# Optional dependencies
try:
    import scanpy as sc
    HAS_SCANPY = True
except ImportError:
    HAS_SCANPY = False
    warnings.warn("scanpy not installed. Some features may be unavailable.")

try:
    import bbknn
    HAS_BBKNN = True
except ImportError:
    HAS_BBKNN = False
    warnings.warn("bbknn not installed. BBKNN enhancement unavailable.")

try:
    from combat.pycombat import pycombat
    HAS_COMBAT = True
except ImportError:
    HAS_COMBAT = False
    warnings.warn("combat-python not installed. ComBat correction unavailable.")


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ============================================================================
# Constants and Configuration
# ============================================================================

BEER_VERSION = "0.2.0"
DEFAULT_CORRELATION_METHOD = "spearman"


# ============================================================================
# Data Classes for Type Safety
# ============================================================================

@dataclass
class BEERConfig:
    """Configuration parameters for BEER algorithm."""
    n_pcs: int = 50
    n_groups: int = 30
    n_variable_genes: int = 2000
    n_rounds: int = 1
    correlation_method: str = DEFAULT_CORRELATION_METHOD
    use_combat: bool = True
    random_seed: int = 123
    rank_correlation_threshold: float = 0.7
    linear_correlation_threshold: float = 0.7
    rank_ratio: float = 0.5
    linear_ratio: float = 0.5


@dataclass
class BEERResult:
    """Container for BEER analysis results."""
    adata: ad.AnnData
    selected_pcs: np.ndarray
    rank_correlations: np.ndarray
    linear_correlations: np.ndarray
    p_values: np.ndarray
    fdr_values: np.ndarray
    mutual_pairs: np.ndarray
    group_labels: np.ndarray
    config: BEERConfig


# ============================================================================
# Utility Functions
# ============================================================================

def check_anndata_format(adata: ad.AnnData, batch_key: str) -> None:
    """
    Validate AnnData object format and batch key.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    batch_key : str
        Key in adata.obs containing batch labels.

    Raises
    ------
    ValueError
        If batch_key not found or data format is invalid.
    """
    if batch_key not in adata.obs.columns:
        raise ValueError(f"Batch key '{batch_key}' not found in adata.obs")

    if adata.X is None:
        raise ValueError("adata.X is None. Provide expression matrix.")

    n_batches = adata.obs[batch_key].nunique()
    if n_batches < 2:
        raise ValueError(f"Need at least 2 batches, found {n_batches}")

    logger.info(f"Validated AnnData: {adata.n_obs} cells, {adata.n_vars} genes, {n_batches} batches")


def to_dense(X: Union[np.ndarray, sparse.spmatrix]) -> np.ndarray:
    """
    Convert sparse matrix to dense numpy array efficiently.

    Parameters
    ----------
    X : array-like or sparse matrix
        Input matrix.

    Returns
    -------
    np.ndarray
        Dense numpy array.
    """
    if sparse.issparse(X):
        return X.toarray()
    return np.asarray(X)


def normalize_counts(
    X: Union[np.ndarray, sparse.spmatrix],
    target_sum: float = 1e4,
    log_transform: bool = True,
    axis: int = 1
) -> np.ndarray:
    """
    Normalize count matrix (library size normalization + log transform).

    Parameters
    ----------
    X : array-like
        Count matrix. If axis=1, assumes (cells x genes). If axis=0, assumes (genes x cells).
    target_sum : float
        Target sum for normalization (default: 10,000).
    log_transform : bool
        Apply log1p transformation after normalization.
    axis : int
        Axis along which to normalize.
        - axis=1: Normalize cells (rows), expects (cells x genes) [default]
        - axis=0: Normalize cells (columns), expects (genes x cells)

    Returns
    -------
    np.ndarray
        Normalized expression matrix (same shape as input).
    """
    X_dense = to_dense(X)

    # Normalize each cell to target_sum
    cell_sums = X_dense.sum(axis=axis, keepdims=True)
    cell_sums[cell_sums == 0] = 1  # Avoid division by zero
    X_norm = (X_dense / cell_sums) * target_sum

    if log_transform:
        X_norm = np.log1p(X_norm)

    return X_norm


def select_variable_genes(
    X: np.ndarray,
    n_top_genes: int = 2000,
    flavor: str = 'seurat'
) -> np.ndarray:
    """
    Select highly variable genes.

    Parameters
    ----------
    X : np.ndarray
        Normalized expression matrix (genes x cells).
    n_top_genes : int
        Number of top variable genes to select.
    flavor : str
        Method for variable gene selection ('seurat' or 'cell_ranger').

    Returns
    -------
    np.ndarray
        Boolean mask of selected genes.
    """
    # Calculate mean and variance
    mean = np.mean(X, axis=1)
    var = np.var(X, axis=1, ddof=1)

    # Avoid division by zero
    mean = np.clip(mean, 1e-10, None)

    if flavor == 'seurat':
        # Seurat v3 method: variance/mean ratio
        dispersion = var / mean
    elif flavor == 'cell_ranger':
        # Cell Ranger method: coefficient of variation
        dispersion = np.sqrt(var) / mean
    else:
        raise ValueError(f"Unknown flavor: {flavor}")

    # Select top genes by dispersion
    top_idx = np.argsort(dispersion)[::-1][:n_top_genes]
    mask = np.zeros(len(X), dtype=bool)
    mask[top_idx] = True

    return mask


def aggregate_by_group(
    X: np.ndarray,
    groups: np.ndarray,
    method: str = 'sum'
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Aggregate expression matrix by group labels.

    Parameters
    ----------
    X : np.ndarray
        Expression matrix (genes x cells).
    groups : np.ndarray
        Group labels for each cell.
    method : str
        Aggregation method ('sum' or 'mean').

    Returns
    -------
    X_agg : np.ndarray
        Aggregated expression matrix (genes x groups).
    group_names : np.ndarray
        Unique group names.
    """
    unique_groups = np.unique(groups)
    n_genes = X.shape[0]
    n_groups = len(unique_groups)

    X_agg = np.zeros((n_genes, n_groups))

    for i, group in enumerate(unique_groups):
        mask = groups == group
        if method == 'sum':
            X_agg[:, i] = np.sum(X[:, mask], axis=1)
        elif method == 'mean':
            X_agg[:, i] = np.mean(X[:, mask], axis=1)
        else:
            raise ValueError(f"Unknown aggregation method: {method}")

    return X_agg, unique_groups


# ============================================================================
# Core BEER Algorithm Components
# ============================================================================

class MutualNearestNeighbors:
    """Find mutual nearest neighbor (MN) pairs between batches."""

    def __init__(self, correlation_method: str = 'spearman'):
        """
        Initialize MNN finder.

        Parameters
        ----------
        correlation_method : str
            Correlation method ('spearman' or 'pearson').
        """
        self.correlation_method = correlation_method

    def _compute_correlation(self, X: np.ndarray) -> np.ndarray:
        """Compute correlation matrix."""
        if self.correlation_method == 'spearman':
            # Use rank correlation
            X_ranked = np.apply_along_axis(stats.rankdata, 0, X)
            corr = np.corrcoef(X_ranked.T)
        elif self.correlation_method == 'pearson':
            corr = np.corrcoef(X.T)
        else:
            raise ValueError(f"Unknown correlation method: {self.correlation_method}")

        return corr

    def _extract_batch_from_group(self, group_name: str) -> str:
        """Extract batch name from group identifier."""
        return group_name.split('_')[0]

    def _find_mn_pairs_in_matrix(self, corr_mat: np.ndarray) -> List[Tuple[int, int]]:
        """
        Find mutual nearest neighbors in correlation matrix.

        A pair (i, j) is MN if:
        - corr[i, j] is maximum in row i
        - corr[i, j] is maximum in column j
        """
        pairs = []
        n_rows, n_cols = corr_mat.shape

        for i in range(n_rows):
            for j in range(n_cols):
                if corr_mat[i, j] == -99999:
                    continue

                # Check if mutual maximum
                if (corr_mat[i, j] == np.max(corr_mat[i, :]) and
                    corr_mat[i, j] == np.max(corr_mat[:, j])):
                    pairs.append((i, j))

        return pairs

    def find_pairs(
        self,
        X: np.ndarray,
        group_labels: np.ndarray,
        batch_labels: np.ndarray,
        n_rounds: int = 1
    ) -> np.ndarray:
        """
        Find mutual nearest neighbor pairs across batches.

        Parameters
        ----------
        X : np.ndarray
            Aggregated expression matrix (genes x groups).
        group_labels : np.ndarray
            Group names.
        batch_labels : np.ndarray
            Batch labels for each group.
        n_rounds : int
            Number of MN detection rounds.

        Returns
        -------
        np.ndarray
            MN pairs (n_pairs x 2), each row contains indices of paired groups.
        """
        logger.info(f"Finding mutual nearest neighbors ({n_rounds} rounds)...")

        # Compute correlation matrix
        corr_mat = self._compute_correlation(X)
        unique_batches = np.unique(batch_labels)

        all_pairs = []

        for round_idx in range(n_rounds):
            logger.info(f"MN detection round {round_idx + 1}/{n_rounds}")

            # Mask already found pairs
            if len(all_pairs) > 0:
                for i, j in all_pairs:
                    batch_i = batch_labels[i]
                    batch_j = batch_labels[j]

                    # Mask cross-batch correlations for found pairs
                    mask_i = batch_labels == batch_i
                    mask_j = batch_labels == batch_j

                    corr_mat[i, mask_j] = -99999
                    corr_mat[j, mask_i] = -99999
                    corr_mat[mask_j, i] = -99999
                    corr_mat[mask_i, j] = -99999

            # Find pairs between all batch combinations
            for b1_idx in range(len(unique_batches) - 1):
                for b2_idx in range(b1_idx + 1, len(unique_batches)):
                    b1 = unique_batches[b1_idx]
                    b2 = unique_batches[b2_idx]

                    mask_b1 = batch_labels == b1
                    mask_b2 = batch_labels == b2

                    # Extract submatrix
                    sub_corr = corr_mat[np.ix_(mask_b1, mask_b2)]

                    # Find MN pairs in submatrix
                    sub_pairs = self._find_mn_pairs_in_matrix(sub_corr)

                    # Convert back to global indices
                    idx_b1 = np.where(mask_b1)[0]
                    idx_b2 = np.where(mask_b2)[0]

                    for i, j in sub_pairs:
                        all_pairs.append((idx_b1[i], idx_b2[j]))

        # Remove duplicates
        all_pairs = list(set(all_pairs))

        logger.info(f"Found {len(all_pairs)} mutual nearest neighbor pairs")

        return np.array(all_pairs)


class PCEvaluator:
    """Evaluate principal components for batch correlation."""

    def __init__(self, correlation_method: str = 'spearman'):
        """
        Initialize PC evaluator.

        Parameters
        ----------
        correlation_method : str
            Correlation method ('spearman' or 'pearson').
        """
        self.correlation_method = correlation_method

    def evaluate(
        self,
        pca_embeddings: np.ndarray,
        group_labels: np.ndarray,
        mn_pairs: np.ndarray
    ) -> Dict[str, np.ndarray]:
        """
        Evaluate each PC for batch correlation using MN pairs.

        Parameters
        ----------
        pca_embeddings : np.ndarray
            PCA embeddings (n_cells x n_pcs).
        group_labels : np.ndarray
            Group labels for each cell.
        mn_pairs : np.ndarray
            Mutual nearest neighbor pairs (n_pairs x 2).

        Returns
        -------
        dict
            Dictionary containing correlation statistics for each PC.
        """
        logger.info("Evaluating principal components for batch effects...")

        n_pcs = pca_embeddings.shape[1]
        n_pairs = len(mn_pairs)

        rank_correlations = np.zeros(n_pcs)
        linear_correlations = np.zeros(n_pcs)
        rank_pvalues = np.zeros(n_pcs)
        linear_pvalues = np.zeros(n_pcs)
        linear_coef_pvalues_intercept = np.zeros(n_pcs)
        linear_coef_pvalues_slope = np.zeros(n_pcs)

        unique_groups = np.unique(group_labels)

        for pc_idx in range(n_pcs):
            # Extract quantiles for each group in this PC
            lst1_quantiles = []
            lst2_quantiles = []

            for pair_idx in range(n_pairs):
                group1_idx = mn_pairs[pair_idx, 0]
                group2_idx = mn_pairs[pair_idx, 1]

                group1 = unique_groups[group1_idx]
                group2 = unique_groups[group2_idx]

                # Get cells in each group
                cells1 = group_labels == group1
                cells2 = group_labels == group2

                # Calculate quantiles (5 quantiles: 0%, 25%, 50%, 75%, 100%)
                q1 = np.percentile(pca_embeddings[cells1, pc_idx], [0, 25, 50, 75, 100])
                q2 = np.percentile(pca_embeddings[cells2, pc_idx], [0, 25, 50, 75, 100])

                lst1_quantiles.extend(q1)
                lst2_quantiles.extend(q2)

            lst1_quantiles = np.array(lst1_quantiles)
            lst2_quantiles = np.array(lst2_quantiles)

            # Rank correlation (Spearman)
            rank_cor, rank_pv = stats.spearmanr(lst1_quantiles, lst2_quantiles)
            rank_correlations[pc_idx] = rank_cor
            rank_pvalues[pc_idx] = rank_pv

            # Linear correlation (Pearson)
            linear_cor, linear_pv = stats.pearsonr(lst1_quantiles, lst2_quantiles)
            linear_correlations[pc_idx] = linear_cor
            linear_pvalues[pc_idx] = linear_pv

            # Linear regression p-values
            from scipy.stats import linregress
            slope, intercept, r_value, p_value, std_err = linregress(lst2_quantiles, lst1_quantiles)
            linear_coef_pvalues_slope[pc_idx] = p_value

            if (pc_idx + 1) % 10 == 0:
                logger.info(f"Evaluated PC {pc_idx + 1}/{n_pcs}")

        logger.info("PC evaluation completed!")

        # Calculate FDR
        from scipy.stats import false_discovery_control
        rank_fdr = false_discovery_control(rank_pvalues)
        linear_fdr = false_discovery_control(linear_pvalues)

        return {
            'rank_correlations': rank_correlations,
            'linear_correlations': linear_correlations,
            'rank_pvalues': rank_pvalues,
            'linear_pvalues': linear_pvalues,
            'rank_fdr': rank_fdr,
            'linear_fdr': linear_fdr,
            'linear_coef_pvalues_intercept': linear_coef_pvalues_intercept,
            'linear_coef_pvalues_slope': linear_coef_pvalues_slope
        }


# ============================================================================
# Main BEER Class
# ============================================================================

class BEER:
    """
    BEER: Batch EffEct Remover for single-cell data.

    This class implements the BEER algorithm for removing batch effects from
    single-cell RNA-seq and ATAC-seq data using mutual nearest neighbors and
    PC selection.

    Parameters
    ----------
    n_pcs : int, optional (default: 50)
        Number of principal components to compute.
    n_groups : int, optional (default: 30)
        Number of groups per batch for k-means clustering.
    n_variable_genes : int, optional (default: 2000)
        Number of highly variable genes to use.
    n_rounds : int, optional (default: 1)
        Number of mutual nearest neighbor detection rounds.
    correlation_method : str, optional (default: 'spearman')
        Correlation method ('spearman' or 'pearson').
    use_combat : bool, optional (default: True)
        Apply ComBat correction before PCA.
    random_seed : int, optional (default: 123)
        Random seed for reproducibility.
    rank_correlation_threshold : float, optional (default: 0.7)
        Threshold for rank correlation in PC selection.
    linear_correlation_threshold : float, optional (default: 0.7)
        Threshold for linear correlation in PC selection.

    Attributes
    ----------
    config : BEERConfig
        Configuration parameters.
    result_ : BEERResult
        Results after fitting (available after fit_transform).

    Examples
    --------
    >>> import anndata as ad
    >>> import numpy as np
    >>> from beer_python import BEER
    >>>
    >>> # Create example data
    >>> X = np.random.randn(100, 2000)
    >>> batch = np.array(['batch1'] * 50 + ['batch2'] * 50)
    >>> adata = ad.AnnData(X)
    >>> adata.obs['batch'] = batch
    >>>
    >>> # Run BEER
    >>> beer = BEER(n_pcs=50, n_groups=30)
    >>> result = beer.fit_transform(adata, batch_key='batch')
    >>>
    >>> # Access selected PCs
    >>> selected_pcs = result['selected_pcs']
    >>> corrected_adata = result['adata']
    """

    def __init__(
        self,
        n_pcs: int = 50,
        n_groups: int = 30,
        n_variable_genes: int = 2000,
        n_rounds: int = 1,
        correlation_method: str = DEFAULT_CORRELATION_METHOD,
        use_combat: bool = True,
        random_seed: int = 123,
        rank_correlation_threshold: float = 0.7,
        linear_correlation_threshold: float = 0.7,
        rank_ratio: float = 0.5,
        linear_ratio: float = 0.5
    ):
        self.config = BEERConfig(
            n_pcs=n_pcs,
            n_groups=n_groups,
            n_variable_genes=n_variable_genes,
            n_rounds=n_rounds,
            correlation_method=correlation_method,
            use_combat=use_combat,
            random_seed=random_seed,
            rank_correlation_threshold=rank_correlation_threshold,
            linear_correlation_threshold=linear_correlation_threshold,
            rank_ratio=rank_ratio,
            linear_ratio=linear_ratio
        )

        self.result_ = None
        self.mn_finder = MutualNearestNeighbors(correlation_method=correlation_method)
        self.pc_evaluator = PCEvaluator(correlation_method=correlation_method)

    def _preprocess(
        self,
        adata: ad.AnnData,
        batch_key: str,
        normalize: bool = True,
        use_raw: bool = False
    ) -> Tuple[ad.AnnData, np.ndarray]:
        """
        Preprocess data: normalization, variable gene selection.

        Parameters
        ----------
        adata : AnnData
            Input annotated data matrix.
        batch_key : str
            Key in adata.obs for batch labels.
        normalize : bool
            Whether to normalize counts.
        use_raw : bool
            If True, use adata.raw.X (assumes already normalized).

        Returns
        -------
        adata : AnnData
            Preprocessed AnnData object.
        var_gene_mask : np.ndarray
            Boolean mask for variable genes.
        """
        logger.info("=" * 80)
        logger.info(f"BEER v{BEER_VERSION} - Batch EffEct Remover")
        logger.info("=" * 80)

        # Validate input
        check_anndata_format(adata, batch_key)

        # Set random seed
        np.random.seed(self.config.random_seed)

        # Copy data to avoid modifying original
        adata = adata.copy()

        # Log info
        batches = adata.obs[batch_key].unique()
        logger.info(f"Number of batches: {len(batches)}")
        logger.info(f"Number of cells: {adata.n_obs}")
        logger.info(f"Number of genes: {adata.n_vars}")
        logger.info(f"Groups per batch (n_groups): {self.config.n_groups}")
        logger.info(f"Variable genes (n_variable_genes): {self.config.n_variable_genes}")
        logger.info(f"MN detection rounds: {self.config.n_rounds}")

        # Determine data source
        if use_raw and adata.raw is not None:
            logger.info("Using adata.raw.X (assuming pre-normalized data)")
            X_source = adata.raw.X
            normalize_data = False
        elif not normalize:
            logger.info("Using adata.X as-is (assuming already log-normalized)")
            X_source = adata.X
            normalize_data = False
        else:
            logger.info("Will normalize raw counts from adata.X")
            X_source = adata.X
            normalize_data = True

        # Step 1: Find variable genes per batch
        logger.info("\n--- Step 1: Finding variable genes ---")

        all_var_genes = np.zeros(adata.n_vars, dtype=bool)

        for batch in batches:
            logger.info(f"Processing batch: {batch}")
            batch_mask = adata.obs[batch_key] == batch
            X_batch = X_source[batch_mask, :]  # Shape: (n_cells_in_batch, n_genes)

            # Normalize if needed
            if normalize_data:
                # normalize_counts with axis=1 expects and returns (cells x genes)
                X_batch_norm = normalize_counts(X_batch, axis=1)
                # Transpose to (genes x cells) for select_variable_genes
                X_batch_norm = X_batch_norm.T
            else:
                # Data is already normalized, just convert to dense and transpose
                X_batch_norm = to_dense(X_batch).T

            # Select variable genes for this batch
            # X_batch_norm is now (genes x cells)
            var_mask = select_variable_genes(
                X_batch_norm,
                n_top_genes=self.config.n_variable_genes
            )

            # var_mask should have length adata.n_vars
            if len(var_mask) != adata.n_vars:
                raise ValueError(
                    f"Variable gene mask has wrong shape: expected {adata.n_vars}, "
                    f"got {len(var_mask)}. This suggests a data preprocessing issue."
                )

            all_var_genes |= var_mask

        logger.info(f"Total variable genes: {np.sum(all_var_genes)}")

        # Store variable genes
        adata.var['highly_variable'] = all_var_genes

        # Step 2: Store normalized data in layers
        logger.info("\n--- Step 2: Preparing normalized data layer ---")
        if normalize_data:
            # Normalize with axis=1 (cells x genes) and keep same shape
            X_norm = normalize_counts(X_source, axis=1)
            adata.layers['normalized'] = X_norm  # Store as cells x genes
        else:
            # Data is already normalized, just store it
            adata.layers['normalized'] = to_dense(X_source)
            logger.info("Stored pre-normalized data in 'normalized' layer")

        return adata, all_var_genes

    def _apply_combat(
        self,
        adata: ad.AnnData,
        batch_key: str,
        var_gene_mask: np.ndarray
    ) -> np.ndarray:
        """
        Apply ComBat batch correction.

        Parameters
        ----------
        adata : AnnData
            Input data.
        batch_key : str
            Batch key in adata.obs.
        var_gene_mask : np.ndarray
            Mask for variable genes.

        Returns
        -------
        np.ndarray
            ComBat-corrected expression matrix (cells x genes).
        """
        if not HAS_COMBAT:
            logger.warning("combat-python not installed. Skipping ComBat correction.")
            return None

        logger.info("\n--- Step 3: ComBat batch correction ---")

        # Extract variable genes
        X_var = adata.layers['normalized'][:, var_gene_mask]
        batch = adata.obs[batch_key].values

        # Apply ComBat (expects genes x samples)
        try:
            X_combat = pycombat(X_var.T, batch)
            X_combat = X_combat.T  # Convert back to cells x genes

            # Clip negative values
            X_combat = np.clip(X_combat, 0, None)

            logger.info("ComBat correction completed")
            return X_combat
        except Exception as e:
            logger.error(f"ComBat failed: {e}")
            return None

    def _compute_pca_and_umap(
        self,
        adata: ad.AnnData,
        var_gene_mask: np.ndarray,
        X_combat: Optional[np.ndarray] = None
    ) -> ad.AnnData:
        """
        Compute PCA and UMAP embeddings.

        Parameters
        ----------
        adata : AnnData
            Input data.
        var_gene_mask : np.ndarray
            Mask for variable genes.
        X_combat : np.ndarray, optional
            ComBat-corrected data for scaling.

        Returns
        -------
        AnnData
            Updated with PCA and UMAP embeddings.
        """
        logger.info("\n--- Step 4: Computing PCA ---")

        # Prepare data for PCA
        if X_combat is not None:
            # Use ComBat-corrected data
            X_pca_input = X_combat
        else:
            X_pca_input = adata.layers['normalized'][:, var_gene_mask]

        # Standardize (z-score)
        scaler = StandardScaler(with_mean=True, with_std=True)
        X_scaled = scaler.fit_transform(X_pca_input)

        # Compute PCA
        pca = PCA(
            n_components=self.config.n_pcs,
            random_state=self.config.random_seed
        )
        X_pca = pca.fit_transform(X_scaled)

        # Store in adata
        adata.obsm['X_pca'] = X_pca
        adata.uns['pca'] = {
            'variance': pca.explained_variance_,
            'variance_ratio': pca.explained_variance_ratio_
        }

        logger.info(f"PCA completed: {self.config.n_pcs} components")
        logger.info(f"Variance explained: {pca.explained_variance_ratio_[:5]}")

        # Compute UMAP if scanpy available
        if HAS_SCANPY:
            logger.info("Computing UMAP...")
            sc.pp.neighbors(adata, n_pcs=self.config.n_pcs, random_state=self.config.random_seed)
            sc.tl.umap(adata, random_state=self.config.random_seed)
            logger.info("UMAP completed")
        else:
            logger.warning("scanpy not available. Skipping UMAP computation.")
            # Create dummy UMAP using first 2 PCs
            adata.obsm['X_umap'] = X_pca[:, :2]

        return adata

    def _cluster_cells_by_batch(
        self,
        adata: ad.AnnData,
        batch_key: str
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Cluster cells within each batch using k-means on UMAP.

        Parameters
        ----------
        adata : AnnData
            Data with UMAP embeddings.
        batch_key : str
            Batch key in adata.obs.

        Returns
        -------
        group_labels : np.ndarray
            Group labels for each cell.
        batch_per_group : np.ndarray
            Batch label for each group.
        """
        logger.info("\n--- Step 5: Clustering cells within batches ---")

        umap_coords = adata.obsm['X_umap']
        batch_labels = adata.obs[batch_key].values
        unique_batches = np.unique(batch_labels)

        group_labels = np.empty(adata.n_obs, dtype=object)

        for batch in unique_batches:
            batch_mask = batch_labels == batch
            n_cells_in_batch = np.sum(batch_mask)

            # Adjust n_groups if needed
            n_groups = min(self.config.n_groups, n_cells_in_batch - 1)

            if n_groups != self.config.n_groups:
                logger.info(f"Adjusted n_groups to {n_groups} for batch {batch} (n={n_cells_in_batch} cells)")

            # K-means clustering
            umap_batch = umap_coords[batch_mask, :]
            kmeans = KMeans(
                n_clusters=n_groups,
                random_state=self.config.random_seed,
                n_init=10,
                max_iter=100
            )
            clusters = kmeans.fit_predict(umap_batch)

            # Create group labels: batch_cluster
            group_labels[batch_mask] = [f"{batch}_{c}" for c in clusters]

            logger.info(f"Created {n_groups} groups for batch {batch}")

        # Store in adata
        adata.obs['beer_group'] = group_labels

        # Get batch for each group
        unique_groups = np.unique(group_labels)
        batch_per_group = np.array([g.split('_')[0] for g in unique_groups])

        return group_labels, batch_per_group

    def _select_pcs(
        self,
        rank_correlations: np.ndarray,
        linear_correlations: np.ndarray
    ) -> np.ndarray:
        """
        Select PCs with low batch correlation.

        Parameters
        ----------
        rank_correlations : np.ndarray
            Rank correlations for each PC.
        linear_correlations : np.ndarray
            Linear correlations for each PC.

        Returns
        -------
        np.ndarray
            Indices of selected PCs.
        """
        n_pcs = len(rank_correlations)

        # Selection criteria: high rank OR above threshold
        rank_criterion = (
            (stats.rankdata(rank_correlations) >= n_pcs * self.config.rank_ratio) |
            (rank_correlations > self.config.rank_correlation_threshold)
        )

        linear_criterion = (
            (stats.rankdata(linear_correlations) >= n_pcs * self.config.linear_ratio) |
            (linear_correlations > self.config.linear_correlation_threshold)
        )

        selected = np.where(rank_criterion & linear_criterion)[0]

        logger.info(f"\nSelected {len(selected)}/{n_pcs} PCs with low batch correlation")

        return selected

    def fit_transform(
        self,
        adata: ad.AnnData,
        batch_key: str = 'batch',
        normalize: bool = True,
        remove_genes: Optional[List[str]] = None,
        use_raw: bool = False
    ) -> BEERResult:
        """
        Fit BEER model and transform data.

        Parameters
        ----------
        adata : AnnData
            Input annotated data matrix (cells x genes).
        batch_key : str
            Key in adata.obs containing batch labels.
        normalize : bool
            Whether to normalize counts (library size + log).
            If False and data is already log-normalized, set use_raw=False.
            If False and data is raw counts, will use raw counts without normalization.
        remove_genes : list of str, optional
            Gene names to remove from analysis (e.g., cell cycle genes).
        use_raw : bool
            If True, use adata.raw.X for analysis (assumes already normalized).
            If False and normalize=False, assumes adata.X is already log-normalized.

        Returns
        -------
        BEERResult
            Result object containing corrected data and statistics.
        """
        # Preprocess
        adata, var_gene_mask = self._preprocess(adata, batch_key, normalize, use_raw)

        # Remove specific genes if requested
        if remove_genes is not None:
            gene_names = adata.var_names.values
            remove_mask = np.isin(gene_names, remove_genes)
            var_gene_mask &= ~remove_mask
            logger.info(f"Removed {np.sum(remove_mask)} specified genes")
            logger.info(f"Variable genes after removal: {np.sum(var_gene_mask)}")

        # ComBat correction
        X_combat = None
        if self.config.use_combat:
            X_combat = self._apply_combat(adata, batch_key, var_gene_mask)

        # PCA and UMAP
        adata = self._compute_pca_and_umap(adata, var_gene_mask, X_combat)

        # Cluster cells by batch
        group_labels, batch_per_group = self._cluster_cells_by_batch(adata, batch_key)

        # Step 6: Find mutual nearest neighbors
        logger.info("\n--- Step 6: Finding mutual nearest neighbors ---")

        # Aggregate expression by group
        X_var = adata.layers['normalized'][:, var_gene_mask].T  # genes x cells
        X_agg, unique_groups = aggregate_by_group(X_var, group_labels, method='sum')

        # Find MN pairs
        mn_pairs = self.mn_finder.find_pairs(
            X_agg,
            unique_groups,
            batch_per_group,
            n_rounds=self.config.n_rounds
        )

        # Mark cells in MN pairs
        mn_map = np.full(adata.n_obs, 'NA', dtype=object)
        for i, group in enumerate(unique_groups):
            if i in mn_pairs[:, 0]:
                mn_map[group_labels == group] = 'V1'
            elif i in mn_pairs[:, 1]:
                mn_map[group_labels == group] = 'V2'
        adata.obs['beer_mn_map'] = mn_map

        # Step 7: Evaluate PCs
        logger.info("\n--- Step 7: Evaluating PCs for batch correlation ---")

        eval_results = self.pc_evaluator.evaluate(
            adata.obsm['X_pca'],
            group_labels,
            mn_pairs
        )

        # Select PCs
        selected_pcs = self._select_pcs(
            eval_results['rank_correlations'],
            eval_results['linear_correlations']
        )

        # Create result
        result = BEERResult(
            adata=adata,
            selected_pcs=selected_pcs,
            rank_correlations=eval_results['rank_correlations'],
            linear_correlations=eval_results['linear_correlations'],
            p_values=eval_results['rank_pvalues'],
            fdr_values=eval_results['rank_fdr'],
            mutual_pairs=mn_pairs,
            group_labels=group_labels,
            config=self.config
        )

        self.result_ = result

        logger.info("=" * 80)
        logger.info("BEER completed successfully!")
        logger.info("=" * 80)

        return result

    def refit(
        self,
        n_groups: Optional[int] = None,
        n_pcs: Optional[int] = None,
        n_rounds: Optional[int] = None,
        remove_genes: Optional[List[str]] = None
    ) -> BEERResult:
        """
        Re-run BEER with adjusted parameters (faster than full fit).

        This method reuses previous preprocessing steps.

        Parameters
        ----------
        n_groups : int, optional
            New number of groups per batch.
        n_pcs : int, optional
            New number of PCs.
        n_rounds : int, optional
            New number of MN detection rounds.
        remove_genes : list of str, optional
            Additional genes to remove.

        Returns
        -------
        BEERResult
            Updated result object.
        """
        if self.result_ is None:
            raise ValueError("Must call fit_transform before refit")

        logger.info("=" * 80)
        logger.info("ReBEER - Re-running with adjusted parameters")
        logger.info("=" * 80)

        # Update config
        if n_groups is not None:
            self.config.n_groups = n_groups
        if n_pcs is not None:
            self.config.n_pcs = n_pcs
        if n_rounds is not None:
            self.config.n_rounds = n_rounds

        # Get previous data
        adata = self.result_.adata.copy()
        batch_key = 'batch'  # Assuming batch_key was 'batch'

        # Re-compute PCA if n_pcs changed
        if n_pcs is not None and n_pcs != adata.obsm['X_pca'].shape[1]:
            logger.info(f"Re-computing PCA with {n_pcs} components...")
            var_gene_mask = adata.var['highly_variable'].values
            adata = self._compute_pca_and_umap(adata, var_gene_mask)

        # Re-cluster and evaluate
        group_labels, batch_per_group = self._cluster_cells_by_batch(adata, batch_key)

        # Aggregate and find MN pairs
        X_var = adata.layers['normalized'][:, adata.var['highly_variable']].T
        X_agg, unique_groups = aggregate_by_group(X_var, group_labels, method='sum')

        mn_pairs = self.mn_finder.find_pairs(
            X_agg,
            unique_groups,
            batch_per_group,
            n_rounds=self.config.n_rounds
        )

        # Evaluate PCs
        eval_results = self.pc_evaluator.evaluate(
            adata.obsm['X_pca'],
            group_labels,
            mn_pairs
        )

        selected_pcs = self._select_pcs(
            eval_results['rank_correlations'],
            eval_results['linear_correlations']
        )

        result = BEERResult(
            adata=adata,
            selected_pcs=selected_pcs,
            rank_correlations=eval_results['rank_correlations'],
            linear_correlations=eval_results['linear_correlations'],
            p_values=eval_results['rank_pvalues'],
            fdr_values=eval_results['rank_fdr'],
            mutual_pairs=mn_pairs,
            group_labels=group_labels,
            config=self.config
        )

        self.result_ = result

        logger.info("ReBEER completed!")

        return result


# ============================================================================
# Enhancement Functions
# ============================================================================

def apply_combat_to_pca(
    adata: ad.AnnData,
    batch_key: str = 'batch',
    pca_key: str = 'X_pca'
) -> ad.AnnData:
    """
    Apply ComBat correction to PCA embeddings.

    Parameters
    ----------
    adata : AnnData
        Data with PCA embeddings.
    batch_key : str
        Batch key in adata.obs.
    pca_key : str
        Key for PCA embeddings in adata.obsm.

    Returns
    -------
    AnnData
        Data with ComBat-corrected PCA.
    """
    if not HAS_COMBAT:
        raise ImportError("combat-python required. Install: pip install combat")

    logger.info("Applying ComBat correction to PCA embeddings...")

    X_pca = adata.obsm[pca_key]
    batch = adata.obs[batch_key].values

    # Apply ComBat (expects genes x samples, so transpose)
    X_pca_combat = pycombat(X_pca.T, batch).T

    # Update adata
    adata.obsm[pca_key] = X_pca_combat

    logger.info("ComBat correction completed")

    return adata


def apply_bbknn(
    adata: ad.AnnData,
    selected_pcs: np.ndarray,
    batch_key: str = 'batch',
    neighbors_within_batch: int = 3,
    n_trees: int = 10,
    n_umap_components: int = 2
) -> np.ndarray:
    """
    Apply BBKNN enhancement for better batch mixing.

    Parameters
    ----------
    adata : AnnData
        Data with PCA embeddings.
    selected_pcs : np.ndarray
        Indices of selected PCs to use.
    batch_key : str
        Batch key in adata.obs.
    neighbors_within_batch : int
        Number of neighbors within each batch.
    n_trees : int
        Number of trees for Annoy index.
    n_umap_components : int
        Number of UMAP dimensions.

    Returns
    -------
    np.ndarray
        BBKNN-enhanced UMAP coordinates.
    """
    if not HAS_BBKNN or not HAS_SCANPY:
        raise ImportError("bbknn and scanpy required. Install: pip install bbknn scanpy")

    logger.info("Applying BBKNN enhancement...")

    # Create temporary copy
    adata_temp = adata.copy()

    # Use only selected PCs
    adata_temp.obsm['X_pca'] = adata.obsm['X_pca'][:, selected_pcs]

    # Run BBKNN
    bbknn.bbknn(
        adata_temp,
        batch_key=batch_key,
        neighbors_within_batch=neighbors_within_batch,
        n_pcs=len(selected_pcs),
        annoy_n_trees=n_trees
    )

    # Compute UMAP
    sc.tl.umap(adata_temp, n_components=n_umap_components)

    umap_coords = adata_temp.obsm['X_umap']

    logger.info("BBKNN enhancement completed")

    return umap_coords


# ============================================================================
# Convenience Functions
# ============================================================================

def read_data(
    file_path: str,
    batch_column: str = 'batch',
    sep: str = '\t'
) -> ad.AnnData:
    """
    Read expression data from file and create AnnData object.

    Parameters
    ----------
    file_path : str
        Path to expression matrix file.
    batch_column : str
        Name of batch column in metadata.
    sep : str
        Column separator.

    Returns
    -------
    AnnData
        Annotated data object.
    """
    # Read data
    df = pd.read_csv(file_path, sep=sep, index_col=0)

    # Transpose if needed (expect cells x genes)
    if df.shape[0] > df.shape[1]:
        df = df.T

    # Create AnnData
    adata = ad.AnnData(df.values)
    adata.obs_names = df.index
    adata.var_names = df.columns

    return adata


def plot_correlation_scatter(
    result: BEERResult,
    figsize: Tuple[int, int] = (8, 6),
    save: Optional[str] = None
):
    """
    Plot rank vs linear correlation for PC selection visualization.

    Parameters
    ----------
    result : BEERResult
        BEER result object.
    figsize : tuple
        Figure size.
    save : str, optional
        Path to save figure.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.error("matplotlib required for plotting")
        return

    rank_cor = result.rank_correlations
    linear_cor = result.linear_correlations
    selected = result.selected_pcs

    fig, ax = plt.subplots(figsize=figsize)

    # Plot all PCs
    ax.scatter(rank_cor, linear_cor, c='gray', alpha=0.5, label='Not selected')

    # Highlight selected PCs
    ax.scatter(
        rank_cor[selected],
        linear_cor[selected],
        c='red',
        alpha=0.7,
        label='Selected'
    )

    ax.axhline(
        y=result.config.linear_correlation_threshold,
        color='blue',
        linestyle='--',
        label='Linear threshold'
    )
    ax.axvline(
        x=result.config.rank_correlation_threshold,
        color='green',
        linestyle='--',
        label='Rank threshold'
    )

    ax.set_xlabel('Rank Correlation')
    ax.set_ylabel('Linear Correlation')
    ax.set_title('PC Selection: Batch Correlation')
    ax.legend()
    ax.grid(True, alpha=0.3)

    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        logger.info(f"Figure saved to {save}")

    plt.show()


# ============================================================================
# Main Entry Point (for CLI usage)
# ============================================================================

def main():
    """Command-line interface for BEER."""
    import argparse

    parser = argparse.ArgumentParser(
        description='BEER: Batch EffEct Remover for Single-Cell Data',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('input', help='Input h5ad file')
    parser.add_argument('--output', '-o', default='beer_output.h5ad', help='Output h5ad file')
    parser.add_argument('--batch-key', default='batch', help='Batch column in adata.obs')
    parser.add_argument('--n-pcs', type=int, default=50, help='Number of PCs')
    parser.add_argument('--n-groups', type=int, default=30, help='Groups per batch')
    parser.add_argument('--n-variable-genes', type=int, default=2000, help='Variable genes')
    parser.add_argument('--n-rounds', type=int, default=1, help='MN detection rounds')
    parser.add_argument('--no-combat', action='store_true', help='Disable ComBat correction')
    parser.add_argument('--seed', type=int, default=123, help='Random seed')
    parser.add_argument('--plot', action='store_true', help='Generate correlation plot')

    args = parser.parse_args()

    # Load data
    logger.info(f"Loading data from {args.input}")
    adata = ad.read_h5ad(args.input)

    # Run BEER
    beer = BEER(
        n_pcs=args.n_pcs,
        n_groups=args.n_groups,
        n_variable_genes=args.n_variable_genes,
        n_rounds=args.n_rounds,
        use_combat=not args.no_combat,
        random_seed=args.seed
    )

    result = beer.fit_transform(adata, batch_key=args.batch_key)

    # Save results
    logger.info(f"Saving results to {args.output}")
    result.adata.write_h5ad(args.output)

    # Plot if requested
    if args.plot:
        plot_file = args.output.replace('.h5ad', '_correlation_plot.png')
        plot_correlation_scatter(result, save=plot_file)

    logger.info(f"Selected PCs: {result.selected_pcs}")
    logger.info("BEER analysis completed successfully!")


if __name__ == '__main__':
    main()
