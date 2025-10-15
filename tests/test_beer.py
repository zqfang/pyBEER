"""
Unit Tests for BEER Python Implementation
==========================================

Test suite for the BEER batch effect removal algorithm.

Author: Feng Zhang (original R), Python tests
Date: October 2025
"""

import unittest
import warnings
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad

# Import BEER components
from beer import (
    BEER,
    BEERConfig,
    BEERResult,
    check_anndata_format,
    to_dense,
    normalize_counts,
    select_variable_genes,
    aggregate_by_group,
    MutualNearestNeighbors,
    PCEvaluator,
)

# Suppress warnings during tests
warnings.filterwarnings('ignore')


class TestUtilityFunctions(unittest.TestCase):
    """Test utility functions."""

    def test_to_dense(self):
        """Test sparse to dense conversion."""
        # Create sparse matrix
        X_sparse = sparse.random(100, 50, density=0.3, format='csr')
        X_dense = to_dense(X_sparse)

        self.assertIsInstance(X_dense, np.ndarray)
        self.assertEqual(X_dense.shape, X_sparse.shape)
        np.testing.assert_array_almost_equal(X_dense, X_sparse.toarray())

    def test_to_dense_with_dense(self):
        """Test that dense arrays pass through."""
        X = np.random.randn(100, 50)
        X_out = to_dense(X)

        np.testing.assert_array_equal(X, X_out)

    def test_normalize_counts(self):
        """Test count normalization."""
        X = np.random.poisson(5, size=(100, 50))

        X_norm = normalize_counts(X, target_sum=1e4, log_transform=True)

        self.assertEqual(X_norm.shape, X.shape)
        self.assertTrue(np.all(X_norm >= 0))
        self.assertTrue(np.all(np.isfinite(X_norm)))

    def test_normalize_counts_no_log(self):
        """Test normalization without log transform."""
        X = np.random.poisson(5, size=(100, 50))

        X_norm = normalize_counts(X, target_sum=1e4, log_transform=False)

        # Check that cell sums are approximately target_sum
        cell_sums = np.sum(X_norm, axis=0)
        np.testing.assert_array_almost_equal(cell_sums, 1e4 * np.ones(50), decimal=5)

    def test_select_variable_genes(self):
        """Test variable gene selection."""
        X = np.random.randn(1000, 100)

        var_mask = select_variable_genes(X, n_top_genes=200)

        self.assertEqual(len(var_mask), 1000)
        self.assertEqual(np.sum(var_mask), 200)
        self.assertTrue(var_mask.dtype == bool)

    def test_aggregate_by_group_sum(self):
        """Test aggregation by sum."""
        X = np.random.randn(100, 50)
        groups = np.array(['A'] * 25 + ['B'] * 25)

        X_agg, group_names = aggregate_by_group(X, groups, method='sum')

        self.assertEqual(X_agg.shape, (100, 2))
        self.assertEqual(len(group_names), 2)
        np.testing.assert_array_equal(group_names, ['A', 'B'])

        # Check sum correctness
        np.testing.assert_array_almost_equal(
            X_agg[:, 0],
            np.sum(X[:, :25], axis=1)
        )

    def test_aggregate_by_group_mean(self):
        """Test aggregation by mean."""
        X = np.random.randn(100, 50)
        groups = np.array(['A'] * 25 + ['B'] * 25)

        X_agg, group_names = aggregate_by_group(X, groups, method='mean')

        self.assertEqual(X_agg.shape, (100, 2))

        # Check mean correctness
        np.testing.assert_array_almost_equal(
            X_agg[:, 0],
            np.mean(X[:, :25], axis=1)
        )


class TestAnnDataValidation(unittest.TestCase):
    """Test AnnData validation functions."""

    def test_check_anndata_format_valid(self):
        """Test validation with valid AnnData."""
        X = np.random.randn(100, 50)
        adata = ad.AnnData(X)
        adata.obs['batch'] = np.array(['A'] * 50 + ['B'] * 50)

        # Should not raise
        check_anndata_format(adata, 'batch')

    def test_check_anndata_format_missing_batch(self):
        """Test validation with missing batch key."""
        X = np.random.randn(100, 50)
        adata = ad.AnnData(X)

        with self.assertRaises(ValueError):
            check_anndata_format(adata, 'batch')

    def test_check_anndata_format_one_batch(self):
        """Test validation with only one batch."""
        X = np.random.randn(100, 50)
        adata = ad.AnnData(X)
        adata.obs['batch'] = np.array(['A'] * 100)

        with self.assertRaises(ValueError):
            check_anndata_format(adata, 'batch')


class TestMutualNearestNeighbors(unittest.TestCase):
    """Test mutual nearest neighbor detection."""

    def setUp(self):
        """Set up test fixtures."""
        self.mn_finder = MutualNearestNeighbors(correlation_method='spearman')

    def test_compute_correlation_spearman(self):
        """Test Spearman correlation computation."""
        X = np.random.randn(100, 20)
        corr = self.mn_finder._compute_correlation(X)

        self.assertEqual(corr.shape, (20, 20))
        np.testing.assert_array_almost_equal(np.diag(corr), np.ones(20))

    def test_compute_correlation_pearson(self):
        """Test Pearson correlation computation."""
        mn_finder_pearson = MutualNearestNeighbors(correlation_method='pearson')
        X = np.random.randn(100, 20)
        corr = mn_finder_pearson._compute_correlation(X)

        self.assertEqual(corr.shape, (20, 20))
        np.testing.assert_array_almost_equal(np.diag(corr), np.ones(20))

    def test_extract_batch_from_group(self):
        """Test batch extraction from group name."""
        group_name = "Batch1_5"
        batch = self.mn_finder._extract_batch_from_group(group_name)

        self.assertEqual(batch, "Batch1")

    def test_find_pairs(self):
        """Test MN pair detection."""
        np.random.seed(123)

        # Create simple test case with 2 batches, 4 groups each
        X = np.random.randn(100, 8)
        group_labels = np.array([f'G{i}' for i in range(8)])
        batch_labels = np.array(['A', 'A', 'A', 'A', 'B', 'B', 'B', 'B'])

        pairs = self.mn_finder.find_pairs(
            X,
            group_labels,
            batch_labels,
            n_rounds=1
        )

        self.assertTrue(isinstance(pairs, np.ndarray))
        self.assertTrue(len(pairs) > 0)
        self.assertEqual(pairs.shape[1], 2)


class TestPCEvaluator(unittest.TestCase):
    """Test PC evaluation."""

    def setUp(self):
        """Set up test fixtures."""
        self.evaluator = PCEvaluator(correlation_method='spearman')

    def test_evaluate(self):
        """Test PC evaluation."""
        np.random.seed(123)

        # Create test data
        n_cells = 100
        n_pcs = 10
        pca_embeddings = np.random.randn(n_cells, n_pcs)
        group_labels = np.array(['A_1'] * 25 + ['A_2'] * 25 + ['B_1'] * 25 + ['B_2'] * 25)
        mn_pairs = np.array([[0, 2], [1, 3]])  # A_1 <-> B_1, A_2 <-> B_2

        results = self.evaluator.evaluate(pca_embeddings, group_labels, mn_pairs)

        self.assertIn('rank_correlations', results)
        self.assertIn('linear_correlations', results)
        self.assertIn('rank_pvalues', results)
        self.assertIn('linear_pvalues', results)

        self.assertEqual(len(results['rank_correlations']), n_pcs)
        self.assertTrue(np.all(np.abs(results['rank_correlations']) <= 1))
        self.assertTrue(np.all(np.abs(results['linear_correlations']) <= 1))


class TestBEERConfig(unittest.TestCase):
    """Test BEER configuration."""

    def test_default_config(self):
        """Test default configuration."""
        config = BEERConfig()

        self.assertEqual(config.n_pcs, 50)
        self.assertEqual(config.n_groups, 30)
        self.assertEqual(config.n_variable_genes, 2000)
        self.assertEqual(config.n_rounds, 1)
        self.assertTrue(config.use_combat)

    def test_custom_config(self):
        """Test custom configuration."""
        config = BEERConfig(
            n_pcs=100,
            n_groups=50,
            use_combat=False
        )

        self.assertEqual(config.n_pcs, 100)
        self.assertEqual(config.n_groups, 50)
        self.assertFalse(config.use_combat)


class TestBEER(unittest.TestCase):
    """Test main BEER class."""

    def setUp(self):
        """Set up test fixtures."""
        np.random.seed(123)

        # Create simple test data
        n_cells_per_batch = 100
        n_genes = 500

        X1 = np.random.negative_binomial(5, 0.3, size=(n_cells_per_batch, n_genes))
        X2 = np.random.negative_binomial(5, 0.3, size=(n_cells_per_batch, n_genes))

        X = np.vstack([X1, X2]).astype(float)
        batch = np.array(['Batch1'] * n_cells_per_batch + ['Batch2'] * n_cells_per_batch)

        self.adata = ad.AnnData(X)
        self.adata.obs['batch'] = batch
        self.adata.var_names = [f'Gene_{i}' for i in range(n_genes)]
        self.adata.obs_names = [f'Cell_{i}' for i in range(len(batch))]

    def test_beer_initialization(self):
        """Test BEER initialization."""
        beer = BEER(n_pcs=30, n_groups=20)

        self.assertEqual(beer.config.n_pcs, 30)
        self.assertEqual(beer.config.n_groups, 20)
        self.assertIsNone(beer.result_)

    def test_beer_fit_transform_basic(self):
        """Test basic BEER fit_transform."""
        beer = BEER(
            n_pcs=20,
            n_groups=10,
            n_variable_genes=300,
            use_combat=False  # Skip ComBat for speed
        )

        result = beer.fit_transform(self.adata, batch_key='batch', normalize=True)

        self.assertIsInstance(result, BEERResult)
        self.assertEqual(result.adata.n_obs, self.adata.n_obs)
        self.assertTrue('X_pca' in result.adata.obsm)
        self.assertTrue('X_umap' in result.adata.obsm)
        self.assertTrue(len(result.selected_pcs) > 0)
        self.assertTrue(len(result.selected_pcs) <= 20)

    def test_beer_stores_metadata(self):
        """Test that BEER stores metadata correctly."""
        beer = BEER(n_pcs=20, n_groups=10, use_combat=False)
        result = beer.fit_transform(self.adata, batch_key='batch')

        # Check metadata
        self.assertTrue('beer_group' in result.adata.obs)
        self.assertTrue('beer_mn_map' in result.adata.obs)
        self.assertTrue('highly_variable' in result.adata.var)

    def test_beer_selected_pcs_range(self):
        """Test that selected PCs are in valid range."""
        beer = BEER(n_pcs=20, use_combat=False)
        result = beer.fit_transform(self.adata, batch_key='batch')

        self.assertTrue(np.all(result.selected_pcs >= 0))
        self.assertTrue(np.all(result.selected_pcs < 20))

    def test_beer_refit(self):
        """Test BEER refitting."""
        beer = BEER(n_pcs=20, n_groups=10, use_combat=False)
        result1 = beer.fit_transform(self.adata, batch_key='batch')

        # Refit with different parameters
        result2 = beer.refit(n_groups=15, n_rounds=2)

        self.assertIsInstance(result2, BEERResult)
        self.assertEqual(result2.config.n_groups, 15)
        self.assertEqual(result2.config.n_rounds, 2)

    def test_beer_with_remove_genes(self):
        """Test BEER with gene removal."""
        remove_genes = ['Gene_0', 'Gene_1', 'Gene_2']

        beer = BEER(n_pcs=20, use_combat=False)
        result = beer.fit_transform(
            self.adata,
            batch_key='batch',
            remove_genes=remove_genes
        )

        self.assertIsInstance(result, BEERResult)

    def test_beer_correlation_values(self):
        """Test that correlation values are in valid range."""
        beer = BEER(n_pcs=20, use_combat=False)
        result = beer.fit_transform(self.adata, batch_key='batch')

        self.assertTrue(np.all(np.abs(result.rank_correlations) <= 1))
        self.assertTrue(np.all(np.abs(result.linear_correlations) <= 1))
        self.assertTrue(np.all(result.p_values >= 0))
        self.assertTrue(np.all(result.p_values <= 1))


class TestBEEREdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""

    def test_small_batch_size(self):
        """Test with very small batch size."""
        np.random.seed(123)

        # Create data with small batches
        X = np.random.randn(30, 100)
        batch = np.array(['A'] * 15 + ['B'] * 15)

        adata = ad.AnnData(X)
        adata.obs['batch'] = batch

        # Should adjust n_groups automatically
        beer = BEER(n_pcs=10, n_groups=20, use_combat=False)  # n_groups > batch size
        result = beer.fit_transform(adata, batch_key='batch')

        self.assertIsInstance(result, BEERResult)

    def test_different_batch_sizes(self):
        """Test with unbalanced batch sizes."""
        np.random.seed(123)

        X1 = np.random.randn(100, 200)
        X2 = np.random.randn(50, 200)
        X = np.vstack([X1, X2])

        batch = np.array(['A'] * 100 + ['B'] * 50)

        adata = ad.AnnData(X)
        adata.obs['batch'] = batch

        beer = BEER(n_pcs=10, n_groups=20, use_combat=False)
        result = beer.fit_transform(adata, batch_key='batch')

        self.assertIsInstance(result, BEERResult)

    def test_sparse_input(self):
        """Test with sparse input matrix."""
        np.random.seed(123)

        # Create sparse matrix
        X_sparse = sparse.random(200, 300, density=0.1, format='csr')
        batch = np.array(['A'] * 100 + ['B'] * 100)

        adata = ad.AnnData(X_sparse)
        adata.obs['batch'] = batch

        beer = BEER(n_pcs=10, n_groups=10, use_combat=False)
        result = beer.fit_transform(adata, batch_key='batch')

        self.assertIsInstance(result, BEERResult)


class TestBEERReproducibility(unittest.TestCase):
    """Test reproducibility with random seeds."""

    def test_reproducibility_with_seed(self):
        """Test that results are reproducible with same seed."""
        np.random.seed(123)

        X = np.random.randn(100, 200)
        batch = np.array(['A'] * 50 + ['B'] * 50)

        adata = ad.AnnData(X)
        adata.obs['batch'] = batch

        # Run twice with same seed
        beer1 = BEER(n_pcs=20, n_groups=10, random_seed=123, use_combat=False)
        result1 = beer1.fit_transform(adata, batch_key='batch')

        beer2 = BEER(n_pcs=20, n_groups=10, random_seed=123, use_combat=False)
        result2 = beer2.fit_transform(adata, batch_key='batch')

        # Check that selected PCs are the same
        np.testing.assert_array_equal(result1.selected_pcs, result2.selected_pcs)

        # Check that correlations are the same
        np.testing.assert_array_almost_equal(
            result1.rank_correlations,
            result2.rank_correlations
        )


# ===========================================================================
# Test Suite Runner
# ===========================================================================

def run_tests():
    """Run all tests."""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestUtilityFunctions))
    suite.addTests(loader.loadTestsFromTestCase(TestAnnDataValidation))
    suite.addTests(loader.loadTestsFromTestCase(TestMutualNearestNeighbors))
    suite.addTests(loader.loadTestsFromTestCase(TestPCEvaluator))
    suite.addTests(loader.loadTestsFromTestCase(TestBEERConfig))
    suite.addTests(loader.loadTestsFromTestCase(TestBEER))
    suite.addTests(loader.loadTestsFromTestCase(TestBEEREdgeCases))
    suite.addTests(loader.loadTestsFromTestCase(TestBEERReproducibility))

    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    return result


if __name__ == '__main__':
    print("=" * 80)
    print("BEER Python Implementation - Test Suite")
    print("=" * 80)
    print()

    result = run_tests()

    print()
    print("=" * 80)
    print("Test Summary")
    print("=" * 80)
    print(f"Tests run: {result.testsRun}")
    print(f"Successes: {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print("=" * 80)

    # Exit with appropriate code
    exit(0 if result.wasSuccessful() else 1)
