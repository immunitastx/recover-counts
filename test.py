import pandas as pd
import scanpy as sc
from anndata import AnnData
from recover_counts_from_log_normalized_data import recover_counts
import numpy as np


from unittest import TestCase

class TestRecoverCounts(TestCase):
    def test_recover_countss(self):
        # Load test data
        df_counts = pd.read_csv('./GSM3516673_MSK_LX682_NORMAL_dense.csv.gz', sep=',', index_col=0)
        ad = AnnData(df_counts)

        # Compute true size-factors
        true_size_factors = np.sum(ad.X, axis=1)

        # Normalize
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)
        X_norm = ad.X

        # Recover the counts
        X_counts_recovered, size_factors = recover_counts(X_norm, 1e4, 1e5, log_base=None, verbose=True)

        # Check output matches input
        self.assertTrue(
            tuple(true_size_factors) == tuple(size_factors)
        ) 
