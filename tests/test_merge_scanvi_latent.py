import importlib.util
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "merge_annotations.py"
spec = importlib.util.spec_from_file_location("merge_annotations", SCRIPT)
merge_annotations = importlib.util.module_from_spec(spec)
spec.loader.exec_module(merge_annotations)


class MergeScanviLatentTests(unittest.TestCase):
    def test_copies_latent_obsm_and_metadata(self):
        query = ad.AnnData(obs=pd.DataFrame(index=["cell-1", "cell-2"]))
        latent = ad.AnnData(obs=pd.DataFrame(index=["cell-1", "cell-2"]))
        latent.obsm["X_scanvi"] = np.array([[1.0, 2.0], [3.0, 4.0]])
        latent.obsm["X_scvi"] = np.array([[5.0, 6.0], [7.0, 8.0]])
        latent.uns["scanvi_latent"] = {"obsm_keys": ["X_scanvi", "X_scvi"]}

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "scanvi_latent.h5ad"
            latent.write_h5ad(path)
            merge_annotations.merge_scanvi_latent_obsm(query, str(path))

        np.testing.assert_array_equal(query.obsm["X_scanvi"], latent.obsm["X_scanvi"])
        np.testing.assert_array_equal(query.obsm["X_scvi"], latent.obsm["X_scvi"])
        self.assertEqual(list(query.uns["scanvi_latent"]["obsm_keys"]), ["X_scanvi", "X_scvi"])

    def test_rejects_latent_with_different_cell_order(self):
        query = ad.AnnData(obs=pd.DataFrame(index=["cell-1", "cell-2"]))
        latent = ad.AnnData(obs=pd.DataFrame(index=["cell-2", "cell-1"]))
        latent.obsm["X_scanvi"] = np.array([[1.0, 2.0], [3.0, 4.0]])

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "scanvi_latent.h5ad"
            latent.write_h5ad(path)
            with self.assertRaisesRegex(ValueError, "cell order"):
                merge_annotations.merge_scanvi_latent_obsm(query, str(path))


if __name__ == "__main__":
    unittest.main()
