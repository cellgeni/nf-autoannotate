import importlib.util
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
spec = importlib.util.spec_from_file_location("validation_common", SCRIPTS / "validation_common.py")
validation_common = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validation_common)


class RawCountMatrixTests(unittest.TestCase):
    def test_raw_count_matrix_rejects_fractional_values(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "query.h5ad"
            with h5py.File(path, "w") as handle:
                handle.create_dataset("X", data=np.array([[0.0, 1.5], [2.0, 3.0]]))

            with h5py.File(path, "r") as handle:
                with self.assertRaisesRegex(ValueError, "does not look like raw counts"):
                    validation_common.require_raw_count_matrix(handle, "query h5ad")

    def test_raw_count_matrix_accepts_integer_float_values(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "query.h5ad"
            with h5py.File(path, "w") as handle:
                handle.create_dataset("X", data=np.array([[0.0, 1.0], [2.0, 3.0]]))

            with h5py.File(path, "r") as handle:
                validation_common.require_raw_count_matrix(handle, "query h5ad")


class CellTypistMatrixTests(unittest.TestCase):
    COUNTS = np.array([[5, 0, 15], [1, 9, 0], [0, 4, 4]], dtype=np.float32)

    @staticmethod
    def normalise(counts):
        return np.log1p(counts / counts.sum(axis=1, keepdims=True) * 1e4).astype(np.float32)

    def select(self, x, layers=None, raw=None, sparse_x=False):
        import scipy.sparse as sp

        adata = ad.AnnData(X=sp.csr_matrix(x) if sparse_x else x)
        if layers is not None:
            adata.layers["counts"] = layers
        if raw is not None:
            adata.raw = ad.AnnData(X=raw)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "query.h5ad"
            adata.write_h5ad(path)
            with h5py.File(path, "r") as handle:
                return validation_common.select_celltypist_matrix(handle, "query h5ad")

    def test_normalised_x_is_used_as_it_is(self):
        for sparse_x in (False, True):
            self.assertEqual(self.select(self.normalise(self.COUNTS), sparse_x=sparse_x), ("X", True))

    def test_normalised_x_wins_over_raw_counts_elsewhere(self):
        self.assertEqual(
            self.select(self.normalise(self.COUNTS), layers=self.COUNTS, raw=self.COUNTS),
            ("X", True),
        )

    def test_csc_matrices_are_sampled_without_loading_them_whole(self):
        import scipy.sparse as sp

        rng = np.random.default_rng(0)
        counts = rng.poisson(1.0, size=(700, 30)).astype(np.float32)
        counts[counts.sum(axis=1) == 0, 0] = 1
        for matrix, expected in ((counts, ("X", False)), (self.normalise(counts), ("X", True))):
            adata = ad.AnnData(X=sp.csc_matrix(matrix))
            with tempfile.TemporaryDirectory() as tmpdir:
                path = Path(tmpdir) / "query.h5ad"
                adata.write_h5ad(path)
                with h5py.File(path, "r") as handle:
                    self.assertEqual(validation_common.select_celltypist_matrix(handle, "query"), expected)

    def test_exact_totals_check_accepts_normalised_and_rejects_scaled_matrices(self):
        import scipy.sparse as sp

        normalised = self.normalise(self.COUNTS)
        scaled = normalised * 0.5
        for convert in (np.asarray, sp.csr_matrix, sp.csc_matrix):
            self.assertTrue(validation_common.log1p_totals_are_normalised(convert(normalised)))
            self.assertFalse(validation_common.log1p_totals_are_normalised(convert(scaled)))

    def test_raw_counts_in_x_are_normalised(self):
        self.assertEqual(self.select(self.COUNTS), ("X", False))

    def test_counts_layer_is_used_when_x_is_neither(self):
        scaled = self.COUNTS / 3 + 0.5
        self.assertEqual(self.select(scaled, layers=self.COUNTS), ("layers/counts", False))

    def test_raw_x_is_used_when_there_is_no_counts_layer(self):
        scaled = self.COUNTS / 3 + 0.5
        self.assertEqual(self.select(scaled, raw=self.COUNTS), ("raw/X", False))

    def test_missing_raw_counts_is_an_error(self):
        with self.assertRaisesRegex(ValueError, "layers\\['counts'\\] or adata.raw.X"):
            self.select(self.COUNTS / 3 + 0.5)


class ToolSelectionTests(unittest.TestCase):
    def test_parse_tools_normalizes_case_and_whitespace(self):
        self.assertEqual(
            validation_common.parse_tools(" CellTypist , panhumanpy "),
            {"celltypist", "panhumanpy"},
        )

    def test_parse_tools_rejects_unknown_names(self):
        with self.assertRaisesRegex(ValueError, "Unknown --tools entries: singler"):
            validation_common.parse_tools("celltypist,singler")

    def test_parse_tools_rejects_empty_selection(self):
        with self.assertRaisesRegex(ValueError, "must name at least one"):
            validation_common.parse_tools("  ,  ")

    def test_panhumanpy_only_selection_needs_no_reference(self):
        tools = validation_common.parse_tools("panhumanpy")

        self.assertFalse(validation_common.REFERENCE_TOOLS & tools)

    def test_celltypist_selection_needs_a_reference(self):
        tools = validation_common.parse_tools("celltypist,panhumanpy")

        self.assertTrue(validation_common.REFERENCE_TOOLS & tools)


class GeneIdentifierTests(unittest.TestCase):
    def test_ensembl_identifiers_are_detected(self):
        self.assertEqual(
            validation_common.infer_gene_id_type(["ENSG00000141510", "ENSG00000012048"]),
            "ensembl",
        )

    def test_symbols_are_detected(self):
        self.assertEqual(validation_common.infer_gene_id_type(["TP53", "BRCA1"]), "symbol")


class H5adMetadataTests(unittest.TestCase):
    def write_h5ad(self, path: Path) -> None:
        adata = ad.AnnData(
            X=np.zeros((3, 2), dtype=np.float32),
            obs=pd.DataFrame(
                {"cell_type": pd.Categorical(["A", "B", "A"]), "donor": ["d1", "d2", "d3"]},
                index=["cell-1", "cell-2", "cell-3"],
            ),
            var=pd.DataFrame({"gene_symbol": ["A", "B"]}, index=["GENE1", "GENE2"]),
            layers={"counts": np.zeros((3, 2), dtype=np.float32)},
        )
        adata.write_h5ad(path)

    def test_reads_indices_and_column_and_layer_names(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "query.h5ad"
            self.write_h5ad(path)

            with h5py.File(path, "r") as handle:
                metadata = validation_common.read_h5ad_metadata(handle, "query h5ad")
                cell_types = validation_common.read_obs_column(handle, "cell_type")

        self.assertEqual(list(metadata.obs_names), ["cell-1", "cell-2", "cell-3"])
        self.assertEqual(list(metadata.var_names), ["GENE1", "GENE2"])
        self.assertEqual(metadata.obs_columns, {"cell_type", "donor"})
        self.assertEqual(metadata.var_columns, {"gene_symbol"})
        self.assertEqual(metadata.layers, {"counts"})
        self.assertEqual(metadata.n_obs, 3)
        self.assertEqual(list(cell_types.astype(str)), ["A", "B", "A"])

    def test_missing_obs_table_is_reported(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "query.h5ad"
            with h5py.File(path, "w") as handle:
                handle.create_dataset("X", data=np.zeros((2, 2)))

            with h5py.File(path, "r") as handle:
                with self.assertRaisesRegex(ValueError, "Missing required /obs table in query h5ad"):
                    validation_common.read_h5ad_metadata(handle, "query h5ad")


class ObsColumnTests(unittest.TestCase):
    def test_missing_obs_column_is_reported_with_its_label(self):
        metadata = validation_common.H5adMetadata(
            obs_names=pd.Index(["cell-1", "cell-2"]),
            var_names=pd.Index(["GENE1"]),
            obs_columns={"cell_type"},
            var_columns=set(),
            layers=set(),
        )

        with self.assertRaisesRegex(ValueError, "Missing required obs column 'donor' in reference h5ad"):
            validation_common.require_obs_column(metadata, "donor", "reference h5ad")


if __name__ == "__main__":
    unittest.main()
