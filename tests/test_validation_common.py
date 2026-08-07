import importlib.util
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import h5py
import numpy as np


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


class ObsColumnTests(unittest.TestCase):
    def test_missing_obs_column_is_reported_with_its_label(self):
        adata = ad.AnnData(obs={"cell_type": ["A", "B"]})

        with self.assertRaisesRegex(ValueError, "Missing required obs column 'donor' in reference h5ad"):
            validation_common.require_obs_column(adata, "donor", "reference h5ad")


if __name__ == "__main__":
    unittest.main()
