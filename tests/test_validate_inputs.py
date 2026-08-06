import importlib.util
import json
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path

import anndata as ad
import h5py
import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "validate_inputs.py"
spec = importlib.util.spec_from_file_location("validate_inputs", SCRIPT)
validate_inputs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validate_inputs)


class ValidateInputsTests(unittest.TestCase):
    def test_query_cluster_column_is_required_when_missing_cluster_computation_is_disabled(self):
        query = ad.AnnData(obs={"cell_type": ["A", "B"]})

        with self.assertRaisesRegex(ValueError, "Missing required obs column 'leiden' in query h5ad"):
            validate_inputs.require_query_cluster_column(
                query, "leiden", compute_missing_clusters=False
            )

    def test_query_cluster_column_can_be_missing_when_cluster_computation_is_enabled(self):
        query = ad.AnnData(obs={"cell_type": ["A", "B"]})

        validate_inputs.require_query_cluster_column(query, "leiden", compute_missing_clusters=True)

    def test_raw_count_matrix_rejects_fractional_values(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "query.h5ad"
            with h5py.File(path, "w") as handle:
                handle.create_dataset("X", data=np.array([[0.0, 1.5], [2.0, 3.0]]))

            with h5py.File(path, "r") as handle:
                with self.assertRaisesRegex(ValueError, "does not look like raw counts"):
                    validate_inputs.require_raw_count_matrix(handle, "query h5ad")

    def test_raw_count_matrix_accepts_integer_float_values(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "query.h5ad"
            with h5py.File(path, "w") as handle:
                handle.create_dataset("X", data=np.array([[0.0, 1.0], [2.0, 3.0]]))

            with h5py.File(path, "r") as handle:
                validate_inputs.require_raw_count_matrix(handle, "query h5ad")


class ToolSelectionTests(unittest.TestCase):
    def test_parse_tools_normalizes_case_and_whitespace(self):
        self.assertEqual(
            validate_inputs.parse_tools(" CellTypist , panhumanpy "),
            {"celltypist", "panhumanpy"},
        )

    def test_parse_tools_rejects_unknown_names(self):
        with self.assertRaisesRegex(ValueError, "Unknown --tools entries: singler"):
            validate_inputs.parse_tools("celltypist,singler")

    def test_parse_tools_rejects_empty_selection(self):
        with self.assertRaisesRegex(ValueError, "must name at least one"):
            validate_inputs.parse_tools("  ,  ")

    def test_panhumanpy_only_selection_needs_no_reference(self):
        tools = validate_inputs.parse_tools("panhumanpy")

        self.assertFalse(validate_inputs.REFERENCE_TOOLS & tools)

    def test_celltypist_selection_needs_a_reference(self):
        tools = validate_inputs.parse_tools("celltypist,panhumanpy")

        self.assertTrue(validate_inputs.REFERENCE_TOOLS & tools)


class ValidateInputsCliTests(unittest.TestCase):
    def write_query(self, directory: Path) -> Path:
        path = directory / "query.h5ad"
        query = ad.AnnData(
            X=np.array([[1.0, 0.0], [0.0, 2.0]], dtype=np.float32),
            obs={"leiden": ["0", "1"]},
        )
        query.obs_names = ["cell-1", "cell-2"]
        query.var_names = ["GENE1", "GENE2"]
        query.write_h5ad(path)
        return path

    def run_cli(self, tmpdir: Path, *extra_args: str) -> tuple[dict, Path]:
        query_path = self.write_query(tmpdir)
        manifest_path = tmpdir / "validation_manifest.json"
        shared_genes_path = tmpdir / "shared_genes.txt"
        argv = [
            "validate_inputs.py",
            "--query-h5ad",
            str(query_path),
            "--query-cluster-col",
            "leiden",
            "--output",
            str(manifest_path),
            "--shared-genes-output",
            str(shared_genes_path),
            *extra_args,
        ]
        with unittest.mock.patch.object(sys, "argv", argv):
            validate_inputs.main()
        return json.loads(manifest_path.read_text(encoding="utf-8")), shared_genes_path

    def test_panhumanpy_only_run_skips_reference_and_shared_genes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest, shared_genes_path = self.run_cli(Path(tmpdir), "--tools", "panhumanpy")

        self.assertEqual(manifest["tools"], ["panhumanpy"])
        self.assertIsNone(manifest["ref_h5ad"])
        self.assertIsNone(manifest["ref_label_col"])
        self.assertIsNone(manifest["shared_gene_count"])
        self.assertEqual(manifest["query_cell_count"], 2)
        self.assertFalse(shared_genes_path.exists())

    def test_reference_based_tool_requires_ref_h5ad(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(ValueError, "--ref-h5ad is required"):
                self.run_cli(Path(tmpdir), "--tools", "celltypist")

    def test_reference_filtering_without_reference_tool_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(ValueError, "no reference-based tool is selected"):
                self.run_cli(
                    Path(tmpdir),
                    "--tools",
                    "panhumanpy",
                    "--ref-filter-col",
                    "tissue",
                    "--ref-filter-values",
                    "lung",
                )


if __name__ == "__main__":
    unittest.main()
