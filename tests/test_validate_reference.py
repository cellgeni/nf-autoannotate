import importlib.util
import json
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location("validate_reference", SCRIPTS / "validate_reference.py")
validate_reference = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validate_reference)


def write_reference(path: Path) -> Path:
    cells = [f"ref-{index}" for index in range(6)]
    reference = ad.AnnData(
        X=np.arange(24, dtype=np.float32).reshape(6, 4),
        obs=pd.DataFrame(
            {
                "cell_type": ["A", "A", "A", "B", "B", "B"],
                "tissue": ["lung", "lung", "gut", "gut", "lung", "gut"],
            },
            index=cells,
        ),
    )
    reference.var_names = ["GENE1", "GENE2", "GENE3", "GENE4"]
    reference.write_h5ad(path)
    return path


class ValidateReferenceTests(unittest.TestCase):
    def run_cli(self, tmp: Path, *extra_args: str) -> tuple[dict, list[str]]:
        reference_path = write_reference(tmp / "reference.h5ad")
        manifest_path = tmp / "reference_manifest.json"
        genes_path = tmp / "reference_genes.txt"
        argv = [
            "validate_reference.py",
            "--ref-h5ad", str(reference_path),
            "--ref-label-col", "cell_type",
            "--output", str(manifest_path),
            "--genes-output", str(genes_path),
            *extra_args,
        ]
        with unittest.mock.patch.object(sys, "argv", argv):
            validate_reference.main()
        genes = genes_path.read_text(encoding="utf-8").split() if genes_path.exists() else []
        return json.loads(manifest_path.read_text(encoding="utf-8")), genes

    def test_publishes_the_full_reference_gene_space(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest, genes = self.run_cli(Path(tmpdir))

        # Training is query-independent, so every reference gene is published.
        self.assertEqual(genes, ["GENE1", "GENE2", "GENE3", "GENE4"])
        self.assertEqual(manifest["reference_gene_count"], 4)
        self.assertEqual(manifest["gene_id_type"], "symbol")
        self.assertEqual(manifest["ref_cell_count"], 6)
        self.assertEqual(manifest["filtered_ref_cell_count"], 6)

    def test_reference_filtering_narrows_the_cell_count(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            manifest, _ = self.run_cli(
                Path(tmpdir), "--ref-filter-col", "tissue", "--ref-filter-values", "lung"
            )

        self.assertEqual(manifest["ref_filter_col"], "tissue")
        self.assertEqual(manifest["ref_filter_values"], ["lung"])
        self.assertEqual(manifest["ref_cell_count"], 6)
        self.assertEqual(manifest["filtered_ref_cell_count"], 3)

    def test_unknown_filter_values_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(ValueError, "were not found in 'tissue': brain"):
                self.run_cli(
                    Path(tmpdir), "--ref-filter-col", "tissue", "--ref-filter-values", "brain"
                )

    def test_missing_batch_key_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(ValueError, "Missing required obs column 'donor'"):
                self.run_cli(Path(tmpdir), "--batch-key", "donor")

    def test_filter_column_without_values_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(ValueError, "required together"):
                self.run_cli(Path(tmpdir), "--ref-filter-col", "tissue")


if __name__ == "__main__":
    unittest.main()
