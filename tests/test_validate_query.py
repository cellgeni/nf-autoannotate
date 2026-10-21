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
spec = importlib.util.spec_from_file_location("validate_query", SCRIPTS / "validate_query.py")
validate_query = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validate_query)


REFERENCE_GENES = ["GENE1", "GENE2", "GENE3", "GENE4"]


def write_query(path: Path, genes: list[str], with_cluster: bool = True) -> Path:
    cells = [f"cell-{index}" for index in range(4)]
    obs = pd.DataFrame({"leiden": ["0", "0", "1", "1"]}, index=cells) if with_cluster else pd.DataFrame(index=cells)
    query = ad.AnnData(
        X=np.arange(4 * len(genes), dtype=np.float32).reshape(4, len(genes)),
        obs=obs,
    )
    query.var_names = genes
    query.write_h5ad(path)
    return path


def write_reference_side(tmp: Path, genes: list[str] = REFERENCE_GENES) -> tuple[Path, Path]:
    genes_path = tmp / "reference_genes.txt"
    genes_path.write_text("\n".join(genes) + "\n", encoding="utf-8")
    manifest_path = tmp / "reference_manifest.json"
    manifest_path.write_text(json.dumps({"gene_id_type": "symbol"}), encoding="utf-8")
    return genes_path, manifest_path


class ValidateQueryTests(unittest.TestCase):
    def run_cli(self, tmp: Path, *extra_args: str) -> dict:
        output = tmp / "manifest.json"
        argv = [
            "validate_query.py",
            "--sample-id", "donor1",
            "--query-cluster-col", "leiden",
            "--output", str(output),
            *extra_args,
        ]
        with unittest.mock.patch.object(sys, "argv", argv):
            validate_query.main()
        return json.loads(output.read_text(encoding="utf-8"))

    def test_full_coverage_is_recorded(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query_path = write_query(tmp / "donor1.h5ad", REFERENCE_GENES)
            genes_path, manifest_path = write_reference_side(tmp)

            manifest = self.run_cli(
                tmp,
                "--query-h5ad", str(query_path),
                "--tools", "celltypist",
                "--reference-genes", str(genes_path),
                "--reference-manifest", str(manifest_path),
            )

        self.assertEqual(manifest["sample_id"], "donor1")
        self.assertEqual(manifest["reference_gene_count"], 4)
        self.assertEqual(manifest["reference_genes_covered"], 4)
        self.assertEqual(manifest["reference_gene_coverage_fraction"], 1.0)

    def test_partial_coverage_is_recorded(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query_path = write_query(tmp / "donor1.h5ad", ["GENE1", "GENE2"])
            genes_path, manifest_path = write_reference_side(tmp)

            manifest = self.run_cli(
                tmp,
                "--query-h5ad", str(query_path),
                "--tools", "celltypist",
                "--reference-genes", str(genes_path),
                "--reference-manifest", str(manifest_path),
            )

        self.assertEqual(manifest["reference_genes_covered"], 2)
        self.assertEqual(manifest["reference_gene_coverage_fraction"], 0.5)

    def test_query_sharing_no_genes_with_reference_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query_path = write_query(tmp / "donor1.h5ad", ["OTHER1", "OTHER2"])
            genes_path, manifest_path = write_reference_side(tmp)

            with self.assertRaisesRegex(ValueError, "shares no genes with the reference"):
                self.run_cli(
                    tmp,
                    "--query-h5ad", str(query_path),
                    "--tools", "celltypist",
                    "--reference-genes", str(genes_path),
                    "--reference-manifest", str(manifest_path),
                )

    def test_mismatched_gene_identifier_types_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query_path = write_query(tmp / "donor1.h5ad", ["ENSG00000141510", "ENSG00000012048"])
            genes_path, manifest_path = write_reference_side(tmp, ["ENSG00000141510", "GENE2"])
            # Reference manifest claims symbols while the query is clearly Ensembl.

            with self.assertRaisesRegex(ValueError, "different gene identifier types"):
                self.run_cli(
                    tmp,
                    "--query-h5ad", str(query_path),
                    "--tools", "celltypist",
                    "--reference-genes", str(genes_path),
                    "--reference-manifest", str(manifest_path),
                )

    def test_panhumanpy_only_run_needs_no_reference_inputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query_path = write_query(tmp / "donor1.h5ad", REFERENCE_GENES)

            manifest = self.run_cli(
                tmp,
                "--query-h5ad", str(query_path),
                "--tools", "panhumanpy",
            )

        self.assertEqual(manifest["tools"], ["panhumanpy"])
        self.assertIsNone(manifest["reference_gene_count"])
        self.assertIsNone(manifest["reference_gene_coverage_fraction"])
        self.assertEqual(manifest["query_cell_count"], 4)

    def test_missing_cluster_column_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query_path = write_query(tmp / "donor1.h5ad", REFERENCE_GENES, with_cluster=False)

            with self.assertRaisesRegex(ValueError, "Missing required obs column 'leiden'"):
                self.run_cli(tmp, "--query-h5ad", str(query_path), "--tools", "panhumanpy")

    def test_missing_cluster_column_is_allowed_when_computing_clusters(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query_path = write_query(tmp / "donor1.h5ad", REFERENCE_GENES, with_cluster=False)

            manifest = self.run_cli(
                tmp,
                "--query-h5ad", str(query_path),
                "--tools", "panhumanpy",
                "--compute-missing-clusters",
            )

        self.assertTrue(manifest["compute_missing_clusters"])

    def test_missing_cluster_column_is_allowed_when_the_merge_step_is_skipped(self):
        # Only the merge step reads the cluster column, so a predictions-only run must not
        # be blocked by its absence.
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query_path = write_query(tmp / "donor1.h5ad", REFERENCE_GENES, with_cluster=False)

            manifest = self.run_cli(
                tmp,
                "--query-h5ad", str(query_path),
                "--tools", "celltypist",
                "--skip-merge",
            )

        self.assertFalse(manifest["merge_annotations"])
        self.assertFalse(manifest["compute_missing_clusters"])

    def test_merge_is_recorded_in_the_manifest_by_default(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query_path = write_query(tmp / "donor1.h5ad", REFERENCE_GENES)

            manifest = self.run_cli(
                tmp,
                "--query-h5ad", str(query_path),
                "--tools", "celltypist",
            )

        self.assertTrue(manifest["merge_annotations"])

    def test_celltypist_only_accepts_normalised_x(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query_path = write_query(tmp / "donor1.h5ad", REFERENCE_GENES)
            query = ad.read_h5ad(query_path)
            counts = np.array([[5, 0, 15, 5], [1, 9, 0, 2], [0, 4, 4, 1], [3, 3, 3, 3]], dtype=np.float32)
            query.X = np.log1p(counts / counts.sum(axis=1, keepdims=True) * 1e4).astype(np.float32)
            query.write_h5ad(query_path)

            manifest = self.run_cli(tmp, "--query-h5ad", str(query_path), "--tools", "celltypist")

        self.assertEqual(manifest["sample_id"], "donor1")

    def test_normalised_x_is_rejected_when_another_tool_needs_counts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query_path = write_query(tmp / "donor1.h5ad", REFERENCE_GENES)
            query = ad.read_h5ad(query_path)
            query.X = np.full((4, 4), 2.5, dtype=np.float32)
            query.write_h5ad(query_path)

            with self.assertRaisesRegex(ValueError, "raw counts"):
                self.run_cli(tmp, "--query-h5ad", str(query_path), "--tools", "celltypist,scanvi")

    def test_celltypist_only_without_any_raw_counts_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query_path = write_query(tmp / "donor1.h5ad", REFERENCE_GENES)
            query = ad.read_h5ad(query_path)
            query.X = np.full((4, 4), 2.5, dtype=np.float32)
            query.write_h5ad(query_path)

            with self.assertRaisesRegex(ValueError, "no usable expression matrix for CellTypist"):
                self.run_cli(tmp, "--query-h5ad", str(query_path), "--tools", "celltypist")

    def test_scanvi_options_without_scanvi_selected_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query_path = write_query(tmp / "donor1.h5ad", REFERENCE_GENES)

            with self.assertRaisesRegex(ValueError, "'scanvi' is not among the selected tools"):
                self.run_cli(
                    tmp,
                    "--query-h5ad", str(query_path),
                    "--tools", "panhumanpy",
                    "--batch-key", "donor",
                )


if __name__ == "__main__":
    unittest.main()
