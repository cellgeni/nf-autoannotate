import argparse
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from scripts.merge_annotations import (
    build_cluster_summary,
    celltypist_scanvi_agreement,
    read_annotation_csv,
    require_matching_cell_ids,
    resolve_tools,
)
from scripts.merge_annotations import main as merge_main


ALL_TOOLS = {"celltypist", "scanvi", "panhumanpy"}


def make_args(**overrides):
    defaults = {
        "celltypist_csv": None,
        "scanvi_csv": None,
        "panhumanpy_csv": None,
        "scanvi_latent_h5ad": None,
    }
    return argparse.Namespace(**{**defaults, **overrides})


class MergeAnnotationsTests(unittest.TestCase):
    def test_celltypist_scanvi_agreement_requires_two_labels(self):
        merged = pd.DataFrame(
            {
                "celltypist_majority_voting": ["A", "A", None, None],
                "scanvi_label": ["A", "B", "A", None],
            }
        )

        observed = celltypist_scanvi_agreement(merged, "scanvi_label")

        self.assertEqual(observed.tolist(), [True, False, False, False])

    def test_cluster_modal_match_is_distinct_from_cell_agreement_fraction(self):
        merged = pd.DataFrame(
            {
                "leiden": ["0", "0", "0", "0"],
                "celltypist_majority_voting": ["A", "A", "B", "B"],
                "scanvi_label": ["A", "B", "A", "B"],
                "panhumanpy_final_level_label": ["P", "P", "Q", "Q"],
            }
        )
        merged["celltypist_scanvi_agree"] = celltypist_scanvi_agreement(merged, "scanvi_label")

        summary = build_cluster_summary(merged, "leiden", "scanvi_label", ALL_TOOLS)

        self.assertIs(bool(summary.loc[0, "celltypist_scanvi_cluster_modal_match"]), True)
        self.assertEqual(summary.loc[0, "celltypist_scanvi_cluster_agreement_fraction"], 0.5)

    def test_require_matching_cell_ids_rejects_missing_or_extra_ids(self):
        frame = pd.DataFrame(index=pd.Index(["cell-1", "cell-extra"], name="cell_id"))

        with self.assertRaisesRegex(ValueError, "missing query cells: cell-2"):
            require_matching_cell_ids(frame, pd.Index(["cell-1", "cell-2"]), "CellTypist")

    def test_cluster_summary_only_reports_selected_tools(self):
        merged = pd.DataFrame(
            {
                "leiden": ["0", "0"],
                "panhumanpy_final_level_label": ["P", "Q"],
            }
        )

        summary = build_cluster_summary(merged, "leiden", "scanvi_label", {"panhumanpy"})

        self.assertEqual(
            summary.columns.tolist(),
            ["leiden", "n_cells", "panhumanpy_top_label", "panhumanpy_top_fraction"],
        )

    def test_cluster_summary_omits_agreement_columns_without_both_tools(self):
        merged = pd.DataFrame(
            {
                "leiden": ["0", "0"],
                "celltypist_majority_voting": ["A", "B"],
                "panhumanpy_final_level_label": ["P", "Q"],
            }
        )

        summary = build_cluster_summary(merged, "leiden", "scanvi_label", {"celltypist", "panhumanpy"})

        self.assertIn("celltypist_majority_top_label", summary.columns)
        self.assertNotIn("scanvi_top_label", summary.columns)
        self.assertNotIn("celltypist_scanvi_cluster_modal_match", summary.columns)
        self.assertNotIn("celltypist_scanvi_cluster_agreement_fraction", summary.columns)

    def test_resolve_tools_derives_selection_from_provided_csvs(self):
        self.assertEqual(
            resolve_tools(make_args(celltypist_csv="c.csv", panhumanpy_csv="p.csv")),
            {"celltypist", "panhumanpy"},
        )

    def test_resolve_tools_requires_at_least_one_csv(self):
        with self.assertRaisesRegex(ValueError, "at least one of"):
            resolve_tools(make_args())

    def test_resolve_tools_rejects_scanvi_latent_without_scanvi_predictions(self):
        with self.assertRaisesRegex(ValueError, "requires --scanvi-csv"):
            resolve_tools(make_args(celltypist_csv="c.csv", scanvi_latent_h5ad="latent.h5ad"))

    def test_read_annotation_csv_rejects_duplicate_cell_ids(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_path = Path(tmpdir) / "annotations.csv"
            csv_path.write_text("cell_id,label\ncell-1,A\ncell-1,B\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "duplicate cell_id"):
                read_annotation_csv(str(csv_path), "CellTypist")


class MergeAnnotationsSingleToolTests(unittest.TestCase):
    def test_panhumanpy_only_run_produces_outputs_without_other_tools(self):
        cell_ids = [f"cell-{index}" for index in range(6)]
        query = ad.AnnData(
            X=np.arange(24, dtype=np.float32).reshape(6, 4),
            obs=pd.DataFrame({"leiden": ["0", "0", "0", "1", "1", "1"]}, index=cell_ids),
        )
        query.var_names = ["GENE1", "GENE2", "GENE3", "GENE4"]

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            query_path = tmp / "query.h5ad"
            query.write_h5ad(query_path)

            panhumanpy_path = tmp / "panhumanpy_predictions.csv"
            pd.DataFrame(
                {
                    "cell_id": cell_ids,
                    "panhumanpy_full_hierarchical_label": ["A|B"] * 6,
                    "panhumanpy_level_zero_label": ["A"] * 6,
                    "panhumanpy_final_level_label": ["B"] * 3 + ["C"] * 3,
                    "panhumanpy_confidence": [0.9] * 6,
                }
            ).to_csv(panhumanpy_path, index=False)

            output_path = tmp / "annotated.h5ad"
            summary_path = tmp / "summary.csv"
            markers_path = tmp / "markers.tsv"
            argv = [
                "merge_annotations.py",
                "--query-h5ad", str(query_path),
                "--panhumanpy-csv", str(panhumanpy_path),
                "--cluster-column", "leiden",
                "--scanvi-label-column", "scanvi_label",
                "--scanvi-score-column", "scanvi_confidence",
                "--output", str(output_path),
                "--cluster-summary-output", str(summary_path),
                "--marker-genes-output", str(markers_path),
                "--marker-genes-method", "t-test",
            ]
            with unittest.mock.patch.object(sys, "argv", argv):
                merge_main()

            annotated = ad.read_h5ad(output_path)
            summary = pd.read_csv(summary_path)
            markers = pd.read_csv(markers_path, sep="\t")

        self.assertIn("panhumanpy_final_level_label", annotated.obs.columns)
        self.assertNotIn("celltypist_majority_voting", annotated.obs.columns)
        self.assertNotIn("celltypist_scanvi_agree", annotated.obs.columns)
        self.assertEqual(
            summary.columns.tolist(),
            ["leiden", "n_cells", "panhumanpy_top_label", "panhumanpy_top_fraction"],
        )
        self.assertEqual(sorted(markers["cluster"].astype(str).unique()), ["0", "1"])


if __name__ == "__main__":
    unittest.main()
