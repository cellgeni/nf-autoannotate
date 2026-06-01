import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.merge_annotations import (
    build_cluster_summary,
    celltypist_scanvi_agreement,
    read_annotation_csv,
    require_matching_cell_ids,
)


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

        summary = build_cluster_summary(merged, "leiden", "scanvi_label")

        self.assertIs(bool(summary.loc[0, "celltypist_scanvi_cluster_modal_match"]), True)
        self.assertEqual(summary.loc[0, "celltypist_scanvi_cluster_agreement_fraction"], 0.5)

    def test_require_matching_cell_ids_rejects_missing_or_extra_ids(self):
        frame = pd.DataFrame(index=pd.Index(["cell-1", "cell-extra"], name="cell_id"))

        with self.assertRaisesRegex(ValueError, "missing query cells: cell-2"):
            require_matching_cell_ids(frame, pd.Index(["cell-1", "cell-2"]), "CellTypist")

    def test_read_annotation_csv_rejects_duplicate_cell_ids(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            csv_path = Path(tmpdir) / "annotations.csv"
            csv_path.write_text("cell_id,label\ncell-1,A\ncell-1,B\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "duplicate cell_id"):
                read_annotation_csv(str(csv_path), "CellTypist")


if __name__ == "__main__":
    unittest.main()
