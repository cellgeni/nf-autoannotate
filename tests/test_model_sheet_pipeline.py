"""Workflow tests for the ANNOTATE_WITH_MODELS entry point, with stubbed annotation steps.

The sheet is the only place this entry point's tools, models, and queries are named, so the
tests here cover the parsing and fan-out that follow from it: one model over several
queries, several models over one query, mixed tools in one sheet, TSV input, and the
per-pair output naming that keeps two models over the same query from overwriting each
other. The rejection cases matter just as much, because a sheet mistake that is not caught
during parsing turns into hours of wasted annotation.

Skipped when nextflow is not on PATH. Run with `module load cellgen/nextflow` first.
"""

import importlib.util
import shutil
import tempfile
import unittest
from pathlib import Path

import anndata as ad


REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "stub_pipeline", REPO / "tests" / "stubs" / "stub_pipeline.py"
)
stub_pipeline = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stub_pipeline)


@unittest.skipIf(shutil.which("nextflow") is None, "nextflow is not on PATH")
class ModelSheetPipelineTests(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tmpdir.name)
        (self.tmp / "queries").mkdir()
        (self.tmp / "models").mkdir()
        self.pipeline = stub_pipeline.build_stub_pipeline(self.tmp)

    def tearDown(self):
        self.tmpdir.cleanup()

    def write_queries(self, *sample_ids: str, with_cluster: bool = True) -> None:
        for sample_id in sample_ids:
            stub_pipeline.write_query(
                self.tmp / "queries" / f"{sample_id}.h5ad", sample_id, with_cluster=with_cluster
            )

    def celltypist_model(self, name: str) -> Path:
        model = self.tmp / "models" / name
        model.parent.mkdir(parents=True, exist_ok=True)
        model.write_text("stub", encoding="utf-8")
        return model

    def scanvi_model(self, name: str) -> Path:
        model = self.tmp / "models" / name
        model.mkdir(parents=True, exist_ok=True)
        (model / "scanvi_model_metadata.json").write_text("{}", encoding="utf-8")
        return model

    def write_sheet(self, rows: str, name: str = "model_sheet.csv") -> Path:
        sheet = self.tmp / name
        sheet.write_text(rows, encoding="utf-8")
        return sheet

    def run_sheet(self, sheet: Path, *args: str) -> Path:
        outdir = self.tmp / "results"
        completed = stub_pipeline.run_pipeline(
            self.pipeline,
            self.tmp,
            "-entry", "ANNOTATE_WITH_MODELS",
            "--model_sheet", str(sheet),
            *args,
            outdir=outdir,
        )
        if completed.returncode != 0:
            self.fail(f"pipeline failed:\n{completed.stdout}\n{completed.stderr}")
        return outdir

    def run_sheet_expecting_failure(self, sheet: Path, *args: str) -> str:
        completed = stub_pipeline.run_pipeline(
            self.pipeline,
            self.tmp,
            "-entry", "ANNOTATE_WITH_MODELS",
            "--model_sheet", str(sheet),
            *args,
            outdir=self.tmp / "results",
        )
        if completed.returncode == 0:
            self.fail(f"pipeline was expected to fail but succeeded:\n{completed.stdout}")
        return f"{completed.stdout}\n{completed.stderr}"

    def unit_dirs(self, outdir: Path) -> set[str]:
        return {path.name for path in outdir.iterdir() if path.is_dir()}

    def test_one_model_over_two_queries_and_another_over_a_third(self):
        self.write_queries("obj1", "obj2", "obj3")
        first = self.celltypist_model("run-a/celltypist_model.pkl")
        second = self.celltypist_model("run-b/celltypist_model.pkl")
        sheet = self.write_sheet(
            "tool,model,model_id,query_h5ad\n"
            f"celltypist,{first},model1,{self.tmp / 'queries' / 'obj1.h5ad'}\n"
            f"celltypist,{first},model1,{self.tmp / 'queries' / 'obj2.h5ad'}\n"
            f"celltypist,{second},model2,{self.tmp / 'queries' / 'obj3.h5ad'}\n"
        )

        outdir = self.run_sheet(sheet)

        self.assertEqual(self.unit_dirs(outdir), {"obj1.model1", "obj2.model1", "obj3.model2"})
        for unit_id in ("obj1.model1", "obj2.model1", "obj3.model2"):
            self.assertTrue((outdir / unit_id / f"test.{unit_id}.celltypist_predictions.csv").exists())
            # Merging is opt-in for this entry point.
            self.assertEqual(list(outdir.glob(f"{unit_id}/*.annotated.h5ad")), [])

    def test_two_models_over_one_query_keep_separate_outputs(self):
        self.write_queries("obj1")
        celltypist = self.celltypist_model("run-a/celltypist_model.pkl")
        scanvi = self.scanvi_model("run-b/scanvi_model")
        sheet = self.write_sheet(
            "tool,model,model_id,query_h5ad\n"
            f"celltypist,{celltypist},gut_ct,{self.tmp / 'queries' / 'obj1.h5ad'}\n"
            f"scanvi,{scanvi},lung_scanvi,{self.tmp / 'queries' / 'obj1.h5ad'}\n"
        )

        outdir = self.run_sheet(sheet)

        self.assertEqual(self.unit_dirs(outdir), {"obj1.gut_ct", "obj1.lung_scanvi"})
        self.assertTrue(
            (outdir / "obj1.gut_ct" / "test.obj1.gut_ct.celltypist_predictions.csv").exists()
        )
        self.assertTrue(
            (outdir / "obj1.lung_scanvi" / "test.obj1.lung_scanvi.scanvi_predictions.csv").exists()
        )
        self.assertTrue(
            (outdir / "obj1.lung_scanvi" / "test.obj1.lung_scanvi.scanvi_latent.h5ad").exists()
        )

    def test_tsv_sheet_with_a_glob_and_a_panhumanpy_row(self):
        self.write_queries("obj1", "obj2")
        celltypist = self.celltypist_model("run-a/celltypist_model.pkl")
        sheet = self.write_sheet(
            "tool\tmodel\tmodel_id\tquery_h5ad\n"
            f"celltypist\t{celltypist}\tgut_ct\t{self.tmp / 'queries' / '*.h5ad'}\n"
            f"panhumanpy\t\t\t{self.tmp / 'queries' / 'obj1.h5ad'}\n",
            name="model_sheet.tsv",
        )

        outdir = self.run_sheet(sheet)

        self.assertEqual(
            self.unit_dirs(outdir), {"obj1.gut_ct", "obj2.gut_ct", "obj1.panhumanpy"}
        )
        self.assertTrue(
            (outdir / "obj1.panhumanpy" / "test.obj1.panhumanpy.panhumanpy_predictions.csv").exists()
        )

    def test_merging_writes_one_annotated_object_per_pair(self):
        self.write_queries("obj1")
        first = self.celltypist_model("run-a/celltypist_model.pkl")
        second = self.celltypist_model("run-b/celltypist_model.pkl")
        sheet = self.write_sheet(
            "tool,model,model_id,query_h5ad\n"
            f"celltypist,{first},model1,{self.tmp / 'queries' / 'obj1.h5ad'}\n"
            f"celltypist,{second},model2,{self.tmp / 'queries' / 'obj1.h5ad'}\n"
        )

        outdir = self.run_sheet(sheet, "--merge_annotations", "true")

        for unit_id in ("obj1.model1", "obj1.model2"):
            annotated = outdir / unit_id / f"test.{unit_id}.annotated.h5ad"
            self.assertTrue(annotated.exists())
            self.assertIn("celltypist_majority_voting", ad.read_h5ad(annotated).obs.columns)
            self.assertTrue(
                (outdir / unit_id / f"test.{unit_id}.cluster_annotation_summary.csv").exists()
            )

    def test_model_id_defaults_to_the_model_name(self):
        self.write_queries("obj1")
        model = self.celltypist_model("run-a/celltypist_model.pkl")
        sheet = self.write_sheet(
            "tool,model,query_h5ad\n"
            f"celltypist,{model},{self.tmp / 'queries' / 'obj1.h5ad'}\n"
        )

        outdir = self.run_sheet(sheet)

        self.assertEqual(self.unit_dirs(outdir), {"obj1.celltypist_model"})

    def test_colliding_model_ids_are_rejected(self):
        # Two training runs both produce celltypist_model.pkl, so the default model ID
        # collides and the two pairs would publish over each other.
        self.write_queries("obj1")
        first = self.celltypist_model("run-a/celltypist_model.pkl")
        second = self.celltypist_model("run-b/celltypist_model.pkl")
        sheet = self.write_sheet(
            "tool,model,query_h5ad\n"
            f"celltypist,{first},{self.tmp / 'queries' / 'obj1.h5ad'}\n"
            f"celltypist,{second},{self.tmp / 'queries' / 'obj1.h5ad'}\n"
        )

        output = self.run_sheet_expecting_failure(sheet)

        self.assertIn("distinct model_id", output)

    def test_a_scanvi_row_pointing_at_a_file_is_rejected(self):
        self.write_queries("obj1")
        model = self.celltypist_model("run-a/celltypist_model.pkl")
        sheet = self.write_sheet(
            "tool,model,query_h5ad\n"
            f"scanvi,{model},{self.tmp / 'queries' / 'obj1.h5ad'}\n"
        )

        output = self.run_sheet_expecting_failure(sheet)

        self.assertIn("scanvi_model directory", output)

    def test_a_celltypist_row_without_a_model_is_rejected(self):
        self.write_queries("obj1")
        sheet = self.write_sheet(
            "tool,model,query_h5ad\n"
            f"celltypist,,{self.tmp / 'queries' / 'obj1.h5ad'}\n"
        )

        output = self.run_sheet_expecting_failure(sheet)

        self.assertIn("no model path", output)

    def test_an_unknown_tool_is_rejected(self):
        self.write_queries("obj1")
        sheet = self.write_sheet(
            "tool,model,query_h5ad\n"
            f"scarches,,{self.tmp / 'queries' / 'obj1.h5ad'}\n"
        )

        output = self.run_sheet_expecting_failure(sheet)

        self.assertIn("unknown tool", output)

    def test_globally_supplied_queries_are_rejected(self):
        self.write_queries("obj1")
        model = self.celltypist_model("run-a/celltypist_model.pkl")
        sheet = self.write_sheet(
            "tool,model,query_h5ad\n"
            f"celltypist,{model},{self.tmp / 'queries' / 'obj1.h5ad'}\n"
        )

        output = self.run_sheet_expecting_failure(
            sheet, "--query_h5ad", str(self.tmp / "queries" / "obj1.h5ad")
        )

        self.assertIn("--query_h5ad", output)

    def test_the_sheet_is_rejected_by_the_default_workflow(self):
        self.write_queries("obj1")
        model = self.celltypist_model("run-a/celltypist_model.pkl")
        sheet = self.write_sheet(
            "tool,model,query_h5ad\n"
            f"celltypist,{model},{self.tmp / 'queries' / 'obj1.h5ad'}\n"
        )

        completed = stub_pipeline.run_pipeline(
            self.pipeline,
            self.tmp,
            "--model_sheet", str(sheet),
            "--tools", "panhumanpy",
            "--query_h5ad", str(self.tmp / "queries" / "obj1.h5ad"),
            outdir=self.tmp / "results",
        )

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("ANNOTATE_WITH_MODELS", f"{completed.stdout}\n{completed.stderr}")

    def test_a_query_without_a_cluster_column_is_accepted_without_merging(self):
        self.write_queries("obj1", with_cluster=False)
        model = self.celltypist_model("run-a/celltypist_model.pkl")
        sheet = self.write_sheet(
            "tool,model,query_h5ad\n"
            f"celltypist,{model},{self.tmp / 'queries' / 'obj1.h5ad'}\n"
        )

        outdir = self.run_sheet(sheet)

        self.assertTrue(
            (outdir / "obj1.celltypist_model"
             / "test.obj1.celltypist_model.celltypist_predictions.csv").exists()
        )


if __name__ == "__main__":
    unittest.main()
