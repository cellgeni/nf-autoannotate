"""End-to-end workflow tests with stubbed annotation steps.

These cover the channel plumbing that unit tests cannot reach: sample-ID resolution,
per-query fan-out, grouping annotations back per sample, and the guarantee that every
resolved query reaches MERGE_ANNOTATIONS. A regression here loses a whole sample while
every Python unit test keeps passing, so the assertions compare the exact set of sample
IDs rather than just counting outputs.

The output paths asserted here are the pipeline's published contract, so they are pinned
deliberately: the model-sheet workflow shares these processes, and its own naming must not
move the default workflow's results.

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
class MultiQueryPipelineTests(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tmpdir.name)
        (self.tmp / "queries").mkdir()
        self.pipeline = stub_pipeline.build_stub_pipeline(self.tmp)
        stub_pipeline.write_reference(self.tmp / "reference.h5ad")

    def tearDown(self):
        self.tmpdir.cleanup()

    def run_pipeline(self, *args: str) -> Path:
        outdir = self.tmp / "results"
        completed = stub_pipeline.run_pipeline(self.pipeline, self.tmp, *args, outdir=outdir)
        if completed.returncode != 0:
            self.fail(f"pipeline failed:\n{completed.stdout}\n{completed.stderr}")
        return outdir

    def merged_sample_ids(self, outdir: Path) -> set[str]:
        return {path.name.split(".")[1] for path in outdir.glob("*/test.*.annotated.h5ad")}

    def test_every_query_reaches_merge_for_each_tool_subset(self):
        expected = {"donor1", "donor2", "donor3"}
        for sample_id in sorted(expected):
            stub_pipeline.write_query(self.tmp / "queries" / f"{sample_id}.h5ad", sample_id)

        for tools in ("celltypist,scanvi,panhumanpy", "celltypist", "panhumanpy", "scanvi,panhumanpy"):
            with self.subTest(tools=tools):
                reference_args = []
                if "celltypist" in tools or "scanvi" in tools:
                    reference_args = [
                        "--ref_h5ad", str(self.tmp / "reference.h5ad"),
                        "--ref_label_col", "cell_type",
                    ]
                shutil.rmtree(self.tmp / "results", ignore_errors=True)
                outdir = self.run_pipeline(
                    "--tools", tools,
                    "--query_h5ad", str(self.tmp / "queries" / "*.h5ad"),
                    *reference_args,
                )

                self.assertEqual(self.merged_sample_ids(outdir), expected)

    def test_each_output_contains_only_its_own_cells(self):
        expected = {"donor1", "donor2"}
        for sample_id in sorted(expected):
            stub_pipeline.write_query(self.tmp / "queries" / f"{sample_id}.h5ad", sample_id)

        outdir = self.run_pipeline(
            "--tools", "panhumanpy",
            "--query_h5ad", str(self.tmp / "queries" / "*.h5ad"),
        )

        self.assertEqual(self.merged_sample_ids(outdir), expected)
        for sample_id in expected:
            annotated = ad.read_h5ad(outdir / sample_id / f"test.{sample_id}.annotated.h5ad")
            owners = {str(name).split("-")[0] for name in annotated.obs_names}
            self.assertEqual(owners, {sample_id})

    def test_pretrained_models_do_not_require_a_reference(self):
        stub_pipeline.write_query(self.tmp / "queries" / "donor1.h5ad", "donor1")
        # A stub CellTypist model is enough: the annotate step is stubbed too, so this
        # only checks that the workflow resolves without --ref_h5ad.
        model = self.tmp / "celltypist_model.pkl"
        model.write_text("stub", encoding="utf-8")

        outdir = self.run_pipeline(
            "--tools", "celltypist",
            "--query_h5ad", str(self.tmp / "queries" / "donor1.h5ad"),
            "--celltypist_model", str(model),
        )

        self.assertEqual(self.merged_sample_ids(outdir), {"donor1"})
        self.assertFalse((outdir / "reference_genes.txt").exists())

    def test_predictions_are_published_alongside_the_merged_output(self):
        stub_pipeline.write_query(self.tmp / "queries" / "donor1.h5ad", "donor1")

        outdir = self.run_pipeline(
            "--tools", "panhumanpy",
            "--query_h5ad", str(self.tmp / "queries" / "donor1.h5ad"),
        )

        self.assertTrue((outdir / "donor1" / "test.donor1.annotated.h5ad").exists())
        self.assertTrue((outdir / "donor1" / "test.donor1.panhumanpy_predictions.csv").exists())

    def test_merge_can_be_skipped(self):
        expected = {"donor1", "donor2"}
        for sample_id in sorted(expected):
            stub_pipeline.write_query(self.tmp / "queries" / f"{sample_id}.h5ad", sample_id)

        outdir = self.run_pipeline(
            "--tools", "scanvi,panhumanpy",
            "--query_h5ad", str(self.tmp / "queries" / "*.h5ad"),
            "--ref_h5ad", str(self.tmp / "reference.h5ad"),
            "--ref_label_col", "cell_type",
            "--merge_annotations", "false",
        )

        self.assertEqual(self.merged_sample_ids(outdir), set())
        for sample_id in sorted(expected):
            published = {path.name for path in (outdir / sample_id).iterdir()}
            self.assertEqual(
                published,
                {
                    f"{sample_id}.query_manifest.json",
                    f"test.{sample_id}.scanvi_predictions.csv",
                    f"test.{sample_id}.scanvi_latent.h5ad",
                    f"test.{sample_id}.panhumanpy_predictions.csv",
                },
            )

    def test_skipping_merge_does_not_require_a_cluster_column(self):
        # The cluster column is only read by the merge step, so a run that skips merging
        # must not reject a query for lacking it.
        stub_pipeline.write_query(
            self.tmp / "queries" / "donor1.h5ad", "donor1", with_cluster=False
        )

        outdir = self.run_pipeline(
            "--tools", "panhumanpy",
            "--query_h5ad", str(self.tmp / "queries" / "donor1.h5ad"),
            "--merge_annotations", "false",
        )

        self.assertTrue((outdir / "donor1" / "test.donor1.panhumanpy_predictions.csv").exists())


if __name__ == "__main__":
    unittest.main()
