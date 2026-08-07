"""End-to-end workflow tests with stubbed annotation steps.

These cover the channel plumbing that unit tests cannot reach: sample-ID resolution,
per-query fan-out, grouping annotations back per sample, and the guarantee that every
resolved query reaches MERGE_ANNOTATIONS. A regression here loses a whole sample while
every Python unit test keeps passing, so the assertions compare the exact set of sample
IDs rather than just counting outputs.

Skipped when nextflow is not on PATH. Run with `module load cellgen/nextflow` first.
"""

import importlib.util
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[1]
STUB_SCRIPT = REPO / "tests" / "stubs" / "annotation_stubs.py"
spec = importlib.util.spec_from_file_location("annotation_stubs", STUB_SCRIPT)
annotation_stubs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(annotation_stubs)

REFERENCE_GENES = ["GENE1", "GENE2", "GENE3", "GENE4"]

LOCAL_CONFIG = """
executor { name = 'local' }
singularity.enabled = false
process { container = null }
"""


def write_query(path: Path, sample_id: str, genes=REFERENCE_GENES) -> None:
    cells = [f"{sample_id}-cell-{index}" for index in range(6)]
    query = ad.AnnData(
        X=np.arange(6 * len(genes), dtype=np.float32).reshape(6, len(genes)),
        obs=pd.DataFrame({"leiden": ["0", "0", "0", "1", "1", "1"]}, index=cells),
    )
    query.var_names = genes
    query.write_h5ad(path)


def write_reference(path: Path) -> None:
    cells = [f"ref-{index}" for index in range(10)]
    reference = ad.AnnData(
        X=np.arange(10 * len(REFERENCE_GENES), dtype=np.float32).reshape(10, len(REFERENCE_GENES)),
        obs=pd.DataFrame({"cell_type": ["A"] * 5 + ["B"] * 5}, index=cells),
    )
    reference.var_names = REFERENCE_GENES
    reference.write_h5ad(path)


@unittest.skipIf(shutil.which("nextflow") is None, "nextflow is not on PATH")
class MultiQueryPipelineTests(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tmpdir.name)
        (self.tmp / "queries").mkdir()
        (self.tmp / "local.config").write_text(LOCAL_CONFIG, encoding="utf-8")

        # The stubbed pipeline needs the real scripts/ next to its main.nf.
        self.pipeline = self.tmp / "pipeline"
        self.pipeline.mkdir()
        os.symlink(REPO / "scripts", self.pipeline / "scripts")
        for name in ("nextflow.config", "nextflow_schema.json"):
            shutil.copy(REPO / name, self.pipeline / name)
        (self.pipeline / "main.nf").write_text(
            annotation_stubs.stub_main_nf((REPO / "main.nf").read_text(encoding="utf-8")),
            encoding="utf-8",
        )

        write_reference(self.tmp / "reference.h5ad")

    def tearDown(self):
        self.tmpdir.cleanup()

    def run_pipeline(self, *args: str) -> Path:
        outdir = self.tmp / "results"
        command = [
            "nextflow", "run", str(self.pipeline / "main.nf"),
            "-c", str(self.tmp / "local.config"),
            "--project_tag", "test",
            "--outdir", str(outdir),
            *args,
        ]
        completed = subprocess.run(
            command, cwd=self.tmp, capture_output=True, text=True, timeout=900
        )
        if completed.returncode != 0:
            self.fail(f"pipeline failed:\n{completed.stdout}\n{completed.stderr}")
        return outdir

    def merged_sample_ids(self, outdir: Path) -> set[str]:
        return {path.name.split(".")[1] for path in outdir.glob("*/test.*.annotated.h5ad")}

    def test_every_query_reaches_merge_for_each_tool_subset(self):
        expected = {"donor1", "donor2", "donor3"}
        for sample_id in sorted(expected):
            write_query(self.tmp / "queries" / f"{sample_id}.h5ad", sample_id)

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
            write_query(self.tmp / "queries" / f"{sample_id}.h5ad", sample_id)

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
        write_query(self.tmp / "queries" / "donor1.h5ad", "donor1")
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


if __name__ == "__main__":
    unittest.main()
