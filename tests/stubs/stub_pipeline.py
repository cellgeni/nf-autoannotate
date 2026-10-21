"""Fixture for the tests that run the real workflow with stubbed annotation steps.

Both pipeline-level test modules build the same thing: a copy of main.nf whose annotation
steps are faked, the real scripts/ next to it, a local-executor config, and tiny h5ad
inputs. That setup lives here so the two modules cannot drift apart, and so a change to the
pipeline's layout is made in one place.

Everything except the annotation steps is the real pipeline: parameter validation, sheet
parsing, sample-ID resolution, fan-out, grouping, joining, and publishing.
"""

import importlib.util
import os
import shutil
import subprocess
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[2]

_spec = importlib.util.spec_from_file_location(
    "annotation_stubs", Path(__file__).resolve().parent / "annotation_stubs.py"
)
annotation_stubs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(annotation_stubs)

REFERENCE_GENES = ["GENE1", "GENE2", "GENE3", "GENE4"]

LOCAL_CONFIG = """
executor { name = 'local' }
singularity.enabled = false
process { container = null }
"""


def write_query(path: Path, sample_id: str, genes=REFERENCE_GENES, with_cluster: bool = True) -> Path:
    cells = [f"{sample_id}-cell-{index}" for index in range(6)]
    obs = pd.DataFrame(index=cells)
    if with_cluster:
        obs["leiden"] = ["0", "0", "0", "1", "1", "1"]
    query = ad.AnnData(
        X=np.arange(6 * len(genes), dtype=np.float32).reshape(6, len(genes)),
        obs=obs,
    )
    query.var_names = genes
    query.write_h5ad(path)
    return path


def write_reference(path: Path) -> Path:
    cells = [f"ref-{index}" for index in range(10)]
    reference = ad.AnnData(
        X=np.arange(10 * len(REFERENCE_GENES), dtype=np.float32).reshape(10, len(REFERENCE_GENES)),
        obs=pd.DataFrame({"cell_type": ["A"] * 5 + ["B"] * 5}, index=cells),
    )
    reference.var_names = REFERENCE_GENES
    reference.write_h5ad(path)
    return path


def build_stub_pipeline(tmp: Path) -> Path:
    """Write a stubbed copy of the pipeline under tmp and return its directory."""
    pipeline = tmp / "pipeline"
    pipeline.mkdir()
    # The stubbed pipeline needs the real scripts/ next to its main.nf.
    os.symlink(REPO / "scripts", pipeline / "scripts")
    for name in ("nextflow.config", "nextflow_schema.json"):
        shutil.copy(REPO / name, pipeline / name)
    (pipeline / "main.nf").write_text(
        annotation_stubs.stub_main_nf((REPO / "main.nf").read_text(encoding="utf-8")),
        encoding="utf-8",
    )
    (tmp / "local.config").write_text(LOCAL_CONFIG, encoding="utf-8")
    return pipeline


def run_pipeline(
    pipeline: Path,
    tmp: Path,
    *args: str,
    project_tag: str = "test",
    outdir: Path | None = None,
    timeout: int = 900,
) -> subprocess.CompletedProcess:
    """Run the stubbed pipeline and return the completed process, failed or not."""
    outdir = outdir if outdir is not None else tmp / "results"
    command = [
        "nextflow", "run", str(pipeline / "main.nf"),
        "-c", str(tmp / "local.config"),
        *args,
        "--project_tag", project_tag,
        "--outdir", str(outdir),
    ]
    return subprocess.run(command, cwd=tmp, capture_output=True, text=True, timeout=timeout)
