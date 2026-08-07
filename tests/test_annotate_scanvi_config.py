import argparse
import importlib.util
import sys
import tempfile
import types
import unittest
import unittest.mock
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "annotate_scanvi.py"

# The config helpers tested here do not need scvi-tools, but the script imports
# SCANVI at module import time. Stub it so this unit test stays lightweight.
scvi_module = types.ModuleType("scvi")
scvi_model_module = types.ModuleType("scvi.model")
class _StubSCANVI:
    """Stand-in for scvi.model.SCANVI; tests patch methods onto it."""


scvi_model_module.SCANVI = _StubSCANVI
scvi_module.model = scvi_model_module
sys.modules.setdefault("scvi", scvi_module)
sys.modules.setdefault("scvi.model", scvi_model_module)

spec = importlib.util.spec_from_file_location("annotate_scanvi", SCRIPT)
annotate_scanvi = importlib.util.module_from_spec(spec)
spec.loader.exec_module(annotate_scanvi)


def default_args():
    return argparse.Namespace(
        query_max_epochs=100,
        query_batch_size=None,
        query_train_size=None,
        query_validation_size=None,
        query_early_stopping=None,
        query_early_stopping_patience=None,
        query_check_val_every_n_epoch=10,
        query_plan_weight_decay=0.0,
    )


class QueryTrainConfigTests(unittest.TestCase):
    def test_single_item_lists_are_unwrapped_for_nf_scvi_metrics_compatibility(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.py"
            config_path.write_text(
                "query_train_input = {\n"
                "    \"max_epochs\": [20],\n"
                "    \"batch_size\": [64],\n"
                "    \"plan_weight_decay\": [0.01],\n"
                "}\n"
            )

            config = annotate_scanvi.load_query_train_input(str(config_path), default_args())

        self.assertEqual(config["max_epochs"], 20)
        self.assertEqual(config["batch_size"], 64)
        self.assertEqual(config["plan_kwargs"], {"weight_decay": 0.01})

    def test_sweep_lists_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.py"
            config_path.write_text("query_train_input = {\"max_epochs\": [10, 20]}\n")

            with self.assertRaisesRegex(ValueError, "not a sweep with 2 values"):
                annotate_scanvi.load_query_train_input(str(config_path), default_args())

    def test_missing_query_train_input_uses_nf_autoannotate_defaults(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.py"
            config_path.write_text("model_input = {}\n")

            config = annotate_scanvi.load_query_train_input(str(config_path), default_args())

        self.assertEqual(config["max_epochs"], 100)
        self.assertEqual(config["check_val_every_n_epoch"], 10)
        self.assertEqual(config["plan_kwargs"], {"weight_decay": 0.0})


class ModelGeneListTests(unittest.TestCase):
    def test_gene_list_comes_from_pipeline_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model_dir = Path(tmpdir)
            (model_dir / "scanvi_model_metadata.json").write_text(
                '{"genes": ["GENE1", "GENE2"]}', encoding="utf-8"
            )

            self.assertEqual(
                annotate_scanvi.load_gene_list(str(model_dir), None), ["GENE1", "GENE2"]
            )

    def test_external_model_falls_back_to_the_scvi_tools_registry(self):
        # A reference gene list is not a valid stand-in, because the model registry is a
        # post-HVG subset. scvi-tools reports the registry itself, so use that instead.
        calls = []

        def fake_prepare(query, model_dir, return_reference_var_names=False):
            calls.append((query, model_dir, return_reference_var_names))
            return ["GENE7", "GENE8"]

        with tempfile.TemporaryDirectory() as tmpdir:
            with unittest.mock.patch.object(
                annotate_scanvi.SCANVI, "prepare_query_anndata", fake_prepare, create=True
            ):
                genes = annotate_scanvi.load_gene_list(tmpdir, "query-object")

        self.assertEqual(genes, ["GENE7", "GENE8"])
        self.assertEqual(calls, [("query-object", tmpdir, True)])

    def test_unreadable_model_registry_is_reported_with_context(self):
        def fake_prepare(query, model_dir, return_reference_var_names=False):
            raise RuntimeError("not a model directory")

        with tempfile.TemporaryDirectory() as tmpdir:
            with unittest.mock.patch.object(
                annotate_scanvi.SCANVI, "prepare_query_anndata", fake_prepare, create=True
            ):
                with self.assertRaisesRegex(ValueError, "Could not determine the gene space"):
                    annotate_scanvi.load_gene_list(tmpdir, None)

    def test_metadata_without_a_gene_list_falls_through_to_the_registry(self):
        def fake_prepare(query, model_dir, return_reference_var_names=False):
            return ["GENE9"]

        with tempfile.TemporaryDirectory() as tmpdir:
            model_dir = Path(tmpdir)
            (model_dir / "scanvi_model_metadata.json").write_text("{}", encoding="utf-8")

            with unittest.mock.patch.object(
                annotate_scanvi.SCANVI, "prepare_query_anndata", fake_prepare, create=True
            ):
                self.assertEqual(annotate_scanvi.load_gene_list(str(model_dir), None), ["GENE9"])


if __name__ == "__main__":
    unittest.main()
