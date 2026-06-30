import argparse
import importlib.util
import sys
import tempfile
import types
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "annotate_scanvi.py"

# The config helpers tested here do not need scvi-tools, but the script imports
# SCANVI at module import time. Stub it so this unit test stays lightweight.
scvi_module = types.ModuleType("scvi")
scvi_model_module = types.ModuleType("scvi.model")
scvi_model_module.SCANVI = object
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


if __name__ == "__main__":
    unittest.main()
