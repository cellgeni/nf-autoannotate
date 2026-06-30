import importlib.util
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import h5py
import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "validate_inputs.py"
spec = importlib.util.spec_from_file_location("validate_inputs", SCRIPT)
validate_inputs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validate_inputs)


class ValidateInputsTests(unittest.TestCase):
    def test_query_cluster_column_is_required_when_missing_cluster_computation_is_disabled(self):
        query = ad.AnnData(obs={"cell_type": ["A", "B"]})

        with self.assertRaisesRegex(ValueError, "Missing required obs column 'leiden' in query h5ad"):
            validate_inputs.require_query_cluster_column(
                query, "leiden", compute_missing_clusters=False
            )

    def test_query_cluster_column_can_be_missing_when_cluster_computation_is_enabled(self):
        query = ad.AnnData(obs={"cell_type": ["A", "B"]})

        validate_inputs.require_query_cluster_column(query, "leiden", compute_missing_clusters=True)

    def test_raw_count_matrix_rejects_fractional_values(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "query.h5ad"
            with h5py.File(path, "w") as handle:
                handle.create_dataset("X", data=np.array([[0.0, 1.5], [2.0, 3.0]]))

            with h5py.File(path, "r") as handle:
                with self.assertRaisesRegex(ValueError, "does not look like raw counts"):
                    validate_inputs.require_raw_count_matrix(handle, "query h5ad")

    def test_raw_count_matrix_accepts_integer_float_values(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "query.h5ad"
            with h5py.File(path, "w") as handle:
                handle.create_dataset("X", data=np.array([[0.0, 1.0], [2.0, 3.0]]))

            with h5py.File(path, "r") as handle:
                validate_inputs.require_raw_count_matrix(handle, "query h5ad")


if __name__ == "__main__":
    unittest.main()
