"""Predictor import and dispatch tests, with no model or scientific execution."""
import importlib.abc
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from CoREMOF import prediction


class BlockOptionalImports(importlib.abc.MetaPathFinder):
    def __init__(self):
        self.attempts = []

    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'PACMANCharge', 'keras', 'tensorflow', 'matminer', 'molSimplify'}:
            self.attempts.append(fullname)
            raise ModuleNotFoundError('Disabled optional backend: ' + fullname, name=fullname)
        return None


class PredictionDependencyTests(unittest.TestCase):
    def test_module_imports_with_only_standard_library(self):
        root = Path(__file__).resolve().parents[1]
        script = """
import sys
from CoREMOF import prediction
for name in ('cloudpickle', 'keras', 'tensorflow', 'numpy', 'pandas',
             'PACMANCharge', 'requests', 'molSimplify', 'matminer'):
    assert name not in sys.modules, name
for name in ('pacman', 'cp', 'stability', 'precision', 'recall', 'f1'):
    assert callable(getattr(prediction, name)), name
"""
        result = subprocess.run([sys.executable, '-B', '-S', '-c', script], cwd=root,
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_missing_heat_capacity_models_are_reported_before_featurization(self):
        with tempfile.TemporaryDirectory() as directory:
            blocker = BlockOptionalImports()
            with patch.object(prediction, 'package_directory', directory), patch.object(sys, 'meta_path', [blocker] + sys.meta_path):
                with self.assertRaisesRegex(FileNotFoundError, 'ensemble models are missing'):
                    prediction.cp(Path(directory) / 'source.cif', T=[300])
            self.assertEqual(blocker.attempts, [])

    def test_missing_pacman_backend_is_not_a_false_prediction_success(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / 'source.cif'
            source.write_text('non-scientific input fixture')
            with patch.dict(sys.modules, {'PACMANCharge': None}):
                with self.assertRaises(ModuleNotFoundError):
                    prediction.pacman(source, root / 'output')
            self.assertFalse((root / 'output').exists())
            self.assertEqual(source.read_text(), 'non-scientific input fixture')


if __name__ == '__main__':
    unittest.main()
