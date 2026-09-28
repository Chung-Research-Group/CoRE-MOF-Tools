"""No-science regression tests for upstream CIF mutation and incomplete models."""
import ast
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from CoREMOF import _mofclassifier as api


class BatchTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.inputs = self.root / 'cifs'
        self.inputs.mkdir()
        self.cif = self.inputs / 'example.cif'
        self.cif.write_text('data_original\n')
        self.output = self.root / 'result.json'
        self.mode = 'normal'

    def predict(self, *, root_cifs, model, batch_size):
        self.assertNotEqual(Path(root_cifs[0]), self.cif)
        self.assertEqual(Path(root_cifs[0]).read_text(), 'data_original\n')
        self.assertEqual((model, batch_size), ('core', 64))
        if self.mode == 'rewrite':
            Path(root_cifs[0]).write_text('data_parser_rewrite\n')
        if self.mode == 'exception':
            raise ValueError('parser failure')
        if self.mode == 'source_change':
            self.cif.write_text('changed concurrently\n')
        result = [['example', [0.8] * 100, 0.8]]
        if self.mode == 'nan':
            result[0][1][0] = float('nan')
        if self.mode == 'incomplete':
            result[0][1].pop()
        if self.mode == 'bad_mean':
            result[0][2] = 0.2
        if self.mode == 'unknown_id':
            result[0][0] = 'other'
        if self.mode == 'missing':
            result.clear()
        if self.mode == 'duplicate':
            result *= 2
        return result

    def run_api(self, **options):
        with patch.object(api, '_load_classifier', return_value=SimpleNamespace(predict_batch=self.predict)):
            return api.predict_directory(self.inputs, self.output, 'core', 64, **options)

    def test_original_shape_and_mean_preserved(self):
        result = self.run_api()
        self.assertEqual(result, {'example': [[0.8] * 100, 0.8]})
        self.assertEqual(json.loads(self.output.read_text()), result)
        self.assertEqual(self.cif.read_text(), 'data_original\n')

    def test_upstream_parser_rewrite_is_private_and_reported(self):
        self.mode = 'rewrite'
        with self.assertWarnsRegex(RuntimeWarning, 'private parser inputs'):
            self.run_api()
        self.assertEqual(self.cif.read_text(), 'data_original\n')

    def test_bad_results_are_never_published(self):
        for mode in ('nan', 'incomplete', 'bad_mean', 'unknown_id', 'missing', 'duplicate', 'exception'):
            self.mode = mode
            with self.subTest(mode=mode), self.assertRaises(ValueError):
                self.run_api()
            self.assertFalse(self.output.exists())

    def test_original_change_rejects_output(self):
        self.mode = 'source_change'
        with self.assertRaises(RuntimeError):
            self.run_api()
        self.assertFalse(self.output.exists())

    def test_existing_output_needs_explicit_overwrite(self):
        self.output.write_text('old result\n')
        with self.assertRaises(FileExistsError):
            self.run_api()
        self.assertEqual(self.output.read_text(), 'old result\n')
        self.run_api(overwrite=True)
        self.assertEqual(json.loads(self.output.read_text())['example'][1], 0.8)

    def test_output_cannot_replace_source(self):
        self.output = self.cif
        with self.assertRaises(ValueError):
            self.run_api(overwrite=True)
        self.assertEqual(self.cif.read_text(), 'data_original\n')

    def test_bad_options_fail_before_import(self):
        for model, size in (('unknown', 64), ('core', 0), ('core', True)):
            with self.subTest(model=model, size=size), self.assertRaises(ValueError):
                api.predict_directory(self.inputs, self.output, model, size)

    def test_import_of_curate_does_not_import_classifier(self):
        tree = ast.parse((Path(api.__file__).parent / 'curate.py').read_text())
        self.assertFalse(any(isinstance(node, ast.ImportFrom) and (node.module or '').startswith('MOFClassifier')
                             for node in ast.walk(tree)))
        function = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == 'run_mofclassifier')
        self.assertIn('overwrite', [arg.arg for arg in function.args.kwonlyargs])


class AssetTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for name in ('models', 'models_qsp', 'models_h'):
            (self.root / name).mkdir()
            (self.root / name / 'already-installed').write_text('present\n')
        (self.root / 'atom_init.json').write_text('{}')
        for index in range(1, 101):
            (self.root / 'models' / f'checkpoint_bag_{index}.pth.tar').write_text('model fixture, never loaded\n')

    def test_missing_last_model_prevents_import_and_download(self):
        (self.root / 'models' / 'checkpoint_bag_100.pth.tar').unlink()
        with patch.object(api.importlib.util, 'find_spec', return_value=SimpleNamespace(submodule_search_locations=[str(self.root)])), \
                patch.object(api.importlib, 'import_module') as load:
            with self.assertRaises(FileNotFoundError):
                api._load_classifier('core')
            load.assert_not_called()

    def test_empty_unselected_family_also_prevents_import_download(self):
        (self.root / 'models_h' / 'already-installed').unlink()
        with patch.object(api.importlib.util, 'find_spec', return_value=SimpleNamespace(submodule_search_locations=[str(self.root)])), \
                patch.object(api.importlib, 'import_module') as load:
            with self.assertRaises(FileNotFoundError):
                api._load_classifier('core')
            load.assert_not_called()

    def test_complete_selected_ensemble_allows_lazy_import(self):
        with patch.object(api.importlib.util, 'find_spec', return_value=SimpleNamespace(submodule_search_locations=[str(self.root)])), \
                patch.object(api.importlib, 'import_module') as load:
            api._load_classifier('core')
            load.assert_called_once_with('MOFClassifier.CLscore')


if __name__ == '__main__':
    unittest.main()
