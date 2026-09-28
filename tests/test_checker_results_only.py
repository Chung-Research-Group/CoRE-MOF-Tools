"""Results remain usable without bundling or executing external checkers."""
import ast
import contextlib
import importlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from CoREMOF._checker_execution import CheckerExecutionUnavailableError
from CoREMOF.dataset import CoREMOFDataset
from test_dataset_labels import _make_release

ROOT = Path(__file__).resolve().parents[1]
REMOVED = (
    '_release_checkers_protocol.py', '_release_checkers_worker.py',
    '_release_mosaec_worker.py', '_release_setc_protocol.py', '_release_setc_worker.py',
)


class CheckerResultsOnlyTests(unittest.TestCase):
    def test_no_external_workers_or_reference_tables_in_source(self):
        for name in REMOVED:
            self.assertFalse((ROOT / 'CoREMOF' / name).exists(), name)
        self.assertFalse(list((ROOT / 'CoREMOF/data/mosaec').glob('*')))
        tree = ast.parse((ROOT / 'CoREMOF/mosaec.py').read_text())
        self.assertEqual({n.name for n in tree.body if isinstance(n, ast.FunctionDef)}, {'run', 'check'})

    def test_retired_interfaces_fail_before_io_or_process_creation(self):
        functions = [('release_checkers', 'calculate_release_checkers'),
                     ('release_mosaec', 'calculate_release_mosaec'),
                     ('release_setc', 'calculate_release_setc'), ('mosaec', 'run'), ('mosaec', 'check')]
        for module, name in functions:
            fn = getattr(importlib.import_module('CoREMOF.' + module), name)
            with self.subTest(module=module, function=name), tempfile.TemporaryDirectory() as tmp:
                target = Path(tmp) / 'must_not_exist'
                with patch('subprocess.Popen') as popen, patch('builtins.open') as opened:
                    with self.assertRaisesRegex(CheckerExecutionUnavailableError, 'precomputed|existing release'):
                        fn('/nonexistent.cif', output_dir=target)
                    popen.assert_not_called()
                    opened.assert_not_called()
                self.assertFalse(target.exists())

    def test_legacy_curate_methods_are_also_results_only(self):
        nodes = [n for n in ast.parse((ROOT / 'CoREMOF/curate.py').read_text()).body
                 if getattr(n, 'name', '') in {'mof_check', 'run_MOSAEC'}]
        namespace = {'__package__': 'CoREMOF'}
        exec(compile(ast.Module(body=nodes, type_ignores=[]), 'curate migration notices', 'exec'), namespace)
        functions = [namespace['mof_check'], namespace['run_MOSAEC']]
        for name in ('check', 'Chen_Manz', 'mof_checker'):
            functions.append(getattr(namespace['mof_check'], name))
        for fn in functions:
            with self.subTest(fn=fn.__name__), self.assertRaises(CheckerExecutionUnavailableError):
                fn('/nonexistent.cif')

    def test_minimal_import_does_not_load_external_checkers(self):
        code = '''
import sys
from CoREMOF import release_checkers, release_mosaec, release_setc, mosaec
from CoREMOF.dataset import CoREMOFDataset
for name in ('ccdc', 'mofchecker', 'torch', 'numpy', 'CoREMOF._release_setc_protocol'):
    assert name not in sys.modules, name
'''
        subprocess.run([sys.executable, '-B', '-S', '-c', code], cwd=ROOT, check=True, capture_output=True)

    def test_cli_help_works_but_execution_is_rejected(self):
        for module in ('release_checkers', 'release_mosaec', 'release_setc'):
            command = [sys.executable, '-B', '-S', '-m', 'CoREMOF.' + module]
            help_run = subprocess.run(command + ['--help'], cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(help_run.returncode, 0, help_run.stderr)
            run = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(run.returncode, 2)
            self.assertIn('Results-only', run.stderr)

    def test_all_saved_checker_states_and_custom_consensus_remain_usable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _make_release(root)
            dataset = CoREMOFDataset.from_release(root)
            before = (root / 'metadata/metadata.csv').read_bytes()
            with patch('subprocess.Popen') as process:
                view = dataset.classify('5checker')
                self.assertEqual(dict(view.label_counts()), {'CR': 1, 'NCR': 1, 'AMBIGUOUS': 1, 'UNCHECKED': 1})
                custom = dataset.classify(checkers=('MOFChecker', 'MOSAEC'))
                self.assertEqual(len(custom), 4)
                process.assert_not_called()
            self.assertEqual((root / 'metadata/metadata.csv').read_bytes(), before)

    def test_read_only_example_executes_on_release_fixture(self):
        namespace = {'__name__': 'test_example'}
        exec(compile((ROOT / 'examples/read_checker_results.py').read_text(), 'read_checker_results.py', 'exec'), namespace)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _make_release(root)
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                self.assertEqual(namespace['main']([str(root), '--structure-id', 'ASR-COD-2026-0001']), 0)
            result = json.loads(output.getvalue())
            self.assertEqual(result['structure']['label'], 'CR')
            self.assertEqual(len(result['structure']['checkers']), 5)


if __name__ == '__main__':
    unittest.main()
