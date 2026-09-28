"""Exercise public PACMAN wrapper I/O with a stub predictor, not new science."""
from contextlib import redirect_stdout
import io
from pathlib import Path
import tempfile
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from CoREMOF.prediction import pacman


class PacmanPredictionIsolationTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.source = self.root / 'source.cif'
        self.source.write_bytes(b'original source CIF fixture')
        self.output = self.root / 'results'
        self.received = []

        def predict(**kwargs):
            path = Path(kwargs['cif_file'])
            self.received.append(kwargs)
            self.assertNotEqual(path.parent, self.source.parent)
            self.assertEqual(path.read_bytes(), self.source.read_bytes())
            # Some predictor versions modify the file given to them.
            path.write_bytes(b'modified isolated input')
            path.with_name(path.stem + '_pacman.cif').write_bytes(b'charged CIF fixture')

        self.backend = SimpleNamespace(predict=Mock(side_effect=predict), Energy=Mock(return_value=(-10.0, 2.0)))
        backend_patch = patch.dict(sys.modules, {'PACMANCharge': SimpleNamespace(pmcharge=self.backend)})
        backend_patch.start()
        self.addCleanup(backend_patch.stop)
        self.predict = pacman

    def test_source_is_unchanged_and_backend_options_preserved(self):
        result = self.predict(self.source, self.output, digits=7, neutral=False, keep_connect=True)
        self.assertEqual(result, {'Name': 'source', 'PBE Energy': -10.0, 'Bandgap': 2.0})
        self.assertEqual(self.source.read_bytes(), b'original source CIF fixture')
        self.assertEqual((self.output / 'source_pacman.cif').read_bytes(), b'charged CIF fixture')
        self.assertFalse((self.root / 'source_pacman.cif').exists())
        self.assertFalse(Path(self.received[0]['cif_file']).exists())
        self.assertEqual({key: value for key, value in self.received[0].items() if key != 'cif_file'},
            {'charge_type': 'DDEC6', 'digits': 7, 'atom_type': True, 'neutral': False, 'keep_connect': True})
        self.assertEqual(self.backend.Energy.call_args.kwargs['cif_file'], self.received[0]['cif_file'])

    def test_existing_output_is_preserved_without_starting_predictor(self):
        self.output.mkdir()
        destination = self.output / 'source_pacman.cif'
        destination.write_bytes(b'prior result')
        with self.assertRaises(FileExistsError):
            self.predict(self.source, self.output)
        self.backend.predict.assert_not_called()
        self.assertEqual(destination.read_bytes(), b'prior result')

    def test_backend_failure_does_not_publish_partial_result(self):
        self.backend.Energy.side_effect = RuntimeError('energy model failed')
        with redirect_stdout(io.StringIO()):
            self.assertIsNone(self.predict(self.source, self.output))
        self.assertFalse(self.output.exists())
        self.assertEqual(self.source.read_bytes(), b'original source CIF fixture')
        self.assertFalse(Path(self.received[0]['cif_file']).exists())

    def test_missing_or_empty_output_is_not_published(self):
        for contents in (None, b''):
            with self.subTest(contents=contents):
                def predict(cif_file, **kwargs):
                    if contents is not None:
                        path = Path(cif_file)
                        path.with_name(path.stem + '_pacman.cif').write_bytes(contents)
                self.backend.predict.side_effect = predict
                with redirect_stdout(io.StringIO()):
                    self.assertIsNone(self.predict(self.source, self.output))
                self.assertFalse(self.output.exists())

    def test_symlinked_output_is_rejected_before_backend_runs(self):
        actual = self.root / 'actual'
        actual.mkdir()
        self.output.symlink_to(actual, target_is_directory=True)
        with self.assertRaises(ValueError):
            self.predict(self.source, self.output)
        self.backend.predict.assert_not_called()

    def test_uppercase_cif_suffix_has_a_normalized_isolated_input(self):
        renamed = self.root / 'source.CIF'
        self.source.rename(renamed)
        self.source = renamed
        self.assertIsNotNone(self.predict(renamed, self.output))
        self.assertTrue(self.received[0]['cif_file'].endswith('source.cif'))


if __name__ == '__main__':
    unittest.main()
