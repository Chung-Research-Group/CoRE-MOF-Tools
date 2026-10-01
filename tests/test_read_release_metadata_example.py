"""Read-only access example integration with disposable release fixtures."""
import contextlib
import hashlib
import importlib.util
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from test_dataset_labels import _make_release

_PATH = Path(__file__).resolve().parents[1] / 'examples/read_release_metadata.py'
_SPEC = importlib.util.spec_from_file_location('metadata_access_example', _PATH)
example = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(example)


class MetadataAccessExampleTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / 'release'
        _make_release(self.root)

    def invoke(self, *options):
        stream = io.StringIO()
        with patch.object(sys, 'argv', [str(_PATH), str(self.root), *options]), \
                contextlib.redirect_stdout(stream):
            example.main()
        return json.loads(stream.getvalue())

    def test_metadata_only_version_and_no_write(self):
        before = {p.relative_to(self.root): p.read_bytes()
                  for p in self.root.rglob('*') if p.is_file()}
        result = self.invoke('--expected-version', 'vtest')
        self.assertEqual(result['structures'], 4)
        self.assertEqual(result['dataset_version'], 'vtest')
        self.assertFalse(result['cif_bytes_verified'])
        self.assertFalse(result['new_calculations'])
        self.assertFalse(result['publication_authorization_granted_by_example'])
        after = {p.relative_to(self.root): p.read_bytes()
                 for p in self.root.rglob('*') if p.is_file()}
        self.assertEqual(before, after)

    def test_wrong_version_refused(self):
        with self.assertRaisesRegex(ValueError, 'explicitly requested'):
            self.invoke('--expected-version', 'vwrong')

    def test_projection_requires_hash(self):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            self.invoke('--expected-version', 'vtest', '--projection-contract', 'contract.json')

    def test_trusted_ledger_preserves_zeros_and_nulls(self):
        payload = self.root / 'targets.json'
        payload.write_text('{"finite_zero": 0, "missing": null}\n')
        original = payload.read_bytes()
        self.write_ledger('targets.json', original)
        expected = hashlib.sha256((self.root / 'SHA256SUMS').read_bytes()).hexdigest()
        self.assertEqual(example.verify_metadata_ledger(self.root, expected), 1)
        self.assertEqual(original, payload.read_bytes())

    def write_ledger(self, name, content):
        (self.root / 'SHA256SUMS').write_text(hashlib.sha256(content).hexdigest() + '  ' + name + '\n')

    def test_tampered_ledger_or_payload_refused(self):
        self.write_ledger('dataset_info.json', (self.root / 'dataset_info.json').read_bytes())
        pin = hashlib.sha256((self.root / 'SHA256SUMS').read_bytes()).hexdigest()
        with self.assertRaisesRegex(ValueError, 'independently received'):
            example.verify_metadata_ledger(self.root, '0' * 64)
        (self.root / 'dataset_info.json').write_text('{}')
        with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
            example.verify_metadata_ledger(self.root, pin)

    def test_unsafe_ledger_path_refused(self):
        for name in ('../outside', '/outside', 'a/../b', 'a\\b', './a', 'a//b'):
            with self.subTest(name=name):
                self.write_ledger(name, b'fixture')
                pin = hashlib.sha256((self.root / 'SHA256SUMS').read_bytes()).hexdigest()
                with self.assertRaisesRegex(ValueError, 'Invalid metadata ledger'):
                    example.verify_metadata_ledger(self.root, pin)

    def test_parent_symlink_refused(self):
        with tempfile.TemporaryDirectory() as outside:
            file = Path(outside) / 'payload'
            file.write_bytes(b'fixture')
            (self.root / 'linked').symlink_to(outside, target_is_directory=True)
            self.write_ledger('linked/payload', file.read_bytes())
            pin = hashlib.sha256((self.root / 'SHA256SUMS').read_bytes()).hexdigest()
            with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
                example.verify_metadata_ledger(self.root, pin)


if __name__ == '__main__':
    unittest.main()
