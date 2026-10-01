"""A failed legacy download must not leave an anonymous partial cache file."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch

try:
    from CoREMOF import structure
except ImportError:
    structure = None


@unittest.skipIf(structure is None, 'legacy structure API requires requests and gemmi')
class LegacyDownloadCleanupTests(unittest.TestCase):
    def test_stream_failure_removes_partial_temporary_file(self):
        response = Mock()
        response.__enter__ = Mock(return_value=response)
        response.__exit__ = Mock(return_value=False)

        def chunks(**kwargs):
            yield b'partial'
            raise OSError('connection interrupted')

        response.iter_content.side_effect = chunks
        with tempfile.TemporaryDirectory() as root:
            with patch.object(structure, 'package_directory', root), patch.object(structure.requests, 'get', return_value=response):
                with self.assertRaisesRegex(OSError, 'interrupted'):
                    structure._ensure_data_file('data/CR.json')
            self.assertEqual(list((Path(root) / 'data').iterdir()), [])

    def test_success_publishes_complete_file(self):
        response = Mock()
        response.__enter__ = Mock(return_value=response)
        response.__exit__ = Mock(return_value=False)
        response.iter_content.return_value = [b'{', b'"complete":true', b'}']
        with tempfile.TemporaryDirectory() as root:
            with patch.object(structure, 'package_directory', root), patch.object(structure.requests, 'get', return_value=response):
                path = structure._ensure_data_file('data/CR.json')
            self.assertEqual(path.read_bytes(), b'{"complete":true}')
            self.assertEqual(list(path.parent.iterdir()), [path])


if __name__ == '__main__':
    unittest.main()
