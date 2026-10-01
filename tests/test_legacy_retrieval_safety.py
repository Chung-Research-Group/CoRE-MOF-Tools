"""Offline retrieval/transaction tests, with no CSD data or network access."""

import json
from pathlib import Path
import stat
import sys
import tempfile
import types
import unittest
from unittest.mock import Mock, patch
import warnings
import zipfile

from CoREMOF import _legacy_retrieval as retrieval
from CoREMOF import _transactions as transactions

try:
    from CoREMOF import structure
except ImportError:
    structure = None


class RetrievalFilesystemTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.output = self.root / 'output'

    def archive(self, name, members):
        path = self.root / name
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            with zipfile.ZipFile(path, 'w') as bundle:
                for member, content in members:
                    bundle.writestr(member, content)
        return path

    def test_two_archives_same_basename_preserve_exact_bytes(self):
        cr = self.archive('cr.zip', [('CR/a.cif', b'data_cr\n')])
        ncr = self.archive('ncr.zip', [('NCR/a.cif', b'data_ncr\n')])
        retrieval.extract_archives([(cr, None), (ncr, None)], self.output)
        self.assertEqual((self.output / 'CR/a.cif').read_bytes(), b'data_cr\n')
        self.assertEqual((self.output / 'NCR/a.cif').read_bytes(), b'data_ncr\n')
        self.assertEqual(list(self.root.glob('.coremof-download-*')), [])

    def test_rejects_unsafe_members_before_any_publication(self):
        for member in ('../escape.cif', '/absolute.cif', 'CR/../a.cif', 'C:/a.cif',
                       'CR\\a.cif', 'CR//a.cif', './a.cif'):
            with self.subTest(member=member):
                path = self.archive('bad.zip', [('CR/good.cif', b'ok'), (member, b'bad')])
                with self.assertRaises(ValueError):
                    retrieval.extract_archives([(path, None)], self.output)
                self.assertFalse(self.output.exists())

    def test_nul_truncation_is_rejected(self):
        path = self.archive('nul.zip', [('a~truncated.cif', b'bad')])
        path.write_bytes(path.read_bytes().replace(b'a~truncated.cif', b'a\x00truncated.cif'))
        with self.assertRaises(ValueError):
            retrieval.extract_archives([(path, None)], self.output)
        self.assertFalse(self.output.exists())

    def test_rejects_symlink_and_special_members(self):
        for kind in (stat.S_IFLNK, stat.S_IFIFO, stat.S_IFCHR):
            entry = zipfile.ZipInfo('bad.cif')
            entry.create_system = 3
            entry.external_attr = (kind | 0o600) << 16
            path = self.archive('bad.zip', [(entry, b'outside')])
            with self.assertRaises(ValueError):
                retrieval.extract_archives([(path, None)], self.output)
        self.assertFalse(self.output.exists())

    def test_rejects_duplicate_and_file_directory_collision(self):
        for members in ([('a.cif', b'a'), ('a.cif', b'b')],
                        [('a', b'a'), ('a/b.cif', b'b')],
                        [('a/', b''), ('a', b'b')]):
            path = self.archive('bad.zip', members)
            with self.assertRaises(ValueError):
                retrieval.extract_archives([(path, None)], self.output)
            self.assertFalse(self.output.exists())

    def test_cross_archive_conflicts_are_preflighted(self):
        first = self.archive('one.zip', [('CR/a.cif', b'a')])
        second = self.archive('two.zip', [('CR/a.cif', b'b')])
        with self.assertRaises(ValueError):
            retrieval.extract_archives([(first, None), (second, None)], self.output)
        self.assertFalse(self.output.exists())
        second = self.archive('two.zip', [('CR/a.cif/nested.cif', b'b')])
        with self.assertRaises(ValueError):
            retrieval.extract_archives([(first, None), (second, None)], self.output)

    def test_bad_second_archive_leaves_no_partial_first_archive(self):
        first = self.archive('one.zip', [('CR/good.cif', b'ok')])
        second = self.archive('two.zip', [('../bad.cif', b'bad')])
        with self.assertRaises(ValueError):
            retrieval.extract_archives([(first, None), (second, None)], self.output)
        self.assertFalse(self.output.exists())

    def test_crc_failure_leaves_no_outputs(self):
        path = self.archive('bad.zip', [('CR/a.cif', b'good'), ('CR/b.cif', b'unique-content')])
        path.write_bytes(path.read_bytes().replace(b'unique-content', b'broken-content'))
        with self.assertRaises(zipfile.BadZipFile):
            retrieval.extract_archives([(path, None)], self.output)
        self.assertFalse(self.output.exists())
        self.assertEqual(list(self.root.glob('.coremof-download-*')), [])

    def test_existing_files_require_explicit_overwrite_preserve_unrelated(self):
        self.output.mkdir()
        target = self.output / 'a.cif'
        target.write_bytes(b'old')
        other = self.output / 'notes.txt'
        other.write_bytes(b'keep')
        path = self.archive('source.zip', [('a.cif', b'new')])
        with self.assertRaises(FileExistsError):
            retrieval.extract_archives([(path, None)], self.output)
        self.assertEqual(target.read_bytes(), b'old')
        retrieval.extract_archives([(path, None)], self.output, overwrite=True)
        self.assertEqual(target.read_bytes(), b'new')
        self.assertEqual(other.read_bytes(), b'keep')

    def test_existing_symlink_parent_or_file_rejected(self):
        outside = self.root / 'outside'
        outside.mkdir()
        path = self.archive('source.zip', [('CR/a.cif', b'new')])
        self.output.symlink_to(outside, target_is_directory=True)
        with self.assertRaises(ValueError):
            retrieval.extract_archives([(path, None)], self.output)
        self.output.unlink()
        (self.output / 'CR').mkdir(parents=True)
        (self.output / 'CR/a.cif').symlink_to(outside / 'missing.cif')
        with self.assertRaises(ValueError):
            retrieval.extract_archives([(path, None)], self.output, overwrite=True)
        self.assertEqual(list(outside.iterdir()), [])

    def test_missing_entry_is_no_op(self):
        path = self.archive('source.zip', [('CR/a.cif', b'a')])
        retrieval.extract_archives([(path, 'missing.cif')], self.output)
        self.assertFalse(self.output.exists())

    def test_rollback_restores_two_equal_basenames(self):
        for part in ('CR', 'NCR'):
            (self.output / part).mkdir(parents=True)
            (self.output / part / 'a.cif').write_text('old-' + part)
        path = self.archive('source.zip', [('CR/a.cif', b'new-cr'),
                                         ('NCR/a.cif', b'new-ncr'), ('third.cif', b'new')])
        real_replace = transactions.os.replace

        def fail_third(source, target):
            if Path(target) == self.output / 'third.cif':
                raise OSError('simulated publication failure')
            return real_replace(source, target)

        with patch.object(transactions.os, 'replace', side_effect=fail_third):
            with self.assertRaisesRegex(OSError, 'simulated'):
                retrieval.extract_archives([(path, None)], self.output, overwrite=True)
        for part in ('CR', 'NCR'):
            self.assertEqual((self.output / part / 'a.cif').read_text(), 'old-' + part)
        self.assertFalse((self.output / 'third.cif').exists())

    def test_create_if_absent_race_preserves_other_writer_and_rolls_back(self):
        path = self.archive('source.zip', [('CR/a.cif', b'a'), ('CR/b.cif', b'b')])
        real_link = transactions.os.link

        def competing_writer(source, target):
            if Path(target) == self.output / 'CR/b.cif':
                Path(target).write_bytes(b'concurrent')
            return real_link(source, target)

        with patch.object(transactions.os, 'link', side_effect=competing_writer):
            with self.assertRaises(FileExistsError):
                retrieval.extract_archives([(path, None)], self.output)
        self.assertEqual((self.output / 'CR/b.cif').read_bytes(), b'concurrent')
        self.assertFalse((self.output / 'CR/a.cif').exists())

    def test_failed_rollback_preserves_backup_and_reports_staging(self):
        self.output.mkdir()
        (self.output / 'a.cif').write_bytes(b'old')
        path = self.archive('source.zip', [('a.cif', b'new'), ('b.cif', b'b')])
        real_replace = transactions.os.replace

        def fail_publish_and_restore(source, target):
            if Path(target).name == 'b.cif' or str(target).endswith('.rollback'):
                raise OSError('simulated failure')
            return real_replace(source, target)

        with patch.object(transactions.os, 'replace', side_effect=fail_publish_and_restore):
            with self.assertRaises(OSError) as caught:
                retrieval.extract_archives([(path, None)], self.output, overwrite=True)
        staging = Path(caught.exception.coremof_preserved_staging_directory)
        self.assertTrue(staging.is_dir())
        self.assertEqual([p.read_bytes() for p in staging.glob('*.previous')], [b'old'])

    def test_safe_refcode_not_a_claim_of_full_identifier_grammar(self):
        self.assertEqual(retrieval.validate_refcode('ABCDEF01'), 'ABCDEF01')
        for refcode in ('../ABCDEF', '/tmp/x', '', None, 'a/b', 'a\\b', 'a\n'):
            with self.assertRaises(ValueError):
                retrieval.validate_refcode(refcode)

    def test_metadata_validation_rejects_non_object(self):
        path = self.root / 'download'
        for value in ([], None, 1, 'error'):
            path.write_text(json.dumps(value))
            with self.assertRaises(ValueError):
                retrieval.validate_cache(path, 'data/CR.json')

    def test_overwrite_flag_must_be_boolean(self):
        with self.assertRaises(TypeError):
            retrieval.extract_archives([], self.output, overwrite='False')
        with self.assertRaises(TypeError):
            retrieval.publish_text(self.output / 'a.cif', 'a', overwrite=1)
        self.assertFalse(self.output.exists())


@unittest.skipIf(structure is None, 'legacy API requires requests and gemmi')
class RetrievalPublicApiTests(unittest.TestCase):
    setUp = RetrievalFilesystemTests.setUp
    archive = RetrievalFilesystemTests.archive

    def response(self, content):
        response = Mock()
        response.__enter__ = Mock(return_value=response)
        response.__exit__ = Mock(return_value=False)
        response.iter_content.return_value = [content]
        return response

    def test_http_success_with_malformed_payload_is_not_cached(self):
        response = self.response(b'<html>error</html>')
        with patch.object(structure, 'package_directory', str(self.root)), patch.object(structure.requests, 'get', return_value=response):
            with self.assertRaises(ValueError):
                structure._ensure_data_file('data/CR.json')
        self.assertEqual(list((self.root / 'data').iterdir()), [])

    def test_cache_created_during_download_is_not_replaced(self):
        response = self.response(b'{"new":1}')

        def chunks(**kwargs):
            (self.root / 'data/CR.json').write_bytes(b'{"winner":2}')
            yield b'{"new":1}'

        response.iter_content.side_effect = chunks
        with patch.object(structure, 'package_directory', str(self.root)), patch.object(structure.requests, 'get', return_value=response):
            path = structure._ensure_data_file('data/CR.json')
        self.assertEqual(path.read_bytes(), b'{"winner":2}')
        self.assertEqual(list(path.parent.iterdir()), [path])

    def test_si_path_object_and_explicit_overwrite(self):
        cr = self.archive('cr.zip', [('CR/a.cif', b'cr')])
        ncr = self.archive('ncr.zip', [('NCR/a.cif', b'ncr')])
        with patch.object(structure, '_ensure_data_file', side_effect=lambda name: cr if name.endswith('/CR.zip') else ncr):
            instance = structure.download_from_SI(self.output)
            self.assertEqual(instance.list_zip(cr), ['CR/a.cif'])
            with self.assertRaises(FileExistsError):
                structure.download_from_SI(self.output)
            structure.download_from_SI(self.output, overwrite=True)
            self.assertIsNone(instance.get_from_SI(cr, 'missing', self.output))

    def fake_ccdc(self, *, failure=None):
        reader = Mock()
        reader.crystal.return_value.to_string.return_value = 'data_exact\n# CIF text\n'
        if failure:
            reader.crystal.side_effect = failure
        ccdc = types.ModuleType('ccdc')
        ccdc.io = types.SimpleNamespace(EntryReader=Mock(return_value=reader))
        return ccdc, reader

    def test_csd_closes_reader_and_preserves_exact_output(self):
        ccdc, reader = self.fake_ccdc()
        with patch.dict(sys.modules, {'ccdc': ccdc}):
            self.assertIsNone(structure.download_from_CSD('ABCDEF', self.output))
        reader.close.assert_called_once_with()
        reader.crystal.assert_called_once_with('ABCDEF')
        self.assertEqual((self.output / 'ABCDEF.cif').read_bytes(), b'data_exact\n# CIF text\n')

    def test_csd_failure_closes_reader_without_output(self):
        ccdc, reader = self.fake_ccdc(failure=RuntimeError('licensed reader failed'))
        with patch.dict(sys.modules, {'ccdc': ccdc}):
            with self.assertRaisesRegex(RuntimeError, 'licensed reader failed'):
                structure.download_from_CSD('ABCDEF', self.output)
        reader.close.assert_called_once_with()
        self.assertFalse(self.output.exists())

    def test_csd_rejects_invalid_refcode_or_existing_output_before_reader(self):
        ccdc, reader = self.fake_ccdc()
        self.output.mkdir()
        target = self.output / 'ABCDEF.cif'
        target.write_bytes(b'keep')
        with patch.dict(sys.modules, {'ccdc': ccdc}):
            with self.assertRaises(ValueError):
                structure.download_from_CSD('../ABCDEF', self.output)
            with self.assertRaises(FileExistsError):
                structure.download_from_CSD('ABCDEF', self.output)
            ccdc.io.EntryReader.assert_not_called()
            structure.download_from_CSD('ABCDEF', self.output, overwrite=True)
        self.assertEqual(target.read_bytes(), b'data_exact\n# CIF text\n')


if __name__ == '__main__':
    unittest.main()
