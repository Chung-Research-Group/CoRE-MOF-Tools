"""Pure-Python OMS parsing, isolation and failure-contract tests."""
import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from CoREMOF import release_oms as api


class ReleaseOMSTests(unittest.TestCase):
    SID = '2013[Cu][nan]3[ASR]5'

    def data(self, count=4):
        return {'structure_id': self.SID, 'protocol_id': api.protocol.PROTOCOL_ID,
                'input': {'cif_sha256': 'a' * 64}, 'execution_status': 'SUCCESS',
                'open_metal_site_props': {'open_metal_site_count': count, 'has_open_metal_sites': count > 0,
                                          'surface_definition_A': 3.5, 'raw_line': f'file.oms #OMS= {count}'}}

    def test_pinned_assets(self):
        self.assertEqual(api._sha(api.protocol.__file__), api.PROTOCOL_SHA256)
        self.assertEqual(api._sha(Path(api.__file__).with_name('data') / 'zeopp_oms_metadata_contract_v1.json'), api.CONTRACT_SHA256)

    def test_zero_is_success(self):
        api._validate(self.data(0), self.SID, 'a' * 64)

    def test_positive_count(self):
        api._validate(self.data(), self.SID, 'a' * 64)

    def test_malformed_or_multiple_lines(self):
        for text in ('file.oms #OMS= -1', 'file.oms #OMS= 2.5', 'file.oms #OMS= 2\nfile.oms #OMS= 2'):
            with self.assertRaises(api.protocol.ZeoppOMSError):
                api.protocol.parse_oms(text)

    def test_boolean_and_count_must_agree(self):
        r = self.data()
        r['open_metal_site_props']['has_open_metal_sites'] = False
        with self.assertRaises(api.ReleaseOMSError):
            api._validate(r, self.SID, 'a' * 64)

    def test_no_value_for_failure(self):
        r = self.data()
        r['execution_status'] = 'ERROR'
        with self.assertRaises(api.ReleaseOMSError):
            api._validate(r, self.SID, 'a' * 64)
        r['open_metal_site_props'] = None
        api._validate(r, self.SID, 'a' * 64)

    def test_surface_can_be_missing_but_not_nonfinite(self):
        r = self.data()
        r['open_metal_site_props']['surface_definition_A'] = None
        api._validate(r, self.SID, 'a' * 64)
        r['open_metal_site_props']['surface_definition_A'] = float('nan')
        with self.assertRaises(api.ReleaseOMSError):
            api._validate(r, self.SID, 'a' * 64)

    def test_surface_conflicting_stdout(self):
        self.assertIsNone(api.protocol.parse_surface_definition('no surface diagnostic'))
        with self.assertRaises(api.protocol.ZeoppOMSError):
            api.protocol.parse_surface_definition('Surface definition = 3.5\nSurface definition = 2.5')

    def fake_run(self, fail=False, timeout=False):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cif = root / 'original.cif'
            cif.write_text('data_original\n')
            original = cif.read_bytes()
            network = root / 'network'
            network.write_text('not executed')
            pending = [timeout]

            class Process:
                returncode = 1 if fail else 0
                pid = 99999

                def wait(self, timeout=None):
                    if pending.pop() if pending else False:
                        raise api.subprocess.TimeoutExpired('oms', timeout)

            def launch(command, **kwargs):
                self.assertTrue(kwargs['start_new_session'])
                self.assertIn('-S', command)
                self.assertNotIn('PYTHONPATH', kwargs['env'])
                flags = dict(zip(command[4::2], command[5::2]))
                with Path(flags['--manifest']).open() as stream:
                    row = next(csv.DictReader(stream))
                self.assertNotEqual(Path(row['canonical_cif_path']), cif)
                self.assertEqual(Path(row['canonical_cif_path']).read_bytes(), original)
                r = self.data()
                r['input']['cif_sha256'] = hashlib.sha256(original).hexdigest()
                folder = Path(flags['--output-root']) / 'records'
                folder.mkdir()
                (folder / (self.SID + '.json')).write_text(json.dumps(r))
                return Process()

            real_check = api._check

            def check(path, digest):
                return Path(path) if Path(path) == network else real_check(path, digest)

            with patch.object(api, '_check', side_effect=check), patch.object(api.subprocess, 'Popen', side_effect=launch), patch.object(api, '_stop') as stop:
                result = api.calculate_release_oms(cif, structure_id=self.SID, source_database='COD', output_dir=root / 'result', network=network)
            self.assertEqual(cif.read_bytes(), original)
            self.assertEqual(stop.called, timeout)
            self.assertFalse(any(root.glob('.zeopp-oms-replay-*')))
            self.assertTrue((root / 'result' / 'receipt.json').exists())
            self.assertFalse(result['automatic_cr_ncr_exclusion_authorized'])
            return result

    def test_isolated_transaction(self):
        self.assertEqual(self.fake_run()['open_metal_site_props']['open_metal_site_count'], 4)

    def test_process_failure_is_unavailable(self):
        self.assertIsNone(self.fake_run(fail=True)['open_metal_site_props'])

    def test_timeout_stops_only_own_process(self):
        self.assertEqual(self.fake_run(timeout=True)['error_type'], 'TIMEOUT')


if __name__ == '__main__':
    unittest.main()
