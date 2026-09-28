"""Offline adapter regressions, without scientific dependencies or network."""
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from CoREMOF import cod_curation as api
from CoREMOF import _cod_curation_runtime as runtime


class CODCurationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.artifact = self.root / 'artifact.cif'
        self.artifact.write_text('data_test\n')
        self.spec = {'bytes': self.artifact.stat().st_size, 'sha256': api._sha(self.artifact)}
        self.record = {'pipeline_id': api.METHOD, 'source_id': 'cod_1234567',
                       'source': {'sha256': '1' * 64}, 'status': 'complete', 'review_reasons': [],
                       'children': [{'status': 'complete', 'category': 'FSR_ASR_IDENTICAL',
                                     'review_reasons': [], 'invariants': {'all_passed': True},
                                     'deliverables': [{'path': str(self.artifact), 'identity': self.spec,
                                                       'atom_count': 1, 'charged': True,
                                                       'pacman': {'charge_count': 1, 'charge_sum': 0.,
                                                                  'used_for_structural_selection': False}}]}]}

    def validate(self, record=None, skip=False, root=None):
        api._validate_record(record or self.record, 'cod_1234567', '1' * 64, root or self.root, skip)

    def test_valid_complete_record(self):
        self.validate()

    def test_wrong_source_method_or_hash(self):
        for key, value in [('source_id', 'other'), ('pipeline_id', 'other'), ('source', {'sha256': '2'*64})]:
            record = copy.deepcopy(self.record)
            record[key] = value
            with self.subTest(key=key), self.assertRaises(api.CODCurationError):
                self.validate(record)

    def test_uncharged_never_complete(self):
        with self.assertRaises(api.CODCurationError):
            self.validate(skip=True)

    def test_failed_invariants_never_complete(self):
        self.record['children'][0]['invariants']['all_passed'] = False
        with self.assertRaises(api.CODCurationError):
            self.validate()

    def test_counterion_must_skip_asr(self):
        child = self.record['children'][0]
        child.update(category='ION_FSR', asr={'executed': True})
        with self.assertRaises(api.CODCurationError):
            self.validate()
        child['asr']['executed'] = False
        self.validate()

    def test_review_reason_required(self):
        self.record['status'] = 'review'
        self.record['children'][0]['status'] = 'review'
        with self.assertRaises(api.CODCurationError):
            self.validate()
        self.record['children'][0]['review_reasons'] = ['partial_occupancy']
        self.validate()

    def test_nonfinite_charge_or_wrong_count_rejected(self):
        for value in (float('nan'), float('inf'), .001, True, '0'):
            record = copy.deepcopy(self.record)
            record['children'][0]['deliverables'][0]['pacman']['charge_sum'] = value
            with self.subTest(value=value), self.assertRaises(api.CODCurationError):
                self.validate(record)
        self.record['children'][0]['deliverables'][0]['pacman']['charge_count'] = 2
        with self.assertRaises(api.CODCurationError):
            self.validate()

    def test_artifact_escape_rejected(self):
        with self.assertRaises(api.CODCurationError):
            self.validate(root=self.root / 'different')

    def test_changed_artifact_rejected(self):
        self.artifact.write_text('data_changed\n')
        with self.assertRaises(api.CODCurationError):
            self.validate()

    def test_stage_exact_private_copies_and_no_overwrite(self):
        profile = {'assets': {'artifact.cif': dict(self.spec, required_for='all')}}
        dest = self.root / 'private'
        result = api._stage_assets(self.root, dest, profile, False)
        self.assertEqual(result, profile['assets'])
        (dest / 'artifact.cif').write_text('different')
        self.assertEqual(api._sha(self.artifact), self.spec['sha256'])
        with self.assertRaises(FileExistsError):
            api._stage_assets(self.root, dest, profile, False)

    def test_missing_asset_never_downloads(self):
        profile = {'assets': {'missing.pth': dict(self.spec, required_for='charging')}}
        with self.assertRaises(api.CODCurationError):
            api._stage_assets(self.root, self.root / 'private', profile, False)
        self.assertEqual(api._stage_assets(self.root, self.root / 'private', profile, True), {})

    def test_profile_paths_cannot_escape(self):
        for name in ('../bad', '/bad'):
            with self.subTest(name=name), self.assertRaises(api.CODCurationError):
                api._stage_assets(self.root, self.root / 'private',
                                  {'assets': {name: dict(self.spec, required_for='all')}}, False)

    def test_profile_has_all_eager_pacman_models(self):
        profile = json.loads(api.PROFILE_PATH.read_text())
        self.assertEqual(profile['profile'], api.METHOD)
        for model in ('cm5', 'bader', 'ddec', 'repeat', 'pbe', 'bandgap'):
            for ext in ('pkl', 'pth'):
                name = f'v4_1/runtime/PACMANCharge/{model}.{ext}'
                self.assertIn(name, profile['assets'])
        self.assertFalse(any('/home/' in name for name in profile['assets']))
        self.assertEqual(profile['scientific_policy']['total_pair_margin_angstrom'], .25)
        self.assertEqual(profile['scientific_policy']['ase_skin_angstrom'], .125)

    def test_parameters_rejected_before_runtime(self):
        common = dict(output_dir=self.root / 'out', python='missing', workflow_root='missing')
        for value in ('../1', '', 'COD123', True, 123):
            with self.subTest(value=value), self.assertRaises(ValueError):
                api.curate_cod_cif_v41(self.artifact, cod_id=value, **common)
        for value in (0, -1, True, 1.5):
            with self.subTest(timeout=value), self.assertRaises(ValueError):
                api.curate_cod_cif_v41(self.artifact, cod_id='1234567', timeout_seconds=value, **common)

    def test_destination_never_overwritten(self):
        with self.assertRaises(FileExistsError):
            api.curate_cod_cif_v41(self.artifact, cod_id='1234567', output_dir=self.artifact,
                                  python='missing', workflow_root='missing')
        self.assertEqual(api._sha(self.artifact), self.spec['sha256'])

    def test_runtime_mismatch_rejected(self):
        profile = json.loads(api.PROFILE_PATH.read_text())
        with patch.object(runtime.sys, 'version_info', (3, 9)), patch.object(runtime.importlib.metadata, 'version', return_value='0'):
            with self.assertRaisesRegex(RuntimeError, 'dependency mismatch'):
                runtime.inspect_runtime(profile)

    def test_offline_guard_denies_connections(self):
        with self.assertRaisesRegex(RuntimeError, 'disabled'):
            runtime.deny_network('host', 443)

    def test_canonical_hash_row_order_independent(self):
        self.assertEqual(api._canonical_sha({'a': 1, 'b': 2}), api._canonical_sha({'b': 2, 'a': 1}))
        self.assertNotEqual(api._canonical_sha({'a': 1}), hashlib.sha256(b'other').hexdigest())


if __name__ == '__main__':
    unittest.main()
