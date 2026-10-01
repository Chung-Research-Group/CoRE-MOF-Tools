"""Dependency-minimal guards for the explicit v4.1 curation replay."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from CoREMOF import release_curation as api
from CoREMOF import _release_curation_worker as worker


def record(root, *, status='complete', skip=False, ion=False):
    candidate = root / 'candidate.cif'
    candidate.write_text('data_fixture\n')
    digest = worker.sha(candidate)
    charge = {'status': 'skipped_by_request' if skip else 'not_run_review' if status == 'review' else 'complete'}
    if charge['status'] == 'complete':
        charge.update(charged_cif=str(candidate), charged_sha256=digest, net_charge_check='passed')
    return {'refcode': 'SOURCE1', 'source_sha256': 'a' * 64,
            'pipeline_version': api.protocol.PIPELINE_VERSION, 'status': status,
            'category': 'ION' if ion else 'FSR_ASR', 'release_eligible': False,
            'curation_stage_eligible': status == 'complete' and not skip,
            'review_reasons': ['unknown removed component'] if status == 'review' else [],
            'invariants': {'all_required_invariants_passed': True},
            'curation': {'asr': {'status': 'skipped'} if ion else {}},
            'deliverables': [{'roles': ['ION_FSR'] if ion else ['FSR', 'ASR'],
                              'uncharged_cif': str(candidate), 'uncharged_sha256': digest,
                              'release_eligible': False,
                              'curation_stage_eligible': status == 'complete' and not skip,
                              'pacman': charge}]}


class CurationContracts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def validate(self, data, skip=False):
        return api._validate(data, 'SOURCE1', 'a' * 64, skip, self.root)

    def test_protocol_is_byte_identical(self):
        self.assertEqual(worker.sha(api.protocol.__file__), worker.PROTOCOL_SHA256)

    def test_corrected_margin_and_sequential_branch_are_frozen(self):
        self.assertEqual(api.protocol.INITIAL_TOTAL_PAIR_MARGIN_ANGSTROM, 0.25)
        self.assertEqual(api.protocol.ADAPTIVE_TOTAL_MARGIN_INCREMENT_ANGSTROM, 0.05)
        self.assertEqual(api.protocol.PACMAN_TORCH_SEED, 0)
        self.assertEqual(api.protocol.build_deliverable_specs('ION', 'fsr', None, [0], None, False)[0]['roles'], ['ION_FSR'])
        specs = api.protocol.build_deliverable_specs('FSR_ASR', 'fsr', 'asr', [0], [0], True)
        self.assertEqual(len(specs), 1)
        self.assertEqual(specs[0]['roles'], ['FSR', 'ASR'])

    def test_complete_validated_charge(self):
        self.validate(record(self.root))

    def test_skip_charge_not_ready(self):
        self.validate(record(self.root, skip=True), skip=True)

    def test_review_not_ncr(self):
        data = record(self.root, status='review')
        self.validate(data)
        self.assertNotIn('NCR', json.dumps(data))

    def test_ion_requires_asr_skipped(self):
        data = record(self.root, ion=True)
        self.validate(data)
        data['curation']['asr']['status'] = 'complete'
        with self.assertRaises(api.ReleaseCurationError):
            self.validate(data)

    def test_curation_never_marks_release_ready(self):
        data = record(self.root)
        data['release_eligible'] = True
        with self.assertRaises(api.ReleaseCurationError):
            self.validate(data)

    def test_unavailable_charge_cannot_be_ready(self):
        data = record(self.root)
        data['deliverables'][0]['pacman']['status'] = 'failed'
        with self.assertRaises(api.ReleaseCurationError):
            self.validate(data)

    def test_bad_atom_mapping_is_not_complete(self):
        data = record(self.root)
        data['invariants']['all_required_invariants_passed'] = False
        with self.assertRaises(api.ReleaseCurationError):
            self.validate(data)

    def test_failed_needs_diagnostic(self):
        data = record(self.root)
        data['status'] = 'failed'
        with self.assertRaises(api.ReleaseCurationError):
            self.validate(data)
        data['errors'] = [{'type': 'PARSE_ERROR'}]
        self.validate(data)

    def test_output_hash_change_rejected(self):
        data = record(self.root)
        (self.root / 'candidate.cif').write_text('modified')
        with self.assertRaises(RuntimeError):
            self.validate(data)

    def test_output_path_escape_rejected(self):
        data = record(self.root)
        data['deliverables'][0]['uncharged_cif'] = str(self.root.parent / 'outside.cif')
        with self.assertRaises(api.ReleaseCurationError):
            self.validate(data)

    def test_network_access_forbidden(self):
        with self.assertRaisesRegex(RuntimeError, 'does not download'):
            worker.deny_network('example.org')

    def test_input_options_before_execution(self):
        for source_id in ('../escape', '', 'bad/name'):
            with self.assertRaises(ValueError):
                api.curate_cif_v41('unused', source_id=source_id, output_dir=self.root / 'out',
                                  python='unused', runtime_coremof_root='unused')
        with self.assertRaises(FileExistsError):
            api.curate_cif_v41('unused', source_id='SOURCE1', output_dir=self.root,
                              python='unused', runtime_coremof_root='unused')

    def run_fake(self, *, rc=0, timeout=False, skip=False):
        source = self.root / 'original.cif'
        source.write_text('data_original\n')
        source_bytes = source.read_bytes()
        runtime = self.root / 'runtime'
        runtime.mkdir()
        holder = {}

        class Process:
            returncode = rc
            pid = 99999

            def wait(self, timeout=None):
                if holder.get('timeout'):
                    holder['timeout'] = False
                    raise api.subprocess.TimeoutExpired('curation', timeout)
                return self.returncode

        def launch(command, **kwargs):
            request = json.loads(Path(command[-1]).read_text())
            self.assertEqual(kwargs['env']['CUDA_VISIBLE_DEVICES'], '')
            self.assertEqual(kwargs['env']['PYTHONHASHSEED'], '0')
            self.assertEqual(kwargs['env']['PYTHONPATH'], str(runtime))
            self.assertTrue(kwargs['start_new_session'])
            self.assertNotEqual(Path(request['input']), source)
            self.assertEqual(Path(request['input']).read_bytes(), source_bytes)
            root = Path(request['output'])
            root.mkdir()
            data = record(root, skip=skip)
            data['source_sha256'] = hashlib.sha256(source_bytes).hexdigest()
            Path(request['record']).write_text(json.dumps(data))
            holder['timeout'] = timeout
            return Process()

        real_check = worker.check

        def check(path, expected):
            path = Path(path)
            if path == source or path.name in {'protocol.py', 'candidate.cif'} or path == Path(api.protocol.__file__):
                return real_check(path, expected)
            return path

        with patch.object(worker, 'check', side_effect=check), patch.object(api.subprocess, 'Popen', side_effect=launch), patch.object(api, '_stop') as stopped:
            result = api.curate_cif_v41(source, source_id='SOURCE1', output_dir=self.root / 'result',
                                       python=self.root / 'python', runtime_coremof_root=runtime,
                                       skip_charges=skip, timeout_seconds=10)
        self.assertEqual(source.read_bytes(), source_bytes)
        self.assertEqual(stopped.called, timeout)
        self.assertFalse(any(self.root.glob('.curation-v41-*')))
        self.assertTrue((self.root / 'result' / 'receipt.json').is_file())
        return result

    def test_transactional_private_input_and_relative_artifacts(self):
        result = self.run_fake()
        self.assertEqual(result['execution_status'], 'COMPLETE')
        self.assertEqual(result['deliverables'][0]['uncharged_cif'], 'artifacts/candidate.cif')
        self.assertFalse(result['release_eligible'])

    def test_skip_charges_remains_ineligible(self):
        result = self.run_fake(skip=True)
        self.assertFalse(result['curation_stage_eligible'])

    def test_process_failure_retains_diagnostic_not_scientific_values(self):
        result = self.run_fake(rc=1)
        self.assertEqual(result['execution_status'], 'ERROR')
        self.assertEqual(result['deliverables'], [])

    def test_timeout_stops_own_group(self):
        self.assertEqual(self.run_fake(timeout=True)['execution_status'], 'TIMEOUT')


if __name__ == '__main__':
    unittest.main()
