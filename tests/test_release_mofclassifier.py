"""Dependency-minimal contract and process tests for recorded classifier replay."""
import json
import math
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from CoREMOF import release_mofclassifier as api

SID = 'FSR-COD-2016-0106'
DIGEST = 'a' * 64


def good_record(digest=DIGEST, score=0.6):
    scores = [score] * 100
    mean = math.fsum(scores) / 100
    return dict(structure_id=SID, canonical_cif_sha256=digest, method_id=api.PROFILE,
                asset_bundle_sha256=api.ASSET_BUNDLE_SHA256, method_config_sha256=api.worker.CONFIG_SHA256,
                canonical_post_run_sha256=digest, canonical_unchanged=True, model_family='core',
                checkpoint_bag_count=100, pass_threshold=0.6, positive_class_index=1,
                mean_aggregation='float64_math_fsum_divide_by_100', bag_scores=scores,
                completed_bag_count=100, execution_status='SUCCESS', private_copy_used=True,
                mean_score=mean, operational_vote='PASS' if mean >= 0.6 else 'FAIL',
                private_copy_initial_sha256=digest, private_copy_model_input_sha256=digest,
                private_copy_rewritten=False, error_type=None, error_message=None)


class RecordTests(unittest.TestCase):
    def test_frozen_protocol_and_configuration_bytes(self):
        self.assertEqual(api.sha(api.protocol.__file__), api.worker.PROTOCOL_SHA256)
        self.assertEqual(api.sha(Path(api.__file__).with_name('data') / 'mofclassifier_fresh_core_v1.json'),
                         api.worker.CONFIG_SHA256)

    def test_inclusive_threshold_and_exact_mean(self):
        for value in (0.0, math.nextafter(0.6, 0.0), 0.6, 1.0):
            record = good_record(score=value)
            api._validate(record, SID, DIGEST)
            self.assertEqual(record['operational_vote'], 'PASS' if record['mean_score'] >= 0.6 else 'FAIL')

    def test_missing_model_does_not_produce_success(self):
        record = good_record()
        record['bag_scores'].pop()
        record['completed_bag_count'] = 99
        with self.assertRaises(api.ReleaseMOFClassifierError):
            api._validate(record, SID, DIGEST)

    def test_nonfinite_boolean_out_of_range_scores_rejected(self):
        for value in (float('nan'), float('inf'), True, -0.1, 1.1, '0.6'):
            record = good_record()
            record['bag_scores'][0] = value
            with self.subTest(value=value), self.assertRaises(api.ReleaseMOFClassifierError):
                api._validate(record, SID, DIGEST)

    def test_mean_rounding_is_not_silently_accepted(self):
        record = good_record()
        record['mean_score'] += 1e-8
        with self.assertRaises(api.ReleaseMOFClassifierError):
            api._validate(record, SID, DIGEST)

    def test_fail_vote_is_not_an_execution_failure(self):
        record = good_record(score=0.2)
        api._validate(record, SID, DIGEST)
        self.assertEqual(record['execution_status'], 'SUCCESS')
        self.assertEqual(record['operational_vote'], 'FAIL')

    def test_unavailable_has_no_mean_or_vote(self):
        record = good_record()
        record.update(execution_status='ERROR', bag_scores=[0.3], completed_bag_count=1,
                      mean_score=None, operational_vote=None, error_type='MODEL_LOAD_ERROR')
        api._validate(record, SID, DIGEST)
        record['operational_vote'] = 'FAIL'
        with self.assertRaises(api.ReleaseMOFClassifierError):
            api._validate(record, SID, DIGEST)

    def test_private_reformat_can_never_mean_source_rewrite(self):
        record = good_record()
        record.update(private_copy_rewritten=True, private_copy_model_input_sha256='b' * 64)
        api._validate(record, SID, DIGEST)
        record['canonical_post_run_sha256'] = 'b' * 64
        with self.assertRaises(api.ReleaseMOFClassifierError):
            api._validate(record, SID, DIGEST)

    def test_identity_and_asset_mismatch_rejected(self):
        for changes in ({'structure_id': 'other'}, {'asset_bundle_sha256': 'b' * 64},
                        {'model_family': 'qsp'}, {'pass_threshold': 0.5}, {'execution_status': 'RUNNING'}):
            record = dict(good_record(), **changes)
            with self.subTest(changes=changes), self.assertRaises(api.ReleaseMOFClassifierError):
                api._validate(record, SID, DIGEST)


class ProcessTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.cif = self.root / 'source.cif'
        self.cif.write_text('data_source\n')
        self.digest = api.sha(self.cif)
        self.output = self.root / 'output'
        self.python = self.root / 'python'
        self.python.write_text('not executed by test\n')
        self.options = dict(structure_id=SID, output_dir=self.output, python=self.python,
                            model_root=self.root, timeout_seconds=2)
        self.mode = 'success'
        self.requests = []

    def fake_check(self, path, digest):
        if digest == api.worker.PYTHON_SHA256:
            return Path(path)
        return self.actual_check(path, digest)

    def fake_process(self, command, **kwargs):
        import csv
        request = json.loads(Path(command[-1]).read_text())
        self.requests.append((request, kwargs))
        with Path(request['manifest']).open() as stream:
            row = next(csv.DictReader(stream))
        private = Path(row['canonical_cif_path'])
        self.assertNotEqual(private, self.cif)
        self.assertEqual(private.read_bytes(), self.cif.read_bytes())
        self.assertTrue(kwargs['start_new_session'])
        self.assertNotIn('PYTHONPATH', kwargs['env'])
        self.assertEqual(kwargs['env']['CUDA_VISIBLE_DEVICES'], '')
        output = Path(request['output'])
        (output / 'records').mkdir()
        if self.mode != 'missing':
            (output / 'records' / ('00000000_' + SID + '.json')).write_text(json.dumps(good_record(self.digest)))
            Path(request['runtime_receipt']).write_text('{}')
        if self.mode == 'mutate':
            self.cif.write_text('modified source\n')
        def wait(timeout):
            if self.mode == 'timeout':
                raise subprocess.TimeoutExpired(command, timeout)
        return SimpleNamespace(returncode=2 if self.mode == 'error' else 0, wait=wait, pid=123456)

    def run_api(self):
        self.actual_check = api.check
        with patch.object(api, 'check', side_effect=self.fake_check), \
                patch.object(api.protocol, 'load_method_contract', return_value=SimpleNamespace(asset_bundle_sha256=api.ASSET_BUNDLE_SHA256)), \
                patch.object(api.subprocess, 'Popen', side_effect=self.fake_process), patch.object(api, '_stop'):
            return api.calculate_release_mofclassifier(self.cif, **self.options)

    def test_source_safe_transactional_success(self):
        result = self.run_api()
        self.assertEqual(result['execution_status'], 'SUCCESS')
        self.assertEqual(result['completed_bag_count'], 100)
        self.assertFalse(result['source_cif_modified'])
        self.assertFalse(result['release_labels_updated'])
        self.assertFalse(result['mechanistic_hard_fail_eligible'])
        self.assertEqual(api.sha(self.cif), self.digest)
        self.assertTrue((self.output / 'receipt.json').is_file())
        self.assertFalse(Path(self.requests[0][0]['manifest']).exists())

    def test_timeout_has_no_invented_score(self):
        self.mode = 'timeout'
        result = self.run_api()
        self.assertEqual(result['execution_status'], 'TIMEOUT')
        self.assertIsNone(result['mean_score'])
        self.assertIsNone(result['operational_vote'])

    def test_process_error_does_not_use_partial_success_record(self):
        self.mode = 'error'
        result = self.run_api()
        self.assertEqual(result['execution_status'], 'PROCESS_ERROR')
        self.assertIsNone(result['operational_vote'])

    def test_missing_record_remains_unavailable(self):
        self.mode = 'missing'
        self.assertEqual(self.run_api()['execution_status'], 'PROCESS_ERROR')

    def test_changed_original_prevents_publication(self):
        self.mode = 'mutate'
        with self.assertRaises(RuntimeError):
            self.run_api()
        self.assertFalse(self.output.exists())

    def test_existing_output_is_preserved(self):
        self.output.mkdir()
        with self.assertRaises(FileExistsError):
            self.run_api()
        self.assertFalse(self.requests)

    def test_invalid_inputs_prevent_runtime_launch(self):
        for change in ({'structure_id': '../escape'}, {'memory_limit_mb': True},
                       {'timeout_seconds': 0}, {'memory_limit_mb': 32}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                api.calculate_release_mofclassifier(self.cif, **dict(self.options, **change))


if __name__ == '__main__':
    unittest.main()
