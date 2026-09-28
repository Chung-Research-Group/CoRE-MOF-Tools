import copy
import csv
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from CoREMOF import attach_targets
from CoREMOF.attachments import frozen_assignment_manifest
from CoREMOF.benchmarks import BenchmarkFeasibilityError, available_group_criteria, build_diversity_index
from CoREMOF.dataset import CoREMOFDataset, ReleaseValidationError
from CoREMOF.labels import CHECKER_COLUMNS
from CoREMOF.parents import ParentResolver, ZEO_NUMERIC_FINGERPRINT_DEFINITION
from CoREMOF.projections import export_source_projection, projection_context
from CoREMOF.targets import TargetSource, TargetDataError, _read_feature_table
from test_benchmarks import EXACT_BENCHMARK_BACKEND
from test_dataset_labels import _make_release, _metadata_rows, _parent_rows, _write_csv


def _fixture(root):
    root.mkdir()
    _make_release(root)
    metadata, parents, manifests = [], [], []
    template = _parent_rows()[0]
    for index in range(25):
        source = 'COD' if index < 24 else 'CSD'
        sid = 'ASR-{}-2026-{:04d}'.format(source, index + 1)
        label = 'CR' if index < 20 else 'NCR'
        row = dict(_metadata_rows()[0], structure_id=sid, source_database=source,
                   source_id='SOURCE-' + sid, cif_file='cifs/' + sid + '.cif',
                   label_3checker=label, label_4checker=label, label_5checker=label)
        row.update({column: 'PASS' if label == 'CR' else 'FAIL' for column in CHECKER_COLUMNS.values()})
        metadata.append(row)
        parent = {'structure_id': sid}
        for key in template:
            if key.endswith('_group'):
                prefix = key[:-6]
                parent[prefix + '_group'] = template[key].split('-')[0] + '-{:08X}'.format(index + 100)
                parent[prefix + '_status'] = 'UNMATCHED'
                parent[prefix + '_size'] = '1'
        parents.append(parent)
        content = ('data_' + sid + '\n').encode()
        (root / 'cifs').mkdir(exist_ok=True)
        (root / row['cif_file']).write_bytes(content)
        manifests.append({'structure_id': sid, 'cif_file': row['cif_file'],
                          'sha256': hashlib.sha256(content).hexdigest(), 'size_bytes': str(len(content))})
    # Two selected CR records are connected through one omitted CSD NCR record.
    for prefix, indices in (('rac', (0, 24)), ('mofid1', (1, 24))):
        group = parents[indices[0]][prefix + '_group']
        for index in indices:
            parents[index].update({prefix + '_group': group, prefix + '_status': 'MATCHED', prefix + '_size': '2'})
    _write_csv(root / 'metadata/metadata.csv', tuple(metadata[0]), metadata)
    _write_csv(root / 'parent_groups/parent_groups.csv', tuple(parents[0]), parents)
    _write_csv(root / 'manifests/cif_manifest.csv', tuple(manifests[0]), manifests)
    info = json.loads((root / 'dataset_info.json').read_text())
    info['structure_count'] = len(metadata)
    (root / 'dataset_info.json').write_text(json.dumps(info))
    return CoREMOFDataset.from_release(root, verify_cif_files=True)


def _slice(dataset, root):
    root.mkdir()
    selected = {row['structure_id'] for row in dataset.metadata_rows if row['source_database'] == 'COD'}
    for name in ('metadata/metadata.csv', 'parent_groups/parent_groups.csv', 'manifests/cif_manifest.csv'):
        with (dataset.release_root / name).open(newline='') as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames
            rows = [row for row in reader if row['structure_id'] in selected]
        _write_csv(root / name, fields, rows)
    shutil.copy2(dataset.release_root / 'parent_groups/parent_group_methods.json', root / 'parent_groups/parent_group_methods.json')
    (root / 'source_bundle.json').write_text('{"standalone_coremof_release":false}')
    (root / 'cifs').mkdir()
    for sid in selected:
        shutil.copy2(dataset.release_root / ('cifs/' + sid + '.cif'), root / ('cifs/' + sid + '.cif'))


def _features(dataset, selected_root):
    ids = dataset.structure_ids
    rac_fields = tuple('rac_{:03d}'.format(index) for index in range(264))
    zeo_fields = tuple(ZEO_NUMERIC_FINGERPRINT_DEFINITION['numeric_fields']) + (
        'n2_channel_dimension', 'structure_periodic_dimension')
    tables = {
        'rac5_features.csv': [dict(structure_id=sid, rac5_available=str(index < 12).lower(),
            **{field: str((index + 1) * (column % 7)) if index < 12 else ''
               for column, field in enumerate(rac_fields)}) for index, sid in enumerate(ids)],
        'zeo_features.csv': [dict(structure_id=sid, n2_he_available=str(index < 20).lower(),
            periodicity_available=str(index < 20).lower(),
            **{field: str(index + 1) if index < 20 else '' for field in zeo_fields})
            for index, sid in enumerate(ids)],
        'topology_features.csv': [dict(structure_id=sid, topology_available='true',
            network_dimension='3', single_node_net='pcu', all_node_net='pcu', single_all_agree='true')
            for sid in ids],
    }
    for name, rows in tables.items():
        _write_csv(dataset.release_root / 'features' / name, tuple(rows[0]), rows)
        _write_csv(selected_root / 'features' / name, tuple(rows[0]),
                   [row for row in rows if row['structure_id'] in ids[:-1]])


class SourceProjectionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.full = _fixture(self.root / 'full')
        self.slice = self.root / 'COD'
        _slice(self.full, self.slice)
        self.path = self.root / 'cod_projection.json'
        self.report = export_source_projection(self.full, self.slice, self.path, sources=('COD',),
            group_profiles=(('priority_main',), ('rac5', 'mofid_v1')), diversity='none')
        self.projected = self.load()

    def tearDown(self):
        self.temp.cleanup()

    def load(self, **kwargs):
        return CoREMOFDataset.from_projection(self.slice, self.path,
            expected_sha256=kwargs.pop('expected_sha256', self.report['sha256']), **kwargs)

    def changed_contract(self, modify):
        value = json.loads(self.path.read_text())
        modify(value)
        path = self.root / 'modified_contract.json'
        path.write_text(json.dumps(value))
        return CoREMOFDataset.from_projection(self.slice, path, expected_sha256=hashlib.sha256(path.read_bytes()).hexdigest())

    def test_classification_parent_scope_and_cif_verification(self):
        self.assertEqual(len(self.projected), 24)
        self.assertEqual(self.projected.classify(5).label_counts()['CR'], 20)
        self.assertEqual(self.projected[self.projected.structure_ids[0]].parent_group('rac5').size, 2)
        self.assertTrue(self.load(verify_cif_files=True).cif_files_verified)
        self.assertEqual(projection_context(self.projected)['complete_release']['structure_count'], 25)

    def test_full_release_loader_remains_strict(self):
        with self.assertRaisesRegex(ReleaseValidationError, 'source-separated review bundle'):
            CoREMOFDataset.from_release(self.slice)
        self.assertFalse((self.slice / 'dataset_info.json').exists())

    def test_full_source_general_assignments_and_hidden_bridge_are_exact(self):
        expected = self.full.classify(5).data_split(sources=('COD',), labels=None, diversity='none', random_state=42)
        actual = self.projected.classify(5).data_split(labels=None, diversity='none', random_state=42)
        self.assertEqual(dict(actual.assignments), dict(expected.assignments))
        first, second = self.projected.structure_ids[:2]
        self.assertEqual(actual.effective_leakage_blocks[first], actual.effective_leakage_blocks[second])
        self.assertTrue(actual.leakage_audit['passed'])
        self.assertEqual(actual.receipt()['source_projection']['complete_release_structure_count'], 25)
        self.assertIn('projections.py', actual.receipt()['implementation']['source_sha256'])

    def test_no_omitted_structure_or_machine_path_in_contract(self):
        text = self.path.read_text()
        self.assertNotIn('ASR-CSD-2026-0025', text)
        self.assertNotIn(str(self.root), text)

    def test_group_profiles_cannot_be_recombined_or_reordered(self):
        available = available_group_criteria(self.projected)
        self.assertTrue(available['priority_main']['available'])
        self.assertFalse(available['zeo']['available'])
        for profile in (('rac5',), ('mofid_v1', 'rac5'), ('priority_main', 'rac5')):
            with self.subTest(profile=profile), self.assertRaisesRegex(ReleaseValidationError, 'ordered grouping profile'):
                self.projected.classify(5).data_split(group_criteria=profile, diversity='none')
        valid = self.projected.classify(5).data_split(group_criteria=('rac5', 'mofid_v1'), diversity='none')
        self.assertTrue(valid.leakage_audit['passed'])

    def test_legacy_parent_reconstruction_is_explicitly_rejected(self):
        with self.assertRaisesRegex(ValueError, 'omitted-source relationships'):
            ParentResolver(self.projected)

    def test_missing_cached_diversity_has_no_silent_fallback(self):
        with self.assertRaisesRegex(ReleaseValidationError, 'no verified representative'):
            self.projected.classify(5).data_split()

    def test_cohorts_check_hidden_member_labels_and_match_full_eligibility(self):
        with self.assertRaisesRegex(BenchmarkFeasibilityError, 'another label'):
            self.projected.classify(5).build_cr_ncr_cohorts(diversity='none')
        options = dict(ncr_pool_fractions=(0, 0.5, 1), seeds=(42, 43), diversity='none',
            cohort_eligibility='complete_release_label_pure_effective_blocks',
            partition_strategy='transition_balanced')
        actual = self.projected.classify(5).build_cr_ncr_benchmark(**options)
        expected = self.full.classify(5).build_cr_ncr_benchmark(**options,
            eligible_structure_ids=self.projected.structure_ids)
        self.assertEqual(actual.fixed_test_ids, expected.fixed_test_ids)
        self.assertEqual([dict(run.assignments) for run in actual.runs], [dict(run.assignments) for run in expected.runs])
        receipt = actual.receipt()['cohort_receipt']
        self.assertEqual(receipt['eligible_pool_counts'], {'C_CR': 18, 'M_NCR': 4})
        self.assertEqual(receipt['complete_release_label_accounting']['counts']['NCR'], 5)
        self.assertEqual(receipt['selected_source_label_accounting']['counts']['NCR'], 4)
        self.assertFalse(set(self.projected.structure_ids[:2]).intersection(actual.fixed_test_ids))
        actual.write(self.root / 'suite')

    def test_wrong_checksum_and_selected_table_drift_fail(self):
        with self.assertRaisesRegex(ReleaseValidationError, 'checksum mismatch'):
            self.load(expected_sha256='0' * 64)
        path = self.slice / 'metadata/metadata.csv'
        path.write_bytes(path.read_bytes() + b'\n')
        with self.assertRaisesRegex(ReleaseValidationError, 'input checksum mismatch'):
            self.load()

    def test_cif_drift_is_rejected_when_requested(self):
        path = self.slice / ('cifs/' + self.projected.structure_ids[0] + '.cif')
        path.write_bytes(b'changed CIF')
        with self.assertRaises(ReleaseValidationError):
            self.load(verify_cif_files=True)

    def test_no_overwrite_or_projection_of_a_projection(self):
        before = self.path.read_bytes()
        with self.assertRaises(FileExistsError):
            export_source_projection(self.full, self.slice, self.path, sources=('COD',), diversity='none')
        self.assertEqual(before, self.path.read_bytes())
        with self.assertRaisesRegex(ReleaseValidationError, 'complete release'):
            export_source_projection(self.projected, self.slice, self.root / 'new.json', sources=('COD',), diversity='none')

    def test_impossible_label_counts_fail_even_with_new_checksum(self):
        def change(value):
            group = value['profiles'][0]['effective_groups'][self.projected.structure_ids[0]]
            value['profiles'][0]['complete_group_label_counts'][group]['CR'] = 1
        with self.assertRaisesRegex(ReleaseValidationError, 'omit selected members'):
            self.changed_contract(change)

    def test_boolean_count_and_unauthorized_publication_fail(self):
        for field, value in (('publication_authorized', True), ('official_split', True)):
            with self.subTest(field=field), self.assertRaises(ReleaseValidationError):
                self.changed_contract(lambda item: item.update({field: value}))
        with self.assertRaises(ReleaseValidationError):
            self.changed_contract(lambda item: item['complete_release'].update(structure_count=True))

    def test_invalid_parent_group_sizes_fail(self):
        def change(value):
            sizes = value['parent_group_sizes']['rac']
            group = self.projected[self.projected.structure_ids[0]].parent_group('rac5').group_id
            sizes[group] = 1
        with self.assertRaisesRegex(ReleaseValidationError, 'declares size'):
            self.changed_contract(change)

    def test_copy_and_private_attribute_replacement_do_not_transfer_authority(self):
        with self.assertRaises(TypeError):
            copy.copy(self.projected)
        object.__setattr__(self.projected, '_authority_extra_state', {})
        with self.assertRaises(ValueError):
            self.projected.classify(5)

    def test_target_attachment_preserves_zero_and_frozen_assignments(self):
        split = self.projected.classify(5).data_split(labels=None, diversity='none')
        target_path = self.root / 'targets.csv'
        _write_csv(target_path, ('structure_id', 'value'), [{'structure_id': sid, 'value': 0}
                                                          for sid in self.projected.structure_ids])
        source = TargetSource(target_path, target_columns=('value',), value_types={'value': 'float'})
        before = split.assignment_digest
        attached = split.attach_targets(source, missing='error')
        self.assertEqual(before, split.assignment_digest)
        frozen = frozen_assignment_manifest(split.assignment_rows(), split.receipt())
        saved = attach_targets(frozen, source, dataset=self.projected, missing='error')
        self.assertEqual(attached.rows(), saved.rows())
        self.assertTrue(all(row['value'] == 0 for row in attached.rows()))

    def test_successful_export_removes_only_its_temporary_staging(self):
        self.assertFalse(list(self.root.glob('.coremof-projection-*')))
        self.assertTrue(Path(str(self.path) + '.sha256').is_file())

    def test_inconsistent_main_union_and_parent_method_bindings_fail(self):
        def change(value):
            value['profiles'][1]['main_union_groups'][self.projected.structure_ids[3]] = 'FAKE'
        with self.assertRaisesRegex(ReleaseValidationError, 'main-union groups differ'):
            self.changed_contract(change)
        with self.assertRaisesRegex(ReleaseValidationError, 'parent methods differ'):
            self.changed_contract(lambda value: value['complete_release']['input_sha256'].update(
                {'parent_groups/parent_group_methods.json': '0' * 64}))

    def test_duplicate_keys_nonfinite_and_extra_contract_fields_fail(self):
        for content in ('{"schema_version":1,"schema_version":2}', '{"value":NaN}'):
            path = self.root / 'invalid.json'
            path.write_text(content)
            with self.assertRaises(ReleaseValidationError):
                CoREMOFDataset.from_projection(self.slice, path,
                    expected_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        with self.assertRaisesRegex(ReleaseValidationError, 'source projection fields'):
            self.changed_contract(lambda value: value.update(unapproved='field'))

    def test_unbound_features_cannot_be_read_after_loading(self):
        _features(self.full, self.slice)
        with self.assertRaisesRegex(TargetDataError, 'verified projection contract'):
            _read_feature_table(self.projected, 'rac5')

    def test_row_order_does_not_change_assignments(self):
        for relative in ('metadata/metadata.csv', 'parent_groups/parent_groups.csv', 'manifests/cif_manifest.csv'):
            path = self.slice / relative
            with path.open(newline='') as handle:
                reader = csv.DictReader(handle)
                fields, rows = reader.fieldnames, list(reader)
            _write_csv(path, fields, list(reversed(rows)))
        output = self.root / 'reordered.json'
        result = export_source_projection(self.full, self.slice, output, sources=('COD',),
            group_profiles=('priority_main',), diversity='none')
        other = CoREMOFDataset.from_projection(self.slice, output, expected_sha256=result['sha256'])
        options = dict(labels=None, diversity='none', random_state=43)
        self.assertEqual(dict(other.classify(5).data_split(**options).assignments),
                         dict(self.projected.classify(5).data_split(**options).assignments))

    @unittest.skipUnless(EXACT_BENCHMARK_BACKEND, 'requires the pinned numerical benchmark environment')
    def test_saved_full_release_strata_and_target_free_tiers_are_reused_without_refitting(self):
        _features(self.full, self.slice)
        path = self.root / 'representative.json'
        result = export_source_projection(self.full, self.slice, path, sources=('COD',),
            group_profiles=('priority_main',), diversity='representative')
        projected = CoREMOFDataset.from_projection(self.slice, path, expected_sha256=result['sha256'])
        full_index = build_diversity_index(self.full)
        expected = self.full.classify(5).data_split(sources=('COD',), labels=None)
        with patch('CoREMOF.benchmarks._load_backend', side_effect=AssertionError('must not refit')):
            index = build_diversity_index(projected)
            actual = projected.classify(5).data_split(labels=None)
        self.assertEqual(dict(index.strata_by_id), {sid: full_index.strata_by_id[sid] for sid in projected.structure_ids})
        self.assertEqual(set(index.tier_by_id.values()), {'rac5', 'zeo', 'no_numeric'})
        self.assertEqual(index.profile['complete_release_index_sha256'], full_index.digest)
        self.assertFalse(index.profile['scientific_feature_imputation'])
        self.assertEqual(dict(actual.assignments), dict(expected.assignments))
        feature = self.slice / 'features/rac5_features.csv'
        feature.write_bytes(feature.read_bytes() + b'\n')
        with self.assertRaisesRegex(TargetDataError, 'verified projection contract'):
            _read_feature_table(projected, 'rac5')


if __name__ == '__main__':
    unittest.main()
