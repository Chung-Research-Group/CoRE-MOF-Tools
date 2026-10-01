"""Small, target-free tests for the explicit frozen-assignment replay example."""
import hashlib
import importlib.util
import io
from pathlib import Path
import tarfile
import tempfile
import unittest


SOURCE = Path(__file__).resolve().parents[1] / 'examples/replay_common_input_benchmark.py'
spec = importlib.util.spec_from_file_location('common_input_replay_example', SOURCE)
example = importlib.util.module_from_spec(spec)
spec.loader.exec_module(example)


def population():
    metadata, groups, topology = {}, {}, {}
    for i in range(120):
        sid = f'structure-{i:03d}'
        metadata[sid] = {'label': 'CR' if i < 100 else 'NCR', 'diversity_tier': 'test',
                         'diversity_stratum': str(i % 4)}
        groups[sid] = {'effective_leakage_block': sid, 'source_database': 'COD', 'structure_variant': 'FSR'}
        topology[sid] = {'topology_available': 'false'}
    return metadata, groups, topology, {f'structure-{i:03d}' for i in range(10)}


def frozen_rows(metadata, groups):
    result = []
    for qkey, q, count in (('0', '0', 0), ('0p5', '0.5', 10), ('1', '1', 20)):
        for i in list(range(100-count)) + list(range(100, 100+count)):
            sid = f'structure-{i:03d}'
            result.append(dict(run_key='seed912_q'+qkey, seed='912',
                requested_ncr_pool_fraction=q, actual_ncr_ratio=str(count/100),
                structure_id=sid, label=metadata[sid]['label'],
                partition='test' if i < 10 else 'validation' if i < 20 else 'train',
                effective_leakage_block=groups[sid]['effective_leakage_block'],
                diversity_tier='test', diversity_stratum=str(i%4)))
    return result


class CommonInputReplayTests(unittest.TestCase):
    def test_assignments_ignore_targets_and_are_order_independent(self):
        metadata, groups, topology, fixed = population()
        archived = frozen_rows(metadata, groups)
        first = example.validate_frozen_assignments(archived, metadata, groups, fixed, seeds=(912,))
        for row in metadata.values():
            row['target'] = object()  # Must never be serialized or consumed.
        second = example.validate_frozen_assignments(archived, dict(reversed(list(metadata.items()))), groups, fixed, seeds=(912,))
        self.assertEqual(first, second)
        self.assertEqual(len(first), 300)
        for q in ('0', '0.5', '1'):
            run = [row for row in first if row['requested_ncr_pool_fraction'] == q]
            self.assertEqual(example.Counter(row['partition'] for row in run), {'train': 80, 'validation': 10, 'test': 10})
            self.assertEqual({r['structure_id'] for r in run if r['partition'] == 'test'}, fixed)
        by_id = {}
        for row in first:
            by_id.setdefault(row['structure_id'], set()).add(row['partition'])
        self.assertTrue(all(len(parts) == 1 for parts in by_id.values()))

    def test_renaming_retains_membership_and_row_order(self):
        metadata, groups, topology, fixed = population()
        archived = frozen_rows(metadata, groups)
        mapping = {sid: f'2026[Cu][nan]3[ASR]{120-i}' for i, sid in enumerate(metadata)}
        renamed = [dict(row, structure_id=mapping[row['structure_id']]) for row in archived]
        new_metadata = {mapping[sid]: value for sid, value in metadata.items()}
        new_groups = {mapping[sid]: value for sid, value in groups.items()}
        rows = example.validate_frozen_assignments(renamed, new_metadata, new_groups,
                                                  {mapping[sid] for sid in fixed}, seeds=(912,))
        self.assertEqual([r['structure_id'] for r in rows], [r['structure_id'] for r in renamed])
        self.assertEqual([r['partition'] for r in rows], [r['partition'] for r in archived])

    def test_changed_assignment_rejected(self):
        metadata, groups, topology, fixed = population()
        archived = frozen_rows(metadata, groups)
        archived[25]['partition'] = 'validation'
        with self.assertRaisesRegex(ValueError, 'partition sizes'):
            example.validate_frozen_assignments(archived, metadata, groups, fixed, seeds=(912,))

    def test_shared_test_group_is_rejected(self):
        metadata, groups, topology, fixed = population()
        groups['structure-020']['effective_leakage_block'] = groups['structure-000']['effective_leakage_block']
        with self.assertRaisesRegex(ValueError, 'fixed-test group'):
            example.validate_frozen_assignments(frozen_rows(metadata, groups), metadata, groups, fixed, seeds=(912,))

    def test_unchecked_label_is_not_treated_as_ncr(self):
        metadata, groups, topology, fixed = population()
        metadata['structure-105']['label'] = 'UNCHECKED'
        with self.assertRaisesRegex(ValueError, 'invalid label'):
            example.validate_frozen_assignments(frozen_rows(metadata, groups), metadata, groups, fixed, seeds=(912,))

    def make_archive(self, root, corrupt=False, duplicate=False, link=False,
                     prefix=None, dataset_id=None):
        archive = root / 'handoff.tar.gz'
        prefix = example.PREFIX if prefix is None else prefix
        files = {key: name.replace(example.DATASET_ID, dataset_id or example.DATASET_ID)
                 for key, name in example.FILES.items()}
        data = b'private input fixture\n'
        expected = hashlib.sha256(data).hexdigest()
        ledger = ''.join(expected + '  ' + name + '\n' for name in files.values()).encode()
        with tarfile.open(archive, 'w:gz') as bundle:
            for index, name in enumerate(files.values()):
                contents = b'altered' if corrupt and index == 0 else data
                item = tarfile.TarInfo(prefix + name)
                item.size = len(contents)
                if link and index == 0:
                    item.type = tarfile.SYMTYPE
                    item.linkname = '/etc/passwd'
                    item.size = 0
                    bundle.addfile(item)
                else:
                    bundle.addfile(item, io.BytesIO(contents))
                    if duplicate and index == 0:
                        bundle.addfile(item, io.BytesIO(contents))
            item = tarfile.TarInfo(prefix + 'SHA256SUMS')
            item.size = len(ledger)
            bundle.addfile(item, io.BytesIO(ledger))
        return archive, example.sha(archive)

    def test_archive_input_read_checks_hashes_without_extracting(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive, digest = self.make_archive(root)
            inputs, hashes = example.load_inputs(archive, digest)
            self.assertEqual(set(inputs), set(example.FILES))
            self.assertEqual(set(hashes), set(example.FILES.values()))
            self.assertEqual(list(root.iterdir()), [archive])
            with self.assertRaisesRegex(ValueError, 'Archive SHA-256 differs'):
                example.load_inputs(archive, '0' * 64)

    def test_archive_roles_accept_renamed_handoff_and_dataset_directories(self):
        with tempfile.TemporaryDirectory() as directory:
            archive, digest = self.make_archive(Path(directory), prefix='fixture-handoff/',
                                                dataset_id='fixture-common-input')
            inputs, hashes = example.load_inputs(archive, digest)
            self.assertEqual(set(inputs), set(example.FILES))
            self.assertTrue(any('/fixture-common-input/' in name for name in hashes))
            self.assertFalse(any(example.DATASET_ID in name for name in hashes))

    def test_corrupt_duplicate_or_linked_input_fails_closed(self):
        for fault in ('corrupt', 'duplicate', 'link'):
            with self.subTest(fault=fault), tempfile.TemporaryDirectory() as directory:
                archive, digest = self.make_archive(Path(directory), **{fault: True})
                with self.assertRaises(ValueError):
                    example.load_inputs(archive, digest)

    def test_existing_output_is_untouched_before_any_input_read(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'existing'
            output.mkdir()
            (output / 'saved').write_bytes(b'unchanged')
            with self.assertRaises(FileExistsError):
                example.replay(Path(directory) / 'nonexistent.tar.gz', '0' * 64, output)
            self.assertEqual((output / 'saved').read_bytes(), b'unchanged')


if __name__ == '__main__':
    unittest.main()
