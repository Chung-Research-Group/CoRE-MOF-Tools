#!/usr/bin/env python3
"""Reproduce the frozen September 12 common-input benchmark assignments.

This is an assignment replay, not a new dataset or scientific recalculation.
It uses archived target-free groups, diversity strata, input-availability
decisions and the original clean test. Targets and predictions are not read.
No code from the input archive is executed and no archive files are extracted.
"""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import shutil
import tarfile
import tempfile
from types import SimpleNamespace

from CoREMOF import benchmarks as algorithm
from CoREMOF._transactions import publish_directory


DATASET_ID = 'coremof_cod_common_input_equal_size_20260912_v1'
ASSIGNMENT_SHA256 = '9e72992970518d039f9631b1f45b516f4ff3603f7945bdcb82dbad283529dcbd'
PREFIX = 'coremof_cod_workflow_handoff_20260913_v1/'
BASE = 'workstation_files/CoREMOF-COD_benchmark_transfer_20260912_v1/dataset/'
DERIVED = 'workstation_files/' + DATASET_ID + '/'
FILES = {
    'groups': BASE + 'full_release_grouping.csv.gz',
    'topology': BASE + 'input_view/features/topology_features.csv',
    'eligibility': BASE + 'eligibility_manifest.csv',
    'original_membership': BASE + 'suite/coremof_cr_ncr_benchmark/membership_manifest.csv',
    'original_receipt': BASE + 'suite/coremof_cr_ncr_benchmark/receipt.json',
    'input_status': DERIVED + 'source_population_model_status.csv',
    'expected_membership': DERIVED + 'membership_manifest.csv',
    'expected_receipt': DERIVED + 'receipt.json',
}


def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(chunk)
    return result.hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def load_inputs(path, expected_sha256):
    """Read one hash-bound handoff by its file roles, independent of edition names.

    Require one top-level directory, one source dataset and one derived dataset.
    Canonical paths, sizes, unique regular members and the original checksum
    ledger remain authoritative; no members are extracted or executed.
    """
    if not re.fullmatch('[0-9a-f]{64}', expected_sha256):
        raise ValueError('Supply the independently verified archive SHA-256')
    if sha(path) != expected_sha256:
        raise ValueError('Archive SHA-256 differs')
    source_suffixes = {key: value[len(BASE):] for key, value in FILES.items()
                       if value.startswith(BASE)}
    derived_suffixes = {key: value[len(DERIVED):] for key, value in FILES.items()
                        if value.startswith(DERIVED)}
    payloads = {}
    total = 0
    with tarfile.open(path, 'r:gz') as archive:
        for member in archive:
            parts = Path(member.name).parts
            if not parts or parts[0] in ('.', '..'):
                continue
            source_match = (len(parts) >= 5 and parts[1] == 'workstation_files'
                            and parts[3] == 'dataset'
                            and '/'.join(parts[4:]) in source_suffixes.values())
            derived_match = (len(parts) == 4 and parts[1] == 'workstation_files'
                             and parts[3] in derived_suffixes.values())
            ledger_match = len(parts) == 2 and parts[1] == 'SHA256SUMS'
            if not (source_match or derived_match or ledger_match):
                continue
            if (Path(member.name).is_absolute() or '..' in parts
                    or Path(member.name).as_posix() != member.name
                    or member.name in payloads or not member.isfile()):
                raise ValueError('Duplicate, unsafe or non-regular replay input: ' + member.name)
            total += member.size
            if member.size > 40 * 1024**2 or total > 100 * 1024**2:
                raise ValueError('Replay input exceeds its bounded size contract')
            payloads[member.name] = archive.extractfile(member).read()
    candidates = [name for name in payloads if name.endswith('/dataset/full_release_grouping.csv.gz')]
    statuses = [name for name in payloads if name.endswith('/source_population_model_status.csv')]
    if len(candidates) != 1 or len(statuses) != 1:
        raise ValueError('Required replay inputs are absent or ambiguous')
    source_root = candidates[0][:-len('full_release_grouping.csv.gz')]
    derived_root = statuses[0][:-len('source_population_model_status.csv')]
    prefix = source_root.split('workstation_files/', 1)[0]
    if derived_root.split('workstation_files/', 1)[0] != prefix:
        raise ValueError('Replay source and derived inputs have different archive roots')
    actual_files = {key: source_root + suffix for key, suffix in source_suffixes.items()}
    actual_files.update({key: derived_root + suffix for key, suffix in derived_suffixes.items()})
    wanted = set(actual_files.values()) | {prefix + 'SHA256SUMS'}
    if set(payloads) != wanted:
        raise ValueError('Required replay inputs are absent or ambiguous')
    ledger = {}
    for line in payloads[prefix + 'SHA256SUMS'].decode().splitlines():
        expected, name = line.split('  ', 1)
        if name in ledger:
            raise ValueError('Duplicate archive ledger entry')
        ledger[name] = expected
    hashes = {}
    for name in actual_files.values():
        logical_name = name[len(prefix):]
        value = hashlib.sha256(payloads[name]).hexdigest()
        if ledger.get(logical_name) != value:
            raise ValueError('Replay input fails its package ledger: ' + name)
        hashes[logical_name] = value
    if sha(path) != expected_sha256:
        raise ValueError('Archive changed during the replay-input read')
    return {key: payloads[name] for key, name in actual_files.items()}, hashes


def rows(data, compressed=False):
    if compressed:
        with gzip.GzipFile(fileobj=io.BytesIO(data)) as stream:
            data = stream.read(40 * 1024**2 + 1)
        if len(data) > 40 * 1024**2:
            raise ValueError('Decompressed grouping table exceeds its size limit')
    return list(csv.DictReader(io.StringIO(data.decode('utf-8-sig'))))


def indexed(table):
    result = {row['structure_id']: row for row in table}
    if len(result) != len(table):
        raise ValueError('Duplicate structure ID in a one-row-per-structure input')
    return result


def reconstruct(metadata, full_groups, topology, fixed_test, seeds=(912, 913, 914, 915)):
    """Use the current package's exact group sampler with saved strata/test IDs."""
    cr = tuple(sorted(s for s, r in metadata.items() if r['label'] == 'CR'))
    ncr = tuple(sorted(s for s, r in metadata.items() if r['label'] == 'NCR'))
    if len(cr) + len(ncr) != len(metadata) or not fixed_test <= set(cr):
        raise ValueError('The selected population or fixed test has an invalid label')
    blocks = {s: full_groups[s]['effective_leakage_block'] for s in metadata}
    test_blocks = {blocks[s] for s in fixed_test}
    if any(blocks[s] in test_blocks for s in set(metadata) - fixed_test):
        raise ValueError('A fixed-test group also contains a train/validation candidate')
    strata = {}
    for s, row in metadata.items():
        top = topology[s]
        category = ('|'.join(str(top.get(k, '')).strip() or '<MISSING>' for k in
                    ('network_dimension', 'single_node_net', 'all_node_net', 'single_all_agree'))
                    if top['topology_available'].lower() == 'true' else '<NO_CURRENT_SUCCESS>')
        strata[s] = (row['label'], full_groups[s]['source_database'], full_groups[s]['structure_variant'],
                     category, row['diversity_tier'], row['diversity_stratum'])
    counts = (0, (len(ncr) + 1) // 2, len(ncr))
    validation_count = len(cr) - (8 * len(cr) + 5) // 10 - len(fixed_test)
    remaining = tuple(sorted(set(cr) - fixed_test))
    pool = tuple(sorted(set(metadata) - fixed_test))
    output = []
    for seed in seeds:
        adds = algorithm._nested_exact_block_memberships(
            ncr, blocks, strata, counts, seed, 'ncr-whole-block-priority', 'NCR addition')
        removes = algorithm._nested_exact_block_memberships(
            remaining, blocks, strata, counts, seed, 'cr-whole-block-removal-priority', 'CR removal')
        levels = []
        for q, n in zip(('0', '0.5', '1'), counts):
            retained = tuple(sorted(set(cr) - set(removes[n])))
            levels.append(SimpleNamespace(requested_ncr_pool_fraction=q, cr_ids=retained,
                ncr_ids=adds[n], structure_ids=tuple(sorted(retained + adds[n]))))
        partitions = algorithm._assign_transition_balanced_blocks(pool, levels, blocks, strata, seed, validation_count)
        for qkey, level in zip(('0', '0p5', '1'), levels):
            run = []
            for s in level.structure_ids:
                row = metadata[s]
                run.append(dict(run_key=f'seed{seed}_q{qkey}', seed=seed,
                    requested_ncr_pool_fraction=level.requested_ncr_pool_fraction,
                    actual_ncr_ratio=len(level.ncr_ids) / len(cr), structure_id=s, label=row['label'],
                    partition='test' if s in fixed_test else partitions[blocks[s]],
                    effective_leakage_block=blocks[s], diversity_tier=row['diversity_tier'], diversity_stratum=row['diversity_stratum']))
            expected = {'train': len(cr) - len(fixed_test) - validation_count,
                        'validation': validation_count, 'test': len(fixed_test)}
            if Counter(row['partition'] for row in run) != expected:
                raise ValueError('The group sampler did not reproduce equal partition sizes')
            output.extend(run)
    return output


def replay(path, expected_sha256, output):
    output = Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    code_paths = [Path(__file__).resolve(), Path(algorithm.__file__).resolve()]
    code_hashes = {str(p): sha(p) for p in code_paths}
    payloads, hashes = load_inputs(path, expected_sha256)
    receipt = json.loads(payloads['expected_receipt'])
    if receipt['dataset_id'] != next(Path(name).parent.name for name in hashes
            if name.endswith('/source_population_model_status.csv')) or receipt['suite_assignment_sha256'] != ASSIGNMENT_SHA256:
        raise ValueError('This example only replays the declared frozen experiment')
    groups = indexed(rows(payloads['groups'], compressed=True))
    topology = indexed(rows(payloads['topology']))
    eligibility = indexed(rows(payloads['eligibility']))
    if len(groups) != 42574 or set(eligibility) != set(groups) or set(topology) != set(groups):
        raise ValueError('Full-release input universes differ')
    metadata = {}
    for row in rows(payloads['original_membership']):
        selected = {k: row[k] for k in ('label', 'diversity_tier', 'diversity_stratum')}
        sid = row['structure_id']
        if sid in metadata and metadata[sid] != selected:
            raise ValueError('Inconsistent archived structure label or diversity stratum')
        metadata[sid] = selected
    statuses = rows(payloads['input_status'])
    status_keys = {(r['model'], r['structure_id']) for r in statuses}
    expected_keys = {(model, sid) for model in ('gbdt', 'matdeeplearn', 'moftransformer') for sid in metadata}
    if len(status_keys) != len(statuses) or status_keys != expected_keys:
        raise ValueError('Model-input status inventory is incomplete or duplicated')
    failed = {r['structure_id'] for r in statuses if r['status'] != 'SUCCESS'}
    rejected_groups = {groups[s]['effective_leakage_block'] for s in failed}
    excluded = {s for s in metadata if groups[s]['effective_leakage_block'] in rejected_groups}
    remaining = {s: r for s, r in metadata.items() if s not in excluded}
    if len(failed) != 14 or len(excluded) != 14 or Counter(r['label'] for r in remaining.values()) != {'CR': 4671, 'NCR': 1153}:
        raise ValueError('Archived common-input eligibility cannot be reproduced')
    full_group_labels = defaultdict(set)
    for sid, record in groups.items():
        full_group_labels[record['effective_leakage_block']].add(eligibility[sid]['label'])
    for sid, record in remaining.items():
        if full_group_labels[groups[sid]['effective_leakage_block']] != {record['label']}:
            raise ValueError('Selected structure belongs to a label-impure full-release group')
    fixed = set(json.loads(payloads['original_receipt'])['fixed_test_ids'])
    if len(fixed) != 468 or fixed & excluded:
        raise ValueError('The original 468-CR test is not preserved')
    assignments = reconstruct(remaining, groups, topology, fixed)
    actual_digest = hashlib.sha256(canonical(assignments)).hexdigest()
    if actual_digest != ASSIGNMENT_SHA256:
        raise ValueError('Current package does not reproduce the frozen assignment digest: ' + actual_digest)
    expected_rows = rows(payloads['expected_membership'])
    if len(expected_rows) != len(assignments) or any(
        {k: str(v) for k, v in row.items()} != expected for row, expected in zip(assignments, expected_rows)
    ):
        raise ValueError('Generated assignments differ from the archived membership table')
    run_counts = {}
    for key in sorted({row['run_key'] for row in assignments}):
        run = [row for row in assignments if row['run_key'] == key]
        by_group = defaultdict(set)
        for row in run:
            by_group[row['effective_leakage_block']].add(row['partition'])
        if any(len(parts) > 1 for parts in by_group.values()):
            raise ValueError('A related-structure group crosses partitions')
        run_counts[key] = dict(Counter(row['partition'] for row in run))
    report = dict(status='PASS_FROZEN_ASSIGNMENT_REPLAY', dataset_id=receipt['dataset_id'],
        assignment_sha256=actual_digest, assignment_rows=len(assignments), runs=run_counts,
        eligible_CR=4671, eligible_NCR=1153, common_test_structures=468,
        input_archive_sha256=expected_sha256, input_file_sha256=hashes, code_sha256=code_hashes,
        source='archived full-release groups, diversity strata and certified input-availability decisions',
        targets_or_predictions_read=False, science_or_training_run=False,
        grouping_or_diversity_recalculated=False, model_input_certification_rerun=False,
        full_release_label_purity_verified=True, crossed_groups=0,
        official_split=False, publication_authorized=False)
    if any(sha(p) != code_hashes[str(p)] for p in code_paths):
        raise ValueError('Replay code changed during execution')
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.coremof-assignment-replay-', dir=output.parent))
    try:
        with (staging / 'membership_manifest.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(assignments[0]))
            writer.writeheader()
            writer.writerows(assignments)
        (staging / 'receipt.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
        (staging / 'SHA256SUMS').write_text(''.join(sha(staging / name) + '  ' + name + '\n'
            for name in ('membership_manifest.csv', 'receipt.json')))
        publish_directory(staging, output, overwrite=False)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--handoff-archive', type=Path, required=True)
    parser.add_argument('--archive-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(replay(args.handoff_archive, args.archive_sha256, args.output), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
