#!/usr/bin/env python3
"""Read a pinned metadata release or authenticated source projection, no science."""
import argparse
import hashlib
import json
import re
from pathlib import Path, PurePosixPath


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(chunk)
    return result.hexdigest()


def verify_metadata_ledger(root, expected):
    """A trusted ledger hash pins all its listed files; it is not a licence."""
    if not re.fullmatch(r'[0-9a-f]{64}', expected):
        raise ValueError('Expected a complete lowercase SHA-256')
    ledger = root / 'SHA256SUMS'
    if ledger.is_symlink() or digest(ledger) != expected:
        raise ValueError('Metadata ledger differs from independently received hash')
    seen = set()
    for line in ledger.read_text().splitlines():
        sha, name = line.split('  ', 1)
        rel = PurePosixPath(name)
        if (not re.fullmatch(r'[0-9a-f]{64}', sha) or rel.is_absolute() or
                '..' in rel.parts or '\\' in name or name != rel.as_posix() or
                name in seen):
            raise ValueError('Invalid metadata ledger row')
        seen.add(name)
        path = root.joinpath(*rel.parts)
        if (any(root.joinpath(*rel.parts[:i]).is_symlink()
                for i in range(1, len(rel.parts) + 1)) or
                not path.resolve().is_relative_to(root.resolve()) or
                not path.is_file() or digest(path) != sha):
            raise ValueError('Metadata checksum mismatch: ' + name)
    return len(seen)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path, help='authorized metadata/source root')
    parser.add_argument('--expected-version', required=True,
                        help='explicit dataset name, e.g. CoRE-MOF-COD')
    parser.add_argument('--expected-ledger-sha256', help='trusted metadata ledger hash')
    parser.add_argument('--projection-contract', type=Path)
    parser.add_argument('--projection-sha256', help='independently received contract hash')
    parser.add_argument('--checkers', nargs='+', default=['5checker'])
    args = parser.parse_args()
    if bool(args.projection_contract) != bool(args.projection_sha256):
        parser.error('source projection requires both contract and trusted SHA-256')
    if args.projection_contract and args.expected_ledger_sha256:
        parser.error('projection uses its own contract, not a complete-release ledger')
    root = args.root.resolve(strict=True)
    verified_files = None
    if args.expected_ledger_sha256:
        verified_files = verify_metadata_ledger(root, args.expected_ledger_sha256)
    # Delay import so --help works without installing the scientific package.
    from CoREMOF.dataset import CoREMOFDataset
    if args.projection_contract:
        dataset = CoREMOFDataset.from_projection(
            root, args.projection_contract, expected_sha256=args.projection_sha256,
            verify_cif_files=False,
        )
    else:
        dataset = CoREMOFDataset.from_release(root, verify_cif_files=False)
    if dataset.dataset_version != args.expected_version:
        raise ValueError('Loaded version differs from explicitly requested version')
    checkers = args.checkers[0] if len(args.checkers) == 1 else tuple(args.checkers)
    view = dataset.classify(checkers=checkers)
    print(json.dumps({'dataset_version': dataset.dataset_version,
                      'structures': len(dataset), 'checkers': args.checkers,
                      'label_counts': dict(view.label_counts()),
                      'metadata_files_verified_against_trusted_ledger': verified_files,
                      'cif_bytes_verified': False, 'new_calculations': False,
                      'benchmark_assignments_changed': False,
                      'release_status': dataset.dataset_info.get('release_status'),
                      'publication_authorization_granted_by_example': False}, indent=2))


if __name__ == '__main__':
    main()
