#!/usr/bin/env python3
"""Explicitly attach an authorized optional target package; no calculations."""
import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('core_root', type=Path)
    parser.add_argument('target_root', type=Path)
    parser.add_argument('--manifest-sha256', required=True,
                        help='independently received target manifest SHA-256')
    args = parser.parse_args()
    from CoREMOF.dataset import CoREMOFDataset
    core = CoREMOFDataset.from_release(args.core_root, verify_cif_files=False)
    attached = core.attach_target_supplement(
        args.target_root, expected_sha256=args.manifest_sha256,
    )
    coverage = {name: sum(attached.target_values(sid)[name] is not None
                          for sid in attached.structure_ids)
                for name in attached.target_columns}
    print(json.dumps({'structures': len(core), 'target_coverage': coverage,
                      'explicit_attachment': True, 'cif_bytes_verified': False,
                      'assignments_changed': False, 'new_calculations': False,
                      'publication_authorization_granted': False}, indent=2))


if __name__ == '__main__':
    main()
