"""Create an additive source-only loader contract from authorized release data.

No calculation, CIF transformation, metadata replacement or publication occurs.
The selected-source tables must already be exact subsets of the complete release.
"""

import argparse
import json

from CoREMOF import export_source_projection
from CoREMOF.dataset import CoREMOFDataset


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--complete-release', required=True)
    parser.add_argument('--source-root', required=True)
    parser.add_argument('--output', required=True, help='new contract JSON path, must not exist')
    parser.add_argument('--sources', nargs='+', choices=('COD', 'CSD', 'SI'), required=True)
    parser.add_argument('--group-profile', action='append', help=(
        'ordered comma-separated criteria; repeat for multiple complete profiles; '
        'default: priority_main and RT,M2T'))
    parser.add_argument('--diversity', choices=('representative', 'none'), default='representative')
    parser.add_argument('--verify-cif-files', action='store_true', help='also verify complete-release CIF bytes')
    args = parser.parse_args()
    profiles = tuple(tuple(part.strip() for part in value.split(','))
                     for value in (args.group_profile or ('priority_main', 'RT,M2T')))
    dataset = CoREMOFDataset.from_release(args.complete_release, verify_cif_files=args.verify_cif_files)
    receipt = export_source_projection(dataset, args.source_root, args.output,
        sources=tuple(args.sources), group_profiles=profiles, diversity=args.diversity)
    print(json.dumps(dict(receipt), indent=2))


if __name__ == '__main__':
    main()
