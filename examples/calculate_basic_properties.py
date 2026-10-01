"""Read cell mass/volume, atom count and symmetry without changing a CIF.

Requires the ASE/pymatgen dependencies of the historical scientific API.
These functions describe the supplied representation, not a repaired crystal.
"""
import argparse
import hashlib
import json
from pathlib import Path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cif', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    from CoREMOF.calculation.mof_features import Mass, Volume, n_atom, SpaceGroup

    before = hashlib.sha256(args.cif.read_bytes()).hexdigest()
    result = {'cif_sha256': before, 'mass': Mass(args.cif), 'volume': Volume(args.cif),
              'atoms': n_atom(args.cif), 'space_group': SpaceGroup(args.cif)}
    if hashlib.sha256(args.cif.read_bytes()).hexdigest() != before:
        raise RuntimeError('Input CIF changed during calculation')
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')


if __name__ == '__main__':
    main()
