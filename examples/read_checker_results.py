"""Read precomputed checker votes and derive CR/NCR without checker execution."""
import argparse
import json

from CoREMOF.dataset import CoREMOFDataset
from CoREMOF.labels import CHECKER_COLUMNS


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('release_root')
    parser.add_argument('--structure-id')
    parser.add_argument('--checkers', nargs='+', default=['MOFClassifier', 'MOFChecker', 'Chen-Manz', 'MOSAEC', 'SETC-GAT'])
    args = parser.parse_args(argv)
    dataset = CoREMOFDataset.from_release(args.release_root)
    view = dataset.classify(checkers=args.checkers)
    result = {'checker_view': view.checker_view, 'label_counts': dict(view.label_counts())}
    if args.structure_id:
        record = dataset[args.structure_id]
        result['structure'] = {'structure_id': record.structure_id,
                               'checkers': {name: record[column] for name, column in CHECKER_COLUMNS.items()},
                               'label': view[args.structure_id].label}
    print(json.dumps(result, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
