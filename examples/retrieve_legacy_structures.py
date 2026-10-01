"""Explicit historical downloads, separate from version-selected CoRE-MOF-COD.

Access rights and a working licensed CSD installation are caller requirements.
Use fetch_release.py with a pinned catalog for an exact database release.
"""
import argparse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', choices=('si', 'csd'))
    parser.add_argument('output', help='Explicit destination directory')
    parser.add_argument('--refcode', help='Required for the csd source')
    parser.add_argument('--overwrite', action='store_true', help='Intentionally replace selected files only')
    args = parser.parse_args()
    if args.source == 'csd' and not args.refcode:
        parser.error('csd requires --refcode')
    if args.source == 'si' and args.refcode:
        parser.error('--refcode is only supported for csd')
    from CoREMOF.structure import download_from_CSD, download_from_SI
    if args.source == 'csd':
        download_from_CSD(args.refcode, args.output, overwrite=args.overwrite)
    else:
        download_from_SI(args.output, overwrite=args.overwrite)


if __name__ == '__main__':
    main()
