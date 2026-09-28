"""Look up a DOI publication date without assigning or changing structure IDs."""
import argparse

from CoREMOF.calculation.get_info import get_publication_date


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('doi', help='DOI identifier without a resolver URL')
    options = parser.parse_args(argv)
    print(get_publication_date(options.doi))


if __name__ == '__main__':
    main()
