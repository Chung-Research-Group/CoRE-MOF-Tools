#!/usr/bin/env python3
"""Private single-CIF COD replay. No download, checker scan or release mutation.

Example::

    python replay_cod_curation.py source.cif --cod-id 7135365 \
        --python /path/to/recorded/python --workflow-root /path/to/COD_new \
        --output-dir new_cod_result

Use --skip-charges for inspection-only proposals. Such results remain REVIEW.
See docs/source/cod_curation_replay.rst for dependencies and numerical limits.
"""
from CoREMOF.cod_curation import main


if __name__ == '__main__':
    raise SystemExit(main())
