"""Run the recorded 100-model CPU ensemble using already-installed assets.

Example:
    python replay_release_mofclassifier.py structure.cif --structure-id "2016[Co][sqc27]3[FSR]3" \
        --output-dir replay --python /path/to/recorded/env/bin/python \
        --model-root /path/to/recorded/env/lib/python3.9/site-packages/MOFClassifier
"""
from CoREMOF.release_mofclassifier import main


if __name__ == '__main__':
    raise SystemExit(main())
