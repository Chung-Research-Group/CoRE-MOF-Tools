#!/usr/bin/env python3
"""Read-only pre-upload check of Git-visible files, not permission to publish data."""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def candidate_paths(root):
    # Check both a fresh add and the existing index: .gitignore cannot remove
    # data or retired checker files that have already been tracked.
    with tempfile.TemporaryDirectory(prefix='coremof-upload-check-') as tmp:
        git = Path(tmp) / 'index.git'
        subprocess.run(['git', 'init', '--bare', '--quiet', str(git)], check=True)
        command = ['git', '--git-dir=' + str(git), '--work-tree=' + str(root),
                   'ls-files', '--others', '--exclude-standard', '-z']
        result = subprocess.run(command, cwd=root, check=True, stdout=subprocess.PIPE)
    names = {os.fsdecode(x) for x in result.stdout.split(b'\0') if x}
    inside = subprocess.run(['git', '-C', str(root), 'rev-parse', '--is-inside-work-tree'],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if inside.returncode == 0 and inside.stdout.strip() == b'true':
        tracked = subprocess.run(['git', '-C', str(root), 'ls-files', '--cached', '-z', '--', '.'],
                                 check=True, stdout=subprocess.PIPE)
        names.update(os.fsdecode(x) for x in tracked.stdout.split(b'\0') if x)
    return sorted(names)


def inspect(root):
    errors, records = [], []
    forbidden_workers = {'_release_checkers_protocol.py', '_release_checkers_worker.py',
                         '_release_mosaec_worker.py', '_release_setc_protocol.py',
                         '_release_setc_worker.py'}
    secret_pattern = re.compile(rb'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----|'
                                rb'gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{50,}')
    for name in candidate_paths(root):
        path = root / name
        rel = Path(name)
        if name in {'GIT_UPLOAD_MANIFEST.json', 'CODE_SHA256SUMS'}:
            continue
        if path.is_symlink():
            errors.append(name + ': symlink is not permitted')
            continue
        if (not path.is_file() or rel.is_absolute() or '..' in rel.parts):
            errors.append(name + ': invalid candidate')
            continue
        size = path.stat().st_size
        if size > 20 * 1024 * 1024:
            errors.append(name + ': exceeds this handoff\'s conservative 20 MiB code-file limit')
            continue
        if ('local' in rel.parts or path.suffix.lower() in
                {'.cif', '.h5', '.hdf5', '.pkl', '.pickle', '.pt', '.pth', '.ckpt', '.zip', '.gz', '.whl'}
                or name in {'CoREMOF/data/CR.json', 'CoREMOF/data/NCR.json'}
                or name.startswith(('CoREMOF/models/stability/', 'CoREMOF/models/cp_app/models/',
                                    'CoREMOF/models/cp_app/ensemble_models_'))):
            errors.append(name + ': local data/model/archive visible to Git')
        if ('mosaec' in rel.parts and 'data' in rel.parts) or path.name in forbidden_workers:
            errors.append(name + ': removed checker implementation visible to Git')
        if '/data/' in name and (name.startswith('2026-CoREMOF-COD/main/') or
                                name.startswith('2026-CoREMOF-COD/si/')):
            errors.append(name + ': private plot input visible to Git')
        data = path.read_bytes()
        if secret_pattern.search(data):
            errors.append(name + ': possible credential/private key')
        if path.suffix == '.py':
            try:
                ast.parse(data.decode('utf-8'), filename=name)
            except (SyntaxError, UnicodeError) as exc:
                errors.append(name + ': ' + str(exc))
        if path.suffix == '.ipynb':
            try:
                notebook = json.loads(data)
                for cell in notebook['cells']:
                    if cell.get('cell_type') == 'code' and (
                            cell.get('outputs') or cell.get('execution_count') is not None):
                        errors.append(name + ': notebook output/execution count must be cleared')
                        break
            except (ValueError, KeyError) as exc:
                errors.append(name + ': invalid notebook: ' + str(exc))
        records.append({'path': name, 'bytes': size, 'sha256': hashlib.sha256(data).hexdigest()})
    return {'status': 'FAIL' if errors else 'PASS_CODE_SURFACE', 'files': records,
            'file_count': len(records), 'bytes': sum(x['bytes'] for x in records),
            'errors': errors, 'publication_authorized': False,
            'scope': 'Local code surface only. Not a software or data licence decision.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--write-manifest', action='store_true',
                        help='write/refresh GIT_UPLOAD_MANIFEST.json with current code hashes')
    args = parser.parse_args()
    root = args.root.resolve()
    result = inspect(root)
    if not args.write_manifest and (root / 'GIT_UPLOAD_MANIFEST.json').is_file():
        try:
            manifest = json.loads((root / 'GIT_UPLOAD_MANIFEST.json').read_text())
            if manifest['files'] != result['files']:
                result['errors'].append('Git-visible files differ from GIT_UPLOAD_MANIFEST.json; '
                                        'review the changes before refreshing the manifest')
        except (ValueError, KeyError, TypeError):
            result['errors'].append('GIT_UPLOAD_MANIFEST.json is invalid')
        if result['errors']:
            result['status'] = 'FAIL'
    if args.write_manifest and not result['errors']:
        (root / 'GIT_UPLOAD_MANIFEST.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'files'}, indent=2))
    return int(bool(result['errors']))


if __name__ == '__main__':
    raise SystemExit(main())

