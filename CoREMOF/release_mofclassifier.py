"""Reproduce the recorded 100-bag MOFClassifier core CPU scores."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

from . import _release_mofclassifier_protocol as protocol
from . import _release_mofclassifier_worker as worker
from ._transactions import publish_directory
from .release_curation import _json, _stop
from ._release_curation_worker import check, sha


PROFILE = 'mofclassifier-0.1.1-core-100bag-cpu-private-copy-v1'
ASSET_BUNDLE_SHA256 = '373348ca0252817016db407931014f7b1887b40d995bf15e8c39a253abda8e2f'


class ReleaseMOFClassifierError(RuntimeError):
    """Recorded model evidence is incomplete or inconsistent."""


def _validate(record, structure_id, cif_sha):
    if (record.get('structure_id') != structure_id or record.get('canonical_cif_sha256') != cif_sha
            or record.get('method_id') != PROFILE or record.get('asset_bundle_sha256') != ASSET_BUNDLE_SHA256
            or record.get('method_config_sha256') != worker.CONFIG_SHA256
            or record.get('canonical_post_run_sha256') != cif_sha or record.get('canonical_unchanged') is not True):
        raise ReleaseMOFClassifierError('Model record identity, assets or source-CIF binding differs')
    if (record.get('model_family') != 'core' or record.get('checkpoint_bag_count') != 100
            or record.get('pass_threshold') != 0.6 or record.get('positive_class_index') != 1
            or record.get('mean_aggregation') != 'float64_math_fsum_divide_by_100'):
        raise ReleaseMOFClassifierError('Unexpected recorded model or score convention')
    scores = record.get('bag_scores')
    count = record.get('completed_bag_count')
    if (not isinstance(scores, list) or type(count) is not int or count != len(scores) or not 0 <= count <= 100
            or any(type(x) not in (int, float) or not math.isfinite(x) or not 0 <= x <= 1 for x in scores)):
        raise ReleaseMOFClassifierError('Invalid or inconsistent individual model scores')
    if record.get('execution_status') == 'SUCCESS':
        if count != 100 or record.get('private_copy_used') is not True:
            raise ReleaseMOFClassifierError('Successful inference lacks all models or private input')
        mean = math.fsum(scores) / 100
        if (record.get('mean_score') != mean or record.get('operational_vote') != ('PASS' if mean >= 0.6 else 'FAIL')
                or record.get('private_copy_initial_sha256') != cif_sha
                or not re.fullmatch('[0-9a-f]{64}', str(record.get('private_copy_model_input_sha256')))):
            raise ReleaseMOFClassifierError('Mean, vote or private-copy evidence disagrees')
    elif record.get('execution_status') == 'ERROR':
        if record.get('mean_score') is not None or record.get('operational_vote') is not None or not record.get('error_type'):
            raise ReleaseMOFClassifierError('Unavailable inference must have null mean and vote')
    else:
        raise ReleaseMOFClassifierError('Nonterminal inference record')


def calculate_release_mofclassifier(cif_path, *, structure_id, output_dir, python,
                                    model_root, timeout_seconds=300, memory_limit_mb=8192):
    """Run all 100 pinned core models on CPU without altering the source CIF.

    The returned CL score uses the recorded binary64 ``math.fsum`` mean and
    PASS threshold >= 0.6. Missing/failed models yield unavailable, not FAIL.
    All weights, source and atom embeddings must already be installed. No
    downloads, training, replacement weights or release-label updates occur.
    The upstream parser fallback is retained only on its isolated private copy.
    """
    if os.name != 'posix':
        raise ReleaseMOFClassifierError('Recorded external runtime requires POSIX')
    from .identifiers import parse_core_id
    parse_core_id(structure_id)
    if type(timeout_seconds) is not int or timeout_seconds < 1 or type(memory_limit_mb) is not int or memory_limit_mb < 1024:
        raise ValueError('Use a positive integer timeout and memory_limit_mb >= 1024')
    target = Path(output_dir).absolute()
    if target.exists() or target.is_symlink():
        raise FileExistsError(target)
    if not target.parent.is_dir():
        raise FileNotFoundError(target.parent)
    python = check(Path(python).absolute(), worker.PYTHON_SHA256)
    runner = check(protocol.__file__, worker.PROTOCOL_SHA256)
    config = check(Path(__file__).with_name('data') / 'mofclassifier_fresh_core_v1.json', worker.CONFIG_SHA256)
    models = Path(model_root).resolve(strict=True)
    contract = protocol.load_method_contract(config, models)
    if contract.asset_bundle_sha256 != ASSET_BUNDLE_SHA256:
        raise ReleaseMOFClassifierError('Model asset bundle differs')
    source = Path(cif_path).resolve(strict=True)
    data = source.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    with tempfile.TemporaryDirectory(prefix='.mofclassifier-replay-', dir=target.parent) as temporary:
        root = Path(temporary)
        for name in ('input', 'output', 'work', 'home', 'tmp'):
            (root / name).mkdir(mode=0o700)
        output = root / 'output'
        cif = root / 'input' / (structure_id + '.cif')
        cif.write_bytes(data)
        row = dict(manifest_schema_version='1.0', row_index='0', persistent_uid='coremof:replay:' + structure_id,
                   structure_id=structure_id, canonical_cif_path=str(cif), canonical_cif_size_bytes=str(len(data)),
                   canonical_cif_sha256=digest)
        manifest = root / 'manifest.csv'
        with manifest.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=protocol.MANIFEST_FIELDS)
            writer.writeheader()
            writer.writerow(row)
        private_worker = root / 'worker.py'
        private_worker.write_bytes(Path(worker.__file__).read_bytes())
        # The protocol calculates its default path from its own location. This
        # explicit invocation always supplies the hash-bound configuration.
        request = dict(protocol=str(runner), config=str(config), model_root=str(models),
                       output=str(output), manifest=str(manifest), manifest_sha256=sha(manifest),
                       private_root=str(root / 'work'), memory_limit_mb=memory_limit_mb,
                       runtime_receipt=str(output / 'runtime.json'))
        _json(root / 'request.json', request)
        env = {'PATH': str(python.parent) + ':/usr/bin:/bin', 'HOME': str(root / 'home'), 'TMPDIR': str(root / 'tmp'),
               'MPLCONFIGDIR': str(root / 'tmp' / 'matplotlib'), 'PYTHONDONTWRITEBYTECODE': '1',
               'PYTHONNOUSERSITE': '1', 'PYTHONHASHSEED': '20260719', 'CUDA_VISIBLE_DEVICES': '',
               'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1',
               'LANG': 'C.UTF-8', 'LC_ALL': 'C.UTF-8', 'LD_LIBRARY_PATH': str(python.parent.parent / 'lib')}
        with (output / 'stdout.txt').open('w') as stdout, (output / 'stderr.txt').open('w') as stderr:
            process = subprocess.Popen([str(python), '-B', '-s', str(private_worker), str(root / 'request.json')],
                                       cwd=root / 'work', env=env, stdout=stdout, stderr=stderr, start_new_session=True)
            timed_out = False
            try:
                process.wait(timeout=timeout_seconds)
            except subprocess.TimeoutExpired:
                timed_out = True
                _stop(process)
            except BaseException:
                _stop(process)
                raise
        check(source, digest)
        check(cif, digest)
        check(runner, worker.PROTOCOL_SHA256)
        check(config, worker.CONFIG_SHA256)
        path = output / 'records' / ('00000000_' + structure_id + '.json')
        if process.returncode == 0 and path.is_file() and (output / 'runtime.json').is_file() and not timed_out:
            record = json.loads(path.read_text())
            _validate(record, structure_id, digest)
            result = {key: record[key] for key in ('execution_status', 'bag_scores', 'completed_bag_count', 'mean_score',
                       'operational_vote', 'error_type', 'error_message', 'private_copy_rewritten', 'private_copy_model_input_sha256')}
        else:
            result = {'execution_status': 'TIMEOUT' if timed_out else 'PROCESS_ERROR', 'mean_score': None,
                      'operational_vote': None, 'bag_scores': [], 'completed_bag_count': None,
                      'error_type': 'TIMEOUT' if timed_out else 'PROCESS_ERROR',
                      'error_message': f'Worker return code {process.returncode}; see stderr.txt',
                      'private_copy_rewritten': None, 'private_copy_model_input_sha256': None}
        result.update(schema_version='coremof-release-mofclassifier/1.0', profile=PROFILE,
                      structure_id=structure_id, cif_sha256=digest, pass_threshold=0.6,
                      source_cif_modified=False, release_labels_updated=False, mechanistic_hard_fail_eligible=False)
        receipt = dict(profile=PROFILE, protocol_sha256=worker.PROTOCOL_SHA256, config_sha256=worker.CONFIG_SHA256,
                       asset_bundle_sha256=ASSET_BUNDLE_SHA256, python_sha256=worker.PYTHON_SHA256,
                       worker_sha256=sha(private_worker), cif_sha256=digest, timeout_seconds=timeout_seconds,
                       memory_limit_mb=memory_limit_mb, returncode=process.returncode,
                       historical_full_runtime_byte_identity_proven=False, release_metadata_promoted=False)
        _json(output / 'record.json', result)
        _json(output / 'receipt.json', receipt)
        publish_directory(output, target, overwrite=False)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cif_path')
    for name in ('structure-id', 'output-dir', 'python', 'model-root'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--timeout-seconds', type=int, default=300)
    parser.add_argument('--memory-limit-mb', type=int, default=8192)
    result = calculate_release_mofclassifier(**vars(parser.parse_args(argv)))
    print(json.dumps({k: result[k] for k in ('structure_id', 'execution_status', 'mean_score', 'operational_vote')}))
    return 0 if result['execution_status'] == 'SUCCESS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
