"""Replay the recorded Zeo++ open-metal-site detector, not the legacy detector."""
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
import sys
import tempfile

from . import _release_zeopp_oms_protocol as protocol
from .release_zeopp import NETWORK_SHA256, _check, _sha, _stop
from ._transactions import publish_directory


PROTOCOL_SHA256 = "d1e5f9c0f2c8244ad92045d8e85609ca1a964e80656b09e9f6785cb95729d533"
CONTRACT_SHA256 = "f6161a909d1f751a67d4fef310ae8e74cdfade6b68016522ee71df12ce31f648"


class ReleaseOMSError(RuntimeError):
    """A recorded-method OMS input or output is inconsistent."""


def _validate(record, structure_id, cif_sha):
    if (record.get('structure_id') != structure_id or record.get('protocol_id') != protocol.PROTOCOL_ID
            or record.get('input', {}).get('cif_sha256') != cif_sha
            or record.get('execution_status') not in {'SUCCESS', 'ERROR'}):
        raise ReleaseOMSError('Unexpected OMS result identity, method or status')
    props = record.get('open_metal_site_props')
    if record['execution_status'] == 'ERROR':
        if props is not None:
            raise ReleaseOMSError('Unavailable OMS calculation contains a value')
        return
    if not isinstance(props, dict):
        raise ReleaseOMSError('Successful OMS result has no properties')
    parsed = protocol.parse_oms(props.get('raw_line', ''))
    count = props.get('open_metal_site_count')
    if (type(count) is not int or count < 0 or count != parsed['open_metal_site_count']
            or props.get('has_open_metal_sites') is not (count > 0)):
        raise ReleaseOMSError('OMS count/boolean disagrees with raw output')
    surface = props.get('surface_definition_A')
    if surface is not None and (isinstance(surface, bool) or not isinstance(surface, (float, int))
                                or not math.isfinite(surface) or surface < 0):
        raise ReleaseOMSError('Invalid surface-definition distance')


def calculate_release_oms(cif_path, *, structure_id, output_dir, network, timeout_seconds=300):
    """Calculate the recorded ``network -oms CIF`` result on an isolated copy.

    A successful zero is retained. An execution failure is unavailable, not a
    zero and not NCR. The detector is not experimental evidence of accessible
    or catalytically active sites. No probes, high-accuracy switch or alternative
    detector are substituted for the exact recorded command.
    """
    if os.name != 'posix':
        raise ReleaseOMSError('The recorded external binary requires POSIX')
    if not isinstance(structure_id, str) or re.fullmatch(
            r'(?:ASR|FSR|ION)-(?:COD|CSD|SI)-(?:[0-9]{4}|UNKN)-[0-9]{4,}', structure_id) is None:
        raise ValueError('Use a public CoRE-MOF structure ID')
    if type(timeout_seconds) is not int or timeout_seconds <= 0:
        raise ValueError('timeout_seconds must be a positive integer')
    destination = Path(output_dir).absolute()
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    if not destination.parent.is_dir():
        raise FileNotFoundError(destination.parent)
    network = _check(Path(network).resolve(strict=True), NETWORK_SHA256)
    runner = _check(protocol.__file__, PROTOCOL_SHA256)
    contract = _check(Path(__file__).with_name('data') / 'zeopp_oms_metadata_contract_v1.json', CONTRACT_SHA256)
    source = Path(cif_path).resolve(strict=True)
    data = source.read_bytes()
    cif_sha = hashlib.sha256(data).hexdigest()
    with tempfile.TemporaryDirectory(prefix='.zeopp-oms-replay-', dir=destination.parent) as temporary:
        root = Path(temporary)
        for name in ('input', 'home', 'tmp', 'work', 'output'):
            (root / name).mkdir(mode=0o700)
        output = root / 'output'
        cif = root / 'input' / (structure_id + '.cif')
        cif.write_bytes(data)
        row = dict.fromkeys(protocol.MANIFEST_FIELDS, '')
        row.update(manifest_schema_version='1.0', row_index='0', structure_id=structure_id,
                   canonical_cif_version='user-supplied-exact-bytes', canonical_cif_path=str(cif),
                   cif_basename=cif.name, cif_size_bytes=str(len(data)), cif_sha256=cif_sha,
                   source_family=structure_id.split('-')[1])
        manifest = root / 'manifest.csv'
        with manifest.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=protocol.MANIFEST_FIELDS)
            writer.writeheader()
            writer.writerow(row)
        command = [sys.executable, '-B', '-S', str(runner), '--manifest', str(manifest),
                   '--manifest-sha256', _sha(manifest), '--row-index', '0', '--output-root', str(output),
                   '--private-root', str(root / 'work'), '--network', str(network),
                   '--network-sha256', NETWORK_SHA256, '--runner-sha256', PROTOCOL_SHA256,
                   '--contract', str(contract), '--contract-sha256', CONTRACT_SHA256,
                   '--timeout-seconds', str(timeout_seconds)]
        env = {'PATH': '/usr/bin:/bin', 'HOME': str(root / 'home'), 'TMPDIR': str(root / 'tmp'),
               'LANG': 'C.UTF-8', 'LC_ALL': 'C.UTF-8', 'OMP_NUM_THREADS': '1',
               'OPENBLAS_NUM_THREADS': '1', 'PYTHONDONTWRITEBYTECODE': '1'}
        with (output / 'stdout.txt').open('w') as stdout, (output / 'stderr.txt').open('w') as stderr:
            process = subprocess.Popen(command, cwd=root / 'work', env=env,
                                       stdout=stdout, stderr=stderr, start_new_session=True)
            timed_out = False
            try:
                process.wait(timeout=timeout_seconds + 30)
            except subprocess.TimeoutExpired:
                timed_out = True
                _stop(process)
            except BaseException:
                _stop(process)
                raise
        record_path = output / 'records' / (structure_id + '.json')
        if process.returncode == 0 and record_path.is_file() and not timed_out:
            record = json.loads(record_path.read_text())
        else:
            record = {'structure_id': structure_id, 'protocol_id': protocol.PROTOCOL_ID,
                      'input': {'cif_sha256': cif_sha}, 'execution_status': 'ERROR',
                      'open_metal_site_props': None, 'error_type': 'TIMEOUT' if timed_out else 'PROCESS_ERROR',
                      'error_message': f'Worker return code {process.returncode}; see stderr.txt'}
            record_path.parent.mkdir(exist_ok=True)
            record_path.write_text(json.dumps(record, indent=2) + '\n')
        _validate(record, structure_id, cif_sha)
        for path, digest in ((source, cif_sha), (cif, cif_sha), (network, NETWORK_SHA256),
                             (runner, PROTOCOL_SHA256), (contract, CONTRACT_SHA256)):
            _check(path, digest)
        props = record.get('open_metal_site_props')
        result = {'schema_version': 'coremof-release-oms-replay/1.0', 'profile': protocol.PROTOCOL_ID,
                  'structure_id': structure_id, 'cif_sha256': cif_sha,
                  'execution_status': record['execution_status'],
                  'open_metal_site_props': {k: v for k, v in props.items() if k != 'raw_line'} if props is not None else None,
                  'error_type': record.get('error_type'), 'error_message': record.get('error_message'),
                  'automatic_cr_ncr_exclusion_authorized': False}
        receipt = {'protocol_sha256': PROTOCOL_SHA256, 'contract_sha256': CONTRACT_SHA256,
                   'network_sha256': NETWORK_SHA256, 'cif_sha256': cif_sha,
                   'timeout_seconds': timeout_seconds, 'returncode': process.returncode,
                   'source_cif_modified': False, 'release_metadata_promoted': False,
                   'historical_full_runtime_byte_identity_proven': False,
                   'raw_record_sha256': _sha(record_path)}
        for name, value in (('record.json', result), ('receipt.json', receipt)):
            (output / name).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
        publish_directory(output, destination, overwrite=False)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cif_path')
    for name in ('structure-id', 'output-dir', 'network'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--timeout-seconds', type=int, default=300)
    result = calculate_release_oms(**vars(parser.parse_args(argv)))
    print(json.dumps({key: result[key] for key in ('structure_id', 'execution_status')}))
    return 0 if result['execution_status'] == 'SUCCESS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
