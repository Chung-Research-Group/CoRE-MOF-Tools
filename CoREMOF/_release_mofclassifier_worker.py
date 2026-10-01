"""Fresh-process launcher for the recorded MOFClassifier ensemble."""
import importlib.metadata
import importlib.util
import json
import platform
from pathlib import Path
import resource
import socket
import sys


PYTHON_SHA256 = '9024c0445314fe47972ac371dd460c4949aa9bb211ba8934f682c608bd209176'
PROTOCOL_SHA256 = 'bf7d571929f6bb9d883a518cb82e069a4b26cea671b1e8b8f5a8028a31a43fc0'
IDENTIFIERS_SHA256 = 'ce2642db0a903b701e8a770b570ac8fb4b47a182c2e53a3e6b02ce82493f07c3'
CONFIG_SHA256 = 'd341fecb4fc275bec143ce2a54cf88953ace1697122f44161e852de6b706fe64'
VERSIONS = {'MOFClassifier': '0.1.1', 'torch': '2.7.0+cu118', 'numpy': '1.26.4',
            'ase': '3.23.0', 'pymatgen': '2024.8.9'}


def deny_network(*args, **kwargs):
    raise RuntimeError('Recorded model replay does not download software or models')


def load_protocol(protocol_path):
    """Load the exact validator beside the protocol in an isolated process."""
    import hashlib
    protocol_path = Path(protocol_path)
    loaded = None
    for name, path, expected in (
        ('identifiers', protocol_path.with_name('identifiers.py'), IDENTIFIERS_SHA256),
        ('recorded_mofclassifier', protocol_path, PROTOCOL_SHA256),
    ):
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise RuntimeError('Recorded model method differs: ' + name)
        spec = importlib.util.spec_from_file_location(name, path)
        loaded = importlib.util.module_from_spec(spec)
        sys.modules[name] = loaded
        spec.loader.exec_module(loaded)
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise RuntimeError('Recorded model method changed during import: ' + name)
    return loaded


def main():
    import hashlib
    request = json.loads(Path(sys.argv[1]).read_text())
    limit = request['memory_limit_mb'] * 1024 * 1024
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    socket.socket.connect = deny_network
    socket.create_connection = deny_network
    if (sys.version_info[:3] != (3, 9, 23)
            or hashlib.sha256(Path(sys.executable).read_bytes()).hexdigest() != PYTHON_SHA256):
        raise RuntimeError('Recorded model Python runtime differs')
    versions = {key: importlib.metadata.version(key) for key in VERSIONS}
    if versions != VERSIONS:
        raise RuntimeError('Recorded model package versions differ: ' + repr(versions))
    for name, expected in (('protocol', PROTOCOL_SHA256), ('config', CONFIG_SHA256)):
        if hashlib.sha256(Path(request[name]).read_bytes()).hexdigest() != expected:
            raise RuntimeError('Recorded model method differs: ' + name)
    protocol = load_protocol(request['protocol'])
    summary = protocol.run(Path(request['manifest']), request['manifest_sha256'], Path(request['config']),
                           Path(request['output']), [0], package_root_override=Path(request['model_root']),
                           private_root=Path(request['private_root']))
    # Revalidate all 100 model assets after inference as well as at load time.
    contract = protocol.load_method_contract(Path(request['config']), Path(request['model_root']))
    import torch
    evidence = {'python_sha256': PYTHON_SHA256, 'versions': versions,
                'machine': platform.machine(), 'processor': platform.processor(),
                'torch_cpu_capability': torch.backends.cpu.get_cpu_capability(),
                'torch_build_configuration': torch.__config__.show(),
                'asset_bundle_sha256': contract.asset_bundle_sha256,
                'assets': [item.fingerprint_value() for item in
                           (contract.source_module, contract.atom_initializer, *contract.checkpoints)],
                'historical_full_runtime_byte_identity_proven': False,
                'summary': summary}
    Path(request['runtime_receipt']).write_text(json.dumps(evidence, indent=2, sort_keys=True) + '\n')


if __name__ == '__main__':
    main()
