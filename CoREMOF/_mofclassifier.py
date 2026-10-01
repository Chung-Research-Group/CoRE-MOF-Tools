"""Input-safe adapter for the optional upstream MOFClassifier batch API."""
import hashlib
import importlib
import importlib.util
import json
import math
from numbers import Real
from pathlib import Path
import tempfile
import warnings

from ._transactions import publish_file_bundle


def _load_classifier(model):
    spec = importlib.util.find_spec('MOFClassifier')
    if spec is None or not spec.submodule_search_locations:
        raise ImportError('Install MOFClassifier==0.1.1 and its model assets before inference')
    root = Path(next(iter(spec.submodule_search_locations)))
    # CLscore downloads assets at import time. Require them first, including
    # the other families whose directory checks run even for a core-only call.
    for name in ('models', 'models_qsp', 'models_h'):
        folder = root / name
        if not folder.is_dir() or not any(folder.iterdir()):
            raise FileNotFoundError('MOFClassifier assets must be installed explicitly: ' + str(folder))
    family = {'core': 'models', 'qsp': 'models_qsp', 'h': 'models_h'}[model]
    required = [root / 'atom_init.json'] + [root / family / f'checkpoint_bag_{i}.pth.tar' for i in range(1, 101)]
    for path in required:
        if not path.is_file() or not path.stat().st_size:
            raise FileNotFoundError('Incomplete MOFClassifier ensemble: ' + str(path))
    return importlib.import_module('MOFClassifier.CLscore')


def _score(value):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError('MOFClassifier scores must be finite numbers in [0, 1]')
    value = float(value)
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError('MOFClassifier scores must be finite numbers in [0, 1]')
    return value


def predict_directory(cif_folder, save_path, model, batch_size, *, overwrite=False):
    """Keep upstream batch/mean conventions, isolate its mutable CIF parser."""
    if model not in {'core', 'qsp', 'h'}:
        raise ValueError('model must be core, qsp or h')
    if type(batch_size) is not int or batch_size < 1:
        raise ValueError('batch_size must be a positive integer')
    if type(overwrite) is not bool:
        raise TypeError('overwrite must be a Boolean')
    folder = Path(cif_folder).resolve(strict=True)
    if not folder.is_dir():
        raise NotADirectoryError(folder)
    sources = sorted(path for path in folder.iterdir() if path.is_file() and path.suffix.lower() == '.cif')
    if not sources:
        raise FileNotFoundError('No CIF files were found in ' + str(folder))
    target = Path(save_path).absolute()
    if target.is_symlink() or target.is_dir() or (target.exists() and not overwrite):
        raise FileExistsError(target)
    if not target.parent.is_dir():
        raise FileNotFoundError(target.parent)
    if target.resolve() in {path.resolve() for path in sources}:
        raise ValueError('The result path must not replace an input CIF')
    expected = {path.stem for path in sources}
    if len(expected) != len(sources):
        raise ValueError('CIF filenames must have unique structure IDs')
    classifier = _load_classifier(model)
    originals = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}
    with tempfile.TemporaryDirectory(prefix='.mofclassifier-', dir=target.parent) as directory:
        root = Path(directory)
        inputs = root / 'inputs'
        inputs.mkdir()
        copies = []
        for source in sources:
            data = source.read_bytes()
            if hashlib.sha256(data).hexdigest() != originals[source]:
                raise RuntimeError('Input CIF changed before inference: ' + str(source))
            private = inputs / source.name
            private.write_bytes(data)
            copies.append(str(private))
        try:
            results = classifier.predict_batch(root_cifs=copies, model=model, batch_size=batch_size)
        finally:
            for source, digest in originals.items():
                if hashlib.sha256(source.read_bytes()).hexdigest() != digest:
                    raise RuntimeError('Input CIF changed during inference: ' + str(source))
        out = {}
        for rid, scores, mean in results:
            if rid not in expected or rid in out:
                raise ValueError('Unexpected or duplicate MOFClassifier structure ID: ' + str(rid))
            scores = [_score(value) for value in scores]
            mean = _score(mean)
            if len(scores) != 100 or not math.isclose(mean, math.fsum(scores) / 100, rel_tol=0, abs_tol=1e-6):
                raise ValueError('MOFClassifier did not return a complete, consistent 100-model ensemble')
            out[rid] = [scores, mean]
        if set(out) != expected:
            raise ValueError('MOFClassifier did not return every input structure')
        rewritten = [source.stem for source in sources
                     if hashlib.sha256((inputs / source.name).read_bytes()).hexdigest() != originals[source]]
        if rewritten:
            warnings.warn('MOFClassifier rewrote private parser inputs, not source CIFs: ' + ', '.join(rewritten),
                          RuntimeWarning, stacklevel=2)
        staged = root / 'results.json'
        staged.write_text(json.dumps(out, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
        publish_file_bundle([staged], [target], overwrite=overwrite)
    return out
