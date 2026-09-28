"""Asset and temperature preflight, without loading scientific models."""

import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from CoREMOF import _heat_capacity as hc


class HeatCapacityAssetsTests(unittest.TestCase):
    def test_temperatures_normalize_without_rounding(self):
        self.assertEqual(hc.normalize_temperatures([300, 350.0]), [300, 350])
        for values in ([], [300, 300.0], [300.2], [True], [0], [-1],
                       [float("nan")], [float("inf")], ["300"]):
            with self.subTest(values=values), self.assertRaises(ValueError):
                hc.normalize_temperatures(values)

    def test_shipped_manifest_is_complete(self):
        path = Path(hc.__file__).parent / "models" / "heat_capacity_assets.json"
        data = json.loads(path.read_text())
        self.assertEqual(data["scikit_learn_version"], "1.4.2")
        self.assertEqual(data["xgboost_version"], "2.0.3")
        self.assertEqual(set(data["temperatures"]), {"300", "350", "400"})
        for assets in data["temperatures"].values():
            self.assertEqual(set(assets), {f"model_{i}" for i in range(100)})
            for record in assets.values():
                self.assertGreater(record["bytes"], 0)
                self.assertEqual(len(bytes.fromhex(record["sha256"])), 32)

    def test_ensemble_identity_and_runtime_checks(self):
        data = b"non-executable fixture, never unpickled"
        manifest = {"scikit_learn_version": "1.4.2", "xgboost_version": "2.0.3", "temperatures": {
            "300": {"model_0": {"bytes": len(data),
                                  "sha256": hashlib.sha256(data).hexdigest()}}}}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            folder = root / "300"
            folder.mkdir()
            asset = folder / "model_0"
            asset.write_bytes(data)
            with patch.object(hc.json, "loads", return_value=manifest), \
                 patch.object(hc.importlib.metadata, "version", side_effect={
                     "scikit-learn": "1.4.2", "xgboost": "2.0.3"}.__getitem__):
                hc.validate_ensemble(root, [300])
                with self.assertRaisesRegex(ValueError, "No bundled"):
                    hc.validate_ensemble(root, [301])
                (folder / "unexpected").touch()
                with self.assertRaisesRegex(ValueError, "unexpected"):
                    hc.validate_ensemble(root, [300])
                (folder / "unexpected").unlink()
                asset.write_bytes(b"corrupt")
                with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                    hc.validate_ensemble(root, [300])
                asset.unlink()
                with self.assertRaisesRegex(ValueError, "missing"):
                    hc.validate_ensemble(root, [300])
                elsewhere = root / "outside"
                elsewhere.write_bytes(data)
                asset.symlink_to(elsewhere)
                with self.assertRaisesRegex(ValueError, "regular file"):
                    hc.validate_ensemble(root, [300])
                asset.unlink()
                asset.write_bytes(data)
                with patch.object(hc.importlib.metadata, "version", side_effect=None, return_value="1.5.0"):
                    with self.assertRaisesRegex(RuntimeError, "separate environment"):
                        hc.validate_ensemble(root, [300])


if __name__ == "__main__":
    unittest.main()
