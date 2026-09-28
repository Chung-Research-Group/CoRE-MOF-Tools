"""Historical predictor guards, independent of TensorFlow and real weights."""
import ast
import hashlib
import json
import sys
from pathlib import Path
import tempfile
import unittest
import types
from unittest.mock import patch

from CoREMOF import prediction
from CoREMOF import _historical_stability as contract


class HistoricalStabilityTests(unittest.TestCase):
    def test_metrics_accept_the_keras3_ops_interface(self):
        epsilon = 1e-7
        operations = types.SimpleNamespace(sum=lambda x:x, round=round,
                                           clip=lambda x,a,b:min(max(x,a),b))
        backend = types.ModuleType('keras.backend')
        backend.epsilon = lambda:epsilon
        keras = types.ModuleType('keras')
        keras.ops, keras.backend = operations, backend
        with patch.dict(sys.modules, {'keras':keras,'keras.backend':backend}):
            expected = 1/(1+epsilon)
            self.assertEqual(prediction.precision(1.,1.),expected)
            self.assertEqual(prediction.recall(1.,1.),expected)
            self.assertEqual(prediction.f1(1.,1.),2*expected*expected/(2*expected+epsilon))

    def test_metrics_keep_the_legacy_backend_interface(self):
        backend = types.ModuleType('keras.backend')
        backend.sum = lambda x:x
        backend.round = round
        backend.clip = lambda x,a,b:min(max(x,a),b)
        backend.epsilon = lambda:1e-7
        keras = types.ModuleType('keras')
        keras.backend = backend
        with patch.dict(sys.modules, {'keras':keras,'keras.backend':backend}):
            self.assertEqual(prediction.precision(0.,0.),0.)
            self.assertEqual(prediction.recall(0.,0.),0.)
            self.assertEqual(prediction.f1(0.,0.),0.)

    def test_seven_known_assets_are_bound(self):
        self.assertEqual(len(contract.ASSET_SHA256), 7)
        self.assertEqual(sum(n.endswith('.h5') for n in contract.ASSET_SHA256), 2)
        self.assertTrue(all(len(h) == 64 for h in contract.ASSET_SHA256.values()))

    def test_original_feature_lists_and_order_are_preserved(self):
        tree = ast.parse(Path(prediction.__file__).read_text())
        expected = {
            'solvent_feature_names': (134, '5644d90b8ff5ebc59e5c4525002915c87acd7e766b39cc0383ac16ac80fbee76'),
            'thermal_feature_names': (134, '5644d90b8ff5ebc59e5c4525002915c87acd7e766b39cc0383ac16ac80fbee76'),
            'water_feature_names': (12, '5b0fa887332bfa370fffcb0b68d5bec4e16af7b22fd077e8c59cabfa95948f88'),
        }
        found = {}
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name):
                name = node.targets[0].id
                if name in expected:
                    values = ast.literal_eval(node.value)
                    found[name] = (len(values), hashlib.sha256(json.dumps(values, separators=(',', ':')).encode()).hexdigest())
        self.assertEqual(found, expected)

    def test_oversized_asset_read_is_bounded(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root/'fake.pkl').write_bytes(b'x'*17)
            with patch.dict(contract.ASSET_SHA256, {'fake.pkl':'0'*64}, clear=True), \
                 patch.object(contract, 'MAX_ASSET_BYTES', 16):
                with self.assertRaisesRegex(ValueError, 'size bound'):
                    contract.copy_verified_models(root, root/'output')

    def test_missing_models_fail_before_scientific_import_or_inference(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cif = root/'input.cif'
            cif.write_bytes(b'fixture')
            with patch.object(prediction, '_stability_prediction') as infer:
                with self.assertRaisesRegex(FileNotFoundError, 'Historical stability models are missing'):
                    prediction.stability(cif, model_directory=root/'missing')
                infer.assert_not_called()

    def test_mismatched_models_are_not_deserialized(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            models = root/'models'
            models.mkdir()
            for name in contract.ASSET_SHA256:
                (models/name).write_bytes(b'untrusted')
            cif = root/'input.cif'
            cif.write_bytes(b'fixture')
            with patch.object(prediction, '_stability_prediction') as infer:
                with self.assertRaisesRegex(ValueError, 'SHA-256 mismatch'):
                    prediction.stability(cif, model_directory=models)
                infer.assert_not_called()

    def test_verified_assets_are_isolated_before_deserialization(self):
        payload = b'synthetic model, never deserialized'
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            models = root/'models'
            models.mkdir()
            (models/'fake.pkl').write_bytes(payload)
            cif = root/'original.CIF'
            cif.write_bytes(b'unchanged CIF fixture')
            observed = []

            def infer(structure, model_directory):
                self.assertIsInstance(structure, str)
                private = Path(structure)
                self.assertNotEqual(private, cif)
                self.assertEqual(private.name, cif.name)
                self.assertEqual(private.read_bytes(), cif.read_bytes())
                self.assertNotEqual(model_directory, models)
                self.assertEqual((model_directory/'fake.pkl').read_bytes(), payload)
                private.write_bytes(b'backend change must remain isolated')
                (models/'fake.pkl').write_bytes(b'changed source after validation')
                self.assertEqual((model_directory/'fake.pkl').read_bytes(), payload)
                observed.append(private)
                return {'thermal stability': 0., 'water probability': 0.}

            with patch.dict(contract.ASSET_SHA256, {'fake.pkl':hashlib.sha256(payload).hexdigest()}, clear=True), \
                 patch.object(prediction, '_stability_prediction', side_effect=infer):
                result = prediction.stability(cif, model_directory=models)
            self.assertEqual(result['thermal stability'], 0.)
            self.assertEqual(cif.read_bytes(), b'unchanged CIF fixture')
            self.assertFalse(observed[0].exists())

    def test_source_directory_is_not_a_cif(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(ValueError, 'one CIF file'):
                prediction.stability(temporary)

    def test_model_destination_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root/'source'
            source.mkdir()
            (source/'fake.pkl').write_bytes(b'fixture')
            target = root/'destination'
            target.mkdir()
            with patch.dict(contract.ASSET_SHA256, {'fake.pkl':hashlib.sha256(b'fixture').hexdigest()}, clear=True):
                with self.assertRaises(FileExistsError):
                    contract.copy_verified_models(source, target)

    def test_nonfinite_and_nonnumeric_inputs_are_rejected(self):
        for value in (float('nan'), float('inf'), -float('inf'), True, False, '1', None):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'finite numeric'):
                contract.validate_values([value], 1, 'input')

    def test_wrong_feature_count_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'exactly 148'):
            contract.validate_values([0.] * 147, 148, 'thermal')

    def test_finite_values_and_zero_are_not_changed(self):
        values = [0., -2.25, 100.5]
        before = values.copy()
        contract.validate_values(values, 3, 'thermal')
        self.assertEqual(values, before)
        contract.validate_values([0., 1.], 2, 'probability', probability=True)

    def test_probabilities_are_not_clipped(self):
        for value in (-0.01, 1.01):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'outside'):
                contract.validate_values([value], 1, 'water', probability=True)


if __name__ == '__main__':
    unittest.main()
