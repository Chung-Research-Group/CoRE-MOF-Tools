from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from CoREMOF import prediction

from CoREMOF.models.cp_app.predictions import (
    predict_Cv_ensemble_structure,
    predict_Cv_ensemble_structure_multitemperatures,
)


class ConstantModel:
    def __init__(self, value):
        self.value = value

    def predict(self, frame):
        return np.full(len(frame), self.value)


class PredictionTests(unittest.TestCase):
    def test_public_cp_uses_private_input_and_normalizes_temperature(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "example.CIF"
            source.write_bytes(b"unchanged synthetic CIF bytes")
            backend = types.ModuleType("CoREMOF.models.cp_app.featurizer")

            def featurize(path, *, verbose, saveto):
                self.assertNotEqual(Path(path), source)
                self.assertEqual(Path(path).name, source.name)
                self.assertEqual(Path(path).read_bytes(), source.read_bytes())
                Path(saveto).write_text("fixture\n")

            def predict(**kwargs):
                self.assertEqual(kwargs["structure_name"], "example.CIF")
                self.assertEqual(kwargs["temperatures"], [300])
                pd.DataFrame({"Cv_gravimetric_300_mean": [0.], "Cv_molar_300_mean": [0.],
                              "Cv_gravimetric_300_std": [0.], "Cv_molar_300_std": [0.]}).to_csv(
                                  kwargs["save_to"], index=False)

            backend.featurize_structure = featurize
            with patch.dict(sys.modules, {backend.__name__: backend}), \
                 patch("CoREMOF._heat_capacity.validate_ensemble") as validator, \
                 patch("CoREMOF.models.cp_app.predictions.predict_Cv_ensemble_structure_multitemperatures",
                       side_effect=predict):
                result = prediction.cp(source, T=[300.0], model_directory=Path(directory)/"models")
            validator.assert_called_once_with(Path(directory)/"models", [300])
            self.assertEqual(result["300_mean"], [0., 0.])
            self.assertEqual(source.read_bytes(), b"unchanged synthetic CIF bytes")

    def frame(self):
        return pd.DataFrame({"structure_name": ["a.cif", "a.cif"],
                             "x": [1., 2.], "site AtomicWeight": [10., 20.]})

    def test_no_atomic_rows_can_be_silently_dropped(self):
        for column in ("structure_name", "x", "site AtomicWeight"):
            with self.subTest(column=column):
                frame = self.frame()
                frame.loc[1, column] = np.nan
                with self.assertRaises(ValueError):
                    predict_Cv_ensemble_structure([ConstantModel(1.)], ["x"], frame, 300)

    def test_invalid_features_or_weights_are_rejected(self):
        for column, values in (("x", [np.inf, -np.inf, "invalid"]),
                               ("site AtomicWeight", [np.inf, -np.inf, 0., -1., "invalid"])):
            for value in values:
                with self.subTest(column=column, value=value):
                    frame = self.frame().astype({column: object})
                    frame.loc[1, column] = value
                    with self.assertRaises(ValueError):
                        predict_Cv_ensemble_structure([ConstantModel(1.)], ["x"], frame, 300)

    def test_model_predictions_must_be_complete_and_finite(self):
        for output in ([1., np.nan], [np.inf, 1.], [1.], [[1.], [2.]], 1.):
            with self.subTest(output=output):
                model = ConstantModel(1.)
                with patch.object(model, "predict", return_value=output):
                    with self.assertRaisesRegex(ValueError, "one finite value per atom"):
                        predict_Cv_ensemble_structure([model], ["x"], self.frame(), 300)

    def test_zero_values_and_valid_legacy_arithmetic_are_preserved(self):
        frame = self.frame()
        original = frame.copy(deep=True)
        result = predict_Cv_ensemble_structure([ConstantModel(0.), ConstantModel(6.)],
                                               ["x"], frame, 300)[0]
        self.assertEqual(result["Cv_molar_300_mean"], 3.)
        self.assertEqual(result["Cv_molar_300_std"], 3.)
        self.assertEqual(result["Cv_gravimetric_300_mean"], 0.2)
        self.assertEqual(result["Cv_gravimetric_300_std"], 0.2)
        pd.testing.assert_frame_equal(frame, original)

    def test_float32_predictions_keep_original_reduction_precision(self):
        frame = self.frame()
        values = np.asarray([10.123456, 5.234567], dtype=np.float32)
        model = ConstantModel(1.)
        with patch.object(model, "predict", return_value=values):
            result = predict_Cv_ensemble_structure([model], ["x"], frame, 300)[0]
        legacy_sum = np.sum(pd.Series(values))
        self.assertEqual(result["Cv_molar_300_mean"], legacy_sum / 2)
        self.assertEqual(result["Cv_gravimetric_300_mean"], legacy_sum / 30.)

    def test_nonfinite_ensemble_dispersion_is_not_returned(self):
        with np.errstate(over="ignore", invalid="ignore"), self.assertRaisesRegex(ValueError, "aggregation"):
            predict_Cv_ensemble_structure([ConstantModel(1e307), ConstantModel(-1e307)],
                                          ["x"], self.frame(), 300)

    def test_duplicate_columns_or_selected_features_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "distinct"):
            predict_Cv_ensemble_structure([ConstantModel(1.)], ["x", "x"], self.frame(), 300)
        frame = pd.concat([self.frame(), self.frame()[["x"]]], axis=1)
        with self.assertRaisesRegex(ValueError, "duplicate column"):
            predict_Cv_ensemble_structure([ConstantModel(1.)], ["x"], frame, 300)

    def test_single_structure_rejects_mixed_input(self):
        frame = pd.DataFrame(
            {
                "structure_name": ["a.cif", "b.cif"],
                "x": [1.0, 2.0],
                "site AtomicWeight": [10.0, 20.0],
            }
        )
        with self.assertRaisesRegex(ValueError, "exactly one structure"):
            predict_Cv_ensemble_structure([ConstantModel(1.0)], ["x"], frame, 300)

    def test_multitemperature_prediction_and_csv_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model_dir = root / "models" / "300"
            model_dir.mkdir(parents=True)
            (model_dir / "model_0").touch()
            (model_dir / "model_1").touch()
            features = root / "features.csv"
            output = root / "prediction.csv"
            pd.DataFrame(
                {
                    "structure_name": ["test.cif", "test.cif"],
                    "x": [1.0, 2.0],
                    "site AtomicWeight": [10.0, 20.0],
                }
            ).to_csv(features, index=False)

            with patch(
                "CoREMOF.models.cp_app.predictions.joblib.load",
                side_effect=[ConstantModel(3.0), ConstantModel(5.0)],
            ):
                result = predict_Cv_ensemble_structure_multitemperatures(
                    str(root / "models"),
                    "test.cif",
                    features_file=str(features),
                    FEATURES=["x"],
                    temperatures=[300],
                    save_to=str(output),
                )

            self.assertAlmostEqual(result.loc[0, "Cv_molar_300_mean"], 4.0)
            self.assertTrue(output.is_file())
            self.assertNotIn("Unnamed: 0", pd.read_csv(output).columns)

    def test_missing_structure_has_actionable_error(self):
        with tempfile.TemporaryDirectory() as directory:
            features = Path(directory, "features.csv")
            pd.DataFrame(
                {
                    "structure_name": ["present.cif"],
                    "x": [1.0],
                    "site AtomicWeight": [10.0],
                }
            ).to_csv(features, index=False)
            with self.assertRaisesRegex(ValueError, "Available structures: present.cif"):
                predict_Cv_ensemble_structure_multitemperatures(
                    directory,
                    "missing.cif",
                    features_file=str(features),
                    FEATURES=["x"],
                    temperatures=[300],
                    save_to=None,
                )


if __name__ == "__main__":
    unittest.main()
