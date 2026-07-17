from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

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
