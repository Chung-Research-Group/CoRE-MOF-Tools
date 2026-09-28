import unittest
from unittest.mock import patch

from CoREMOF.structure import information, read_aif


class StructureDataTests(unittest.TestCase):
    def setUp(self):
        # Exercise the legacy lookup and AIF parser without distributing or
        # downloading a real structure-resolved database in code-only CI.
        self.fixture = {
            "ASR": {"2020[Cu][sql]2[ASR]1": {"GEMC": (
                "data_synthetic_test\n_units_loading 'Molecules/Supercell'\n"
                "loop_\n_adsorp_pressure\n_adsorp_amount\n1 0\n2 3\n"
            )}},
            "FSR": {}, "Ion": {}, "unit": {"fixture": "synthetic"},
        }
        data = patch("CoREMOF.structure._load_json_data", return_value=self.fixture)
        network = patch("CoREMOF.structure.requests.get",
                        side_effect=AssertionError("Unit tests must not download data"))
        data.start()
        network.start()
        self.addCleanup(data.stop)
        self.addCleanup(network.stop)

    def test_known_database_record_can_be_loaded(self):
        record = information("CR-ASR", "2020[Cu][sql]2[ASR]1")
        self.assertIsInstance(record, dict)

    def test_embedded_adsorption_data_can_be_parsed(self):
        record = information("CR-ASR", "2020[Cu][sql]2[ASR]1")
        adsorption = read_aif(record["GEMC"])
        self.assertEqual(len(adsorption["pressure"]), len(adsorption["uptake"]))
        self.assertGreater(len(adsorption["pressure"]), 0)

    def test_invalid_dataset_lists_valid_choices(self):
        with self.assertRaisesRegex(ValueError, "CR-ASR"):
            information("invalid", "anything")

    def test_missing_entry_reports_close_match(self):
        with self.assertRaisesRegex(KeyError, "Close matches"):
            information("CR-ASR", "2020[Cu][sql]2[ASR]")


if __name__ == "__main__":
    unittest.main()
