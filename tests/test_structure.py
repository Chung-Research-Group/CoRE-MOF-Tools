import unittest

from CoREMOF.structure import information, read_aif


class StructureDataTests(unittest.TestCase):
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
