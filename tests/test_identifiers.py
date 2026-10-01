"""CoRE IDs identify records, never their source or their leakage groups."""
import csv
from dataclasses import FrozenInstanceError
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from urllib.parse import quote, unquote

from CoREMOF import CoREID, format_core_id, parse_core_id, resolve_cif_inputs
from CoREMOF.dataset import CoREMOFDataset, ReleaseValidationError
from CoREMOF.identifiers import validate_source_database
from tests.test_dataset_labels import _make_release, _write_csv


class IdentifierTests(unittest.TestCase):
    def test_paper_example_meanings_and_roundtrip(self):
        value = parse_core_id("2013[Cu][nan]3[ASR]5")
        self.assertEqual(value, CoREID(2013, "Cu", "nan", 3, "ASR", 5))
        self.assertEqual(str(value), "2013[Cu][nan]3[ASR]5")
        with self.assertRaises(FrozenInstanceError):
            value.serial = 7

    def test_metalloids_unknown_year_and_full_topology(self):
        for element in ("Si", "Ge", "Sb", "CuCo", "CoCu", "FeZn"):
            for dimension in range(4):
                value = format_core_id(0, element, "pts-x", dimension, "ION", 12345)
                parsed = parse_core_id(value)
                self.assertEqual(parsed.elements, element)
                self.assertEqual(parsed.topology, "pts-x")
                self.assertEqual(parsed.dimension, dimension)
                self.assertEqual(parsed.year, 0)

    def test_reject_unsafe_or_noncanonical_names(self):
        for value in (None, "", " 2013[Cu][pcu]3[ASR]1", "2013[Cu][pcu]3[ASR]1\n",
                      "2013[Cu][pcu]3[ASR]1.cif", "../2013[Cu][pcu]3[ASR]1",
                      "2013[Cu][pcu/abc]3[ASR]1", "2013[Cu][pcu]4[ASR]1",
                      "2013[Cu][pcu]3[ASR]0", "2013[Cu][pcu]3[ASR]01",
                      "2013[Xx][pcu]3[ASR]1", "2013[CuCu][pcu]3[ASR]1",
                      "2013[Cu][pcu]3[CR]1", "2013[Cu][pcu]3[Ion]1"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                parse_core_id(value)
        for year, dim, serial in ((True, 3, 1), (2026, True, 1), (2026, 3, True),
                                   (10000, 3, 1), (2026, 3, 0)):
            with self.assertRaises(ValueError):
                format_core_id(year, "Cu", "pcu", dim, "ASR", serial)

    def test_brackets_are_exact_files_and_url_escaped(self):
        name = "2013[Cu][nan]3[ASR]5.cif"
        self.assertEqual(unquote(quote(name, safe="")), name)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / name
            path.write_text("data_example\n")
            self.assertEqual(path.read_text(), "data_example\n")
            self.assertEqual(resolve_cif_inputs(path), (path,))
            self.assertEqual(resolve_cif_inputs(directory), (path,))

    def test_source_is_explicit_metadata_not_a_name_segment(self):
        self.assertEqual(validate_source_database("CSD"), "CSD")
        for value in (None, "", "csd", "CSD "):
            with self.assertRaises(ValueError):
                validate_source_database(value)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _make_release(root)
            path = root / "metadata/metadata.csv"
            with path.open() as stream:
                rows = list(csv.DictReader(stream))
            rows[0]["source_database"] = "SI"
            _write_csv(path, tuple(rows[0]), rows)
            dataset = CoREMOFDataset.from_release(root)
            self.assertEqual(dataset.records[0].metadata["source_database"], "SI")
            rows[0]["source_database"] = "UNKNOWN"
            _write_csv(path, tuple(rows[0]), rows)
            with self.assertRaisesRegex(ReleaseValidationError, "source_database"):
                CoREMOFDataset.from_release(root)

    def test_standard_library_only(self):
        result = subprocess.run(
            [sys.executable, "-B", "-S", "-c",
             "from CoREMOF import parse_core_id; "
             "assert parse_core_id('2013[Cu][nan]3[ASR]5').dimension == 3"],
            cwd=Path(__file__).resolve().parents[1], capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
