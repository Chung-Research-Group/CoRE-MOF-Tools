import os
from pathlib import Path
import tempfile
import unittest

from CoREMOF import collect_cifs, resolve_cif_inputs


class CifInputResolverTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)

    def tearDown(self):
        self.temporary_directory.cleanup()

    def test_single_file_has_stable_collection_return_type(self):
        structure = self.root / "structure with spaces.CIF"
        structure.write_text("data_test\n", encoding="utf-8")

        self.assertEqual(resolve_cif_inputs(structure), (structure,))
        self.assertEqual(collect_cifs(structure), (structure,))

    def test_directory_is_nonrecursive_and_lexically_sorted(self):
        expected = [self.root / name for name in ("a.cif", "b.cif", "c.CIF")]
        for structure in reversed(expected):
            structure.write_text("data_test\n", encoding="utf-8")
        (self.root / "notes.txt").write_text("ignored\n", encoding="utf-8")
        nested = self.root / "nested"
        nested.mkdir()
        (nested / "not_collected.cif").write_text("data_nested\n", encoding="utf-8")

        self.assertEqual(resolve_cif_inputs(self.root), tuple(expected))

    def test_missing_unsupported_and_empty_inputs_are_rejected(self):
        with self.assertRaises(FileNotFoundError):
            resolve_cif_inputs(self.root / "missing.cif")

        unsupported = self.root / "structure.txt"
        unsupported.write_text("data_test\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, r"\.cif suffix"):
            resolve_cif_inputs(unsupported)

        empty = self.root / "empty"
        empty.mkdir()
        with self.assertRaisesRegex(ValueError, "no direct CIF files"):
            resolve_cif_inputs(empty)

    def test_symlinks_and_special_cif_paths_are_rejected(self):
        structure = self.root / "structure.cif"
        structure.write_text("data_test\n", encoding="utf-8")
        linked_input = self.root / "linked-input.cif"
        linked_input.symlink_to(structure)
        with self.assertRaisesRegex(ValueError, "must not be a symlink"):
            resolve_cif_inputs(linked_input)

        symlink_directory = self.root / "symlink-directory"
        symlink_directory.mkdir()
        (symlink_directory / "linked-child.cif").symlink_to(structure)
        with self.assertRaisesRegex(ValueError, "contains a symlink"):
            resolve_cif_inputs(symlink_directory)

        special_directory = self.root / "special-directory"
        special_directory.mkdir()
        os.mkfifo(special_directory / "pipe.cif")
        with self.assertRaisesRegex(ValueError, "special CIF path"):
            resolve_cif_inputs(special_directory)

    def test_ambiguous_basenames_are_rejected(self):
        (self.root / "sample.cif").write_text("data_one\n", encoding="utf-8")
        (self.root / "SAMPLE.CIF").write_text("data_two\n", encoding="utf-8")

        with self.assertRaisesRegex(ValueError, "basenames collide"):
            resolve_cif_inputs(self.root)


if __name__ == "__main__":
    unittest.main()
