"""Science-free regressions for curation/checker boundary validation.

Load the actual definitions without importing licensed CCDC software or
initializing optional scientific runtimes. Chemistry and predictors are stubs;
the validation, output selection and aggregation code under test is unchanged.
"""

import ast
import collections
from collections.abc import Mapping
import copy
import csv
import itertools
import math
import os
from pathlib import Path
import re
import shutil
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch
import warnings


ROOT = Path(__file__).resolve().parents[1]


def _definitions(relative_path, names, namespace):
    tree = ast.parse((ROOT / relative_path).read_text(encoding="utf-8"))
    selected = [node for node in tree.body if getattr(node, "name", None) in names]
    exec(compile(ast.Module(body=selected, type_ignores=[]), relative_path, "exec"), namespace)
    return SimpleNamespace(**{name: namespace[name] for name in names})


CURATE_GLOBALS = dict(
    collections=collections, Mapping=Mapping, csv=csv, itertools=itertools,
    math=math, os=os, re=re, shutil=shutil, tempfile=tempfile,
    warnings=warnings, np=SimpleNamespace(bool_=bool), pd=SimpleNamespace(DataFrame=object),
)
curate = _definitions(
    "CoREMOF/curate.py", {"clean_pacman", "_validated_pacman_charges"},
    CURATE_GLOBALS,
)

class FakeColumn(list):
    def __init__(self, values, tags):
        super().__init__(values)
        self.tags = tags

    def get_loop(self):
        return SimpleNamespace(tags=self.tags)


class FakeBlock:
    def __init__(self, charges=None):
        self.columns = {
            "_atom_site_label": ["Zn1", "C1", "C2", "C3", "Cl1"],
            "_atom_site_type_symbol": ["Zn", "C", "C", "C", "Cl"],
            "_atom_site_fract_x": ["0", "0.1", "0.2", "0.3", "0.4"],
            "_atom_site_fract_y": ["0"] * 5,
            "_atom_site_fract_z": ["0"] * 5,
            "_atom_site_occupancy": ["1"] * 5,
            "_space_group_symop_operation_xyz": ["x,y,z"],
        }
        if charges is not None:
            self.columns["_atom_site_charge"] = charges
        self.scalars = {
            "_cell_length_a": "10", "_cell_length_b": "10", "_cell_length_c": "10",
            "_cell_angle_alpha": "90", "_cell_angle_beta": "90", "_cell_angle_gamma": "90",
        }

    def find_loop(self, tag):
        return FakeColumn(self.columns.get(tag, []), list(self.columns))

    def find_value(self, tag):
        return self.scalars.get(tag)

    def find_values(self, tag):
        return self.columns.get(tag, [self.scalars[tag]] if tag in self.scalars else [])


class FakeAtoms:
    symbols = ["Zn", "C", "C", "C", "Cl"]

    def __len__(self):
        return len(self.symbols)

    def __getitem__(self, index):
        return tuple(index) if isinstance(index, list) else SimpleNamespace(symbol=self.symbols[index])

    def get_chemical_symbols(self):
        return list(self.symbols)

    def get_scaled_positions(self, wrap=False):
        return [[i / 10, 0, 0] for i in range(len(self))]


class PacmanCurationSafetyTests(unittest.TestCase):
    def setUp(self):
        self.source = FakeBlock()
        self.charged = FakeBlock(["0"] * 5)
        self.atoms = FakeAtoms()
        self.writer = Mock()
        cif = SimpleNamespace(
            as_number=float,
            read_file=lambda path: SimpleNamespace(
                sole_block=lambda: self.charged if "pacman" in str(path) else self.source
            ),
        )
        self.environment = patch.dict(CURATE_GLOBALS, CIF=cif, read=lambda _: self.atoms, write=self.writer)
        self.environment.start()
        self.addCleanup(self.environment.stop)
        self.cleaner = curate.clean_pacman.__new__(curate.clean_pacman)
        self.cleaner.build_ASE_neighborlist = Mock(return_value=None)
        self.cleaner.CustomMatrix = Mock(return_value=None)
        self.cleaner.find_clusters = Mock(return_value=[[0, 1, 2, 3], [4]])

    def validate(self):
        return curate._validated_pacman_charges("original.cif", "original_pacman.cif", self.atoms)

    def test_empty_short_extra_missing_and_nonfinite_charges_never_write_a_cif(self):
        for values in ([], ["0"] * 4, ["0"] * 6, ["0"] * 4 + ["nan"],
                       ["0"] * 4 + ["inf"], ["0"] * 4 + ["-inf"], ["0"] * 4 + ["?"]):
            with self.subTest(charges=values):
                self.charged.columns["_atom_site_charge"] = values
                with self.assertRaises((TypeError, ValueError)):
                    self.cleaner.free_clean("original", "unused", 0.25)
                self.writer.assert_not_called()
                self.cleaner.build_ASE_neighborlist.assert_not_called()

    def test_complete_neutral_and_ionic_components_keep_native_selection(self):
        self.assertEqual(self.validate(), [0.0] * 5)
        result = self.cleaner.free_clean("original", "unused", 0.25)
        self.assertEqual(result, (0.25, ["Cl"], [], []))
        self.assertEqual(self.writer.call_args.args[1], (0, 1, 2, 3))
        self.charged.columns["_atom_site_charge"][-1] = "-1"
        result = self.cleaner.free_clean("original", "unused", 0.25)
        self.assertEqual(result, (0.25, [], ["Cl"], [-1.0]))
        self.assertEqual(self.writer.call_args.args[1], (0, 1, 2, 3, 4))
        self.assertTrue(self.writer.call_args.args[0].endswith("_FSR_ION.cif"))

    def test_atom_order_chemistry_coordinates_occupancy_cell_and_symmetry_must_match(self):
        original = copy.deepcopy(self.charged)
        mutations = (
            lambda b: b.columns["_atom_site_label"].reverse(),
            lambda b: b.columns["_atom_site_type_symbol"].__setitem__(0, "Cu"),
            lambda b: b.columns["_atom_site_fract_x"].__setitem__(0, "0.01"),
            lambda b: b.columns["_atom_site_occupancy"].__setitem__(0, "0.5"),
            lambda b: b.scalars.__setitem__("_cell_length_a", "11"),
            lambda b: b.columns["_space_group_symop_operation_xyz"].append("-x,-y,-z"),
        )
        for mutate in mutations:
            self.charged = copy.deepcopy(original)
            mutate(self.charged)
            with self.subTest(mutation=mutate), self.assertRaises(ValueError):
                self.validate()
        self.writer.assert_not_called()

    def test_symmetry_expanded_or_reordered_ase_atoms_are_not_indexed_by_raw_rows(self):
        self.atoms.get_scaled_positions = lambda **kwargs: list(reversed(FakeAtoms().get_scaled_positions()))
        with self.assertRaisesRegex(ValueError, "unsupported ASE atom mapping"):
            self.validate()
        self.atoms.symbols = FakeAtoms.symbols + ["Zn"]
        with self.assertRaisesRegex(ValueError, "6 values"):
            self.validate()

    def test_charge_column_in_a_separate_loop_has_no_proven_atom_index(self):
        original_find = self.charged.find_loop
        self.charged.find_loop = lambda tag: (
            FakeColumn(["0"] * 5, ["_atom_site_charge"])
            if tag == "_atom_site_charge" else original_find(tag)
        )
        with self.assertRaisesRegex(ValueError, "share one atom-site loop"):
            self.validate()
        self.writer.assert_not_called()

    def test_predictor_receives_an_isolated_copy_and_validation_precedes_public_output(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            source = base / "original.cif"
            source.write_text("original bytes", encoding="utf-8")
            output = base / "outputs"
            output.mkdir()
            self.cleaner.structure = str(source)
            self.cleaner.output = str(output)
            seen = []

            def predict(**kwargs):
                isolated = Path(kwargs["cif_file"])
                seen.append(isolated)
                self.assertNotEqual(isolated, source)
                self.assertEqual(isolated.read_text(), "original bytes")
                isolated.write_text("predictor modifies its input", encoding="utf-8")
                isolated.with_name("input_pacman.cif").write_text("charged bytes", encoding="utf-8")

            with patch.dict(CURATE_GLOBALS, pmcharge=SimpleNamespace(predict=predict)):
                with patch.dict(CURATE_GLOBALS, _validated_pacman_charges=Mock(side_effect=ValueError("bad mapping"))):
                    with self.assertRaisesRegex(ValueError, "bad mapping"):
                        self.cleaner.run_pacman()
                self.assertEqual(list(output.iterdir()), [])
                with patch.dict(CURATE_GLOBALS, _validated_pacman_charges=Mock(return_value=[0.0] * 5)):
                    self.cleaner.run_pacman()
            self.assertEqual(source.read_text(), "original bytes")
            self.assertEqual((output / "original_pacman.cif").read_text(), "charged bytes")
            self.assertTrue(all(not isolated.parent.exists() for isolated in seen))

    def test_asr_is_explicitly_unavailable_and_never_emits_an_asr_cif(self):
        self.cleaner.structure = "original.cif"
        self.cleaner.saveto = False
        self.cleaner.run_fsr = Mock(return_value=(0.25, [], [], []))
        with self.assertWarnsRegex(RuntimeWarning, "ASR_UNSUPPORTED"):
            self.cleaner.process()
        self.assertEqual(self.cleaner.asr_status, "NOT_AVAILABLE")
        with self.assertRaises(NotImplementedError):
            self.cleaner.run_asr()
        with self.assertRaises(NotImplementedError):
            self.cleaner.all_clean("original", "unused", 0.25)
        self.writer.assert_not_called()

    def test_process_propagates_invalid_fsr_instead_of_reporting_completion(self):
        self.cleaner.run_fsr = Mock(side_effect=ValueError("invalid charges"))
        with self.assertRaisesRegex(ValueError, "invalid charges"):
            self.cleaner.process()

    def test_native_p1_cif_charge_rows_match_ase_atoms_without_changing_bytes(self):
        try:
            from ase.io import read as read_atoms
            from gemmi import cif
        except ImportError:
            self.skipTest("optional ASE/gemmi parser integration")
        header = (
            "data_test\n_cell_length_a 10\n_cell_length_b 10\n_cell_length_c 10\n"
            "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
            "_symmetry_space_group_name_H-M 'P 1'\n"
            "loop_\n_space_group_symop_operation_xyz\nx,y,z\n"
            "loop_\n_atom_site_label\n_atom_site_type_symbol\n"
            "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n"
        )
        rows = ("Zn1 Zn 0 0 0 1", "C1 C 0.1 0 0 1", "C2 C 0.2 0 0 1",
                "C3 C 0.3 0 0 1", "Cl1 Cl 0.4 0 0 1")
        source_text = header + "\n".join(rows) + "\n"
        charged_text = header + "_atom_site_charge\n" + "\n".join(
            row + (" -1" if index == 4 else " 0") for index, row in enumerate(rows)
        ) + "\n"
        with tempfile.TemporaryDirectory() as directory:
            source, charged = Path(directory) / "source.cif", Path(directory) / "charged.cif"
            source.write_text(source_text, encoding="utf-8")
            charged.write_text(charged_text, encoding="utf-8")
            original_bytes, charged_bytes = source.read_bytes(), charged.read_bytes()
            atoms = read_atoms(source)
            with patch.dict(CURATE_GLOBALS, CIF=cif):
                result = curate._validated_pacman_charges(source, charged, atoms)
            self.assertEqual(result, [0.0, 0.0, 0.0, 0.0, -1.0])
            self.assertEqual(atoms.get_chemical_symbols(), ["Zn", "C", "C", "C", "Cl"])
            self.assertEqual(source.read_bytes(), original_bytes)
            self.assertEqual(charged.read_bytes(), charged_bytes)






if __name__ == "__main__":
    unittest.main()
