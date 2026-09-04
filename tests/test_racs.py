import hashlib
import json
from pathlib import Path
import sys
from types import ModuleType
import unittest
from unittest.mock import patch


def _scientific_import_stubs():
    modules = {}
    for name in (
        "ase",
        "ase.io",
        "pymatgen",
        "pymatgen.core",
        "pymatgen.core.structure",
        "pymatgen.io",
        "pymatgen.io.ase",
        "pymatgen.symmetry",
        "pymatgen.symmetry.analyzer",
    ):
        module = ModuleType(name)
        module.__path__ = []
        modules[name] = module
    modules["ase.io"].read = object()
    modules["pymatgen.core.structure"].Structure = object
    modules["pymatgen.io.ase"].AseAtomsAdaptor = object
    modules["pymatgen.symmetry.analyzer"].SpacegroupAnalyzer = object
    return modules


with patch.dict(sys.modules, _scientific_import_stubs()):
    from CoREMOF.calculation import mof_features


def _flatten_names(depth):
    groups = mof_features._rac_names_by_group(depth)
    return [name for group_names in groups.values() for name in group_names]


def _fake_molsimplify(get_descriptors):
    modules = {}
    for name in (
        "molSimplify",
        "molSimplify.Informatics",
        "molSimplify.Informatics.MOF",
        "molSimplify.Informatics.MOF.MOF_descriptors",
    ):
        module = ModuleType(name)
        module.__path__ = []
        modules[name] = module
    modules[
        "molSimplify.Informatics.MOF.MOF_descriptors"
    ].get_MOF_descriptors = get_descriptors
    return patch.dict(sys.modules, modules)


class RacDepthTests(unittest.TestCase):
    structure = Path("structure with spaces.cif")

    def _run(self, depth=None, result_transform=None):
        calls = []

        def get_descriptors(structure, requested_depth, **kwargs):
            recorded_kwargs = dict(kwargs)
            recorded_kwargs["_path_was_directory"] = Path(kwargs["path"]).is_dir()
            calls.append((structure, requested_depth, recorded_kwargs))
            names = _flatten_names(requested_depth)
            values = [index + 0.123456 for index in range(len(names))]
            if result_transform is not None:
                names, values = result_transform(names, values)
            return names, values

        with _fake_molsimplify(get_descriptors):
            if depth is None:
                result = mof_features.RACs(self.structure)
            else:
                result = mof_features.RACs(self.structure, depth=depth)
        return result, calls

    def test_default_depth_three_matches_explicit_historical_shape_and_values(self):
        implicit, implicit_calls = self._run()
        explicit, explicit_calls = self._run(depth=3)

        self.assertEqual(implicit, explicit)
        self.assertEqual(implicit_calls[0][1], 3)
        self.assertEqual(explicit_calls[0][1], 3)
        self.assertEqual(list(implicit), ["Metal", "Linker", "Function-group"])
        self.assertEqual([len(implicit[group]) for group in implicit], [60, 68, 48])
        self.assertEqual(next(iter(implicit["Metal"])), "D_mc-I-0-all")
        self.assertEqual(next(reversed(implicit["Metal"])), "mc-chi-3-all")
        self.assertEqual(next(iter(implicit["Linker"])), "D_lc-I-0-all")
        self.assertEqual(next(reversed(implicit["Linker"])), "f-lig-chi-3")
        self.assertEqual(next(iter(implicit["Function-group"])), "D_func-I-0-all")
        self.assertEqual(next(reversed(implicit["Function-group"])), "func-chi-3-all")
        self.assertEqual(implicit["Metal"]["D_mc-I-0-all"], 0.1235)
        serialized = json.dumps(
            implicit, ensure_ascii=False, separators=(",", ":")
        ).encode("utf-8")
        self.assertEqual(len(serialized), 4123)
        self.assertEqual(
            hashlib.sha256(serialized).hexdigest(),
            "09543d1b57cadc3fdfa86f86da9fa2de5530d24a25eeee8080e32015df2b6d23",
        )

    def test_requested_depth_is_forwarded_and_schema_is_derived(self):
        result, calls = self._run(depth=5)

        self.assertEqual(calls[0][1], 5)
        self.assertEqual(sum(len(group) for group in result.values()), 264)
        self.assertIn("D_mc-I-5-all", result["Metal"])
        self.assertIn("f-lig-chi-5", result["Linker"])
        self.assertIn("func-alpha-5-all", result["Function-group"])
        self.assertTrue(calls[0][2]["_path_was_directory"])
        self.assertEqual(calls[0][2]["max_num_atoms"], 6000)

    def test_non_bool_integer_and_nonnegative_depth_are_required(self):
        for invalid in (True, False, 3.0, "3", None):
            with self.subTest(depth=invalid):
                with self.assertRaises(TypeError):
                    mof_features.RACs(self.structure, depth=invalid)
        with self.assertRaises(ValueError):
            mof_features.RACs(self.structure, depth=-1)

    def test_descriptor_order_from_molsimplify_does_not_change_public_order(self):
        result, _calls = self._run(
            depth=1,
            result_transform=lambda names, values: (
                list(reversed(names)),
                list(reversed(values)),
            ),
        )

        self.assertEqual(next(iter(result["Metal"])), "D_mc-I-0-all")
        self.assertEqual(next(reversed(result["Function-group"])), "func-chi-1-all")
        self.assertEqual(sum(len(group) for group in result.values()), 88)

    def test_malformed_or_nonfinite_descriptor_output_fails_closed(self):
        def duplicate_name(names, values):
            names[-1] = names[0]
            return names, values

        with self.assertRaisesRegex(ValueError, "duplicate RAC names"):
            self._run(depth=0, result_transform=duplicate_name)

        def missing_value(names, values):
            return names, values[:-1]

        with self.assertRaisesRegex(ValueError, "different RAC name/value counts"):
            self._run(depth=0, result_transform=missing_value)

        def nonfinite_value(names, values):
            values[0] = float("nan")
            return names, values

        with self.assertRaisesRegex(ValueError, "not finite"):
            self._run(depth=0, result_transform=nonfinite_value)


if __name__ == "__main__":
    unittest.main()
