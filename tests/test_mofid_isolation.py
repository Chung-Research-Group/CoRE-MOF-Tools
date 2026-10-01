"""Wrapper safety without invoking MOFid, Open Babel, Java or a matcher.

The function bodies are compiled from the maintained module. Scientific
dependencies are explicit stubs, so these are not scientific-equivalence tests.
"""
import ast
from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
import tempfile
import threading
from types import SimpleNamespace
import unittest
from unittest.mock import Mock


def wrappers():
    source = Path(__file__).resolve().parents[1] / "CoREMOF/get_mofid.py"
    tree = ast.parse(source.read_text(encoding="utf-8"))
    functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)
                 and node.name in ("run_v1", "run_v2")]
    namespace = {"Path": Path, "tempfile": tempfile}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(source), "exec"), namespace)
    return namespace


class MofidIsolationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="mofid wrapper test ")
        self.root = Path(self.temporary.name)
        self.original_cwd = Path.cwd()
        os.chdir(self.root)
        self.addCleanup(self.temporary.cleanup)
        self.addCleanup(os.chdir, self.original_cwd)
        self.structure = self.root / "MOF with spaces.cif"
        self.structure.write_text("data_untouched\n")
        self.sentinel = self.root / "Output/sentinel"
        self.sentinel.parent.mkdir()
        self.sentinel.write_text("unrelated output\n")
        self.library = self.root / "node library"
        self.library.mkdir()
        (self.library / "Zn1_Type-1.xyz").write_text("reference node\n")
        self.calls = []
        self.ns = wrappers()

        def cif2mofid(structure, output_path):
            self.calls.append((structure, Path(output_path)))
            output = Path(output_path)
            (output / "AllNode").mkdir()
            (output / "AllNode/nodes.cif").write_text("mock extracted nodes\n")
            return {"smiles_linkers": ["C"], "topology": "pcu", "cat": "0"}

        def split_nodes(structure, prefix):
            self.assertTrue(Path(structure).is_file())
            (Path(prefix) / "node0.xyz").write_text("mock node\n")
            return 0

        self.matcher = Mock()
        self.matcher.fit.return_value = True
        self.ns.update(
            cif2mofid=cif2mofid, split_nodes_from_cif=split_nodes,
            xyz2fomula=lambda path: "Zn1", ase_read=lambda path: str(path),
            remove_pbc_cuts=lambda atoms: atoms, convert_ase_pymat=lambda atoms: atoms,
            StructureMatcher=Mock(return_value=self.matcher), ElementComparator=Mock(),
            sf=SimpleNamespace(encoder=lambda smiles: smiles),
        )

    def library_contents(self):
        return {p.name: p.read_bytes() for p in self.library.iterdir()}

    def test_v1_private_output_is_cleaned_and_caller_output_is_preserved(self):
        result = self.ns["run_v1"](self.structure)
        self.assertEqual(result["topology"], "pcu")
        self.assertFalse(self.calls[0][1].exists())
        self.assertEqual(self.sentinel.read_text(), "unrelated output\n")
        self.assertEqual(self.structure.read_text(), "data_untouched\n")

    def test_v1_explicit_outputs_are_retained_and_never_overwrite(self):
        destination = self.root / "saved fragments"
        self.ns["run_v1"](self.structure, output_path=destination)
        self.assertTrue((destination / "AllNode/nodes.cif").is_file())
        with self.assertRaises(FileExistsError):
            self.ns["run_v1"](self.structure, output_path=destination)
        with self.assertRaises(FileExistsError):
            self.ns["run_v1"](self.structure, output_path=self.sentinel.parent)
        self.assertEqual(self.sentinel.read_text(), "unrelated output\n")

    def test_v2_valid_identifier_does_not_mutate_input_or_library(self):
        before = self.library_contents()
        result = self.ns["run_v2"](self.structure, self.library, "example")
        self.assertEqual(result, "[Zn1_Type-1].C MOFid-v2.pcu.cat0;example")
        self.assertEqual(before, self.library_contents())
        self.assertFalse(self.calls[0][1].parent.exists())
        self.assertEqual(self.sentinel.read_text(), "unrelated output\n")
        settings = self.ns["StructureMatcher"].call_args.kwargs
        self.assertEqual((settings["ltol"], settings["stol"]), (0.3, 2))

    def test_v2_unmatched_never_invents_or_moves_a_node(self):
        self.matcher.fit.return_value = False
        before = self.library_contents()
        with self.assertRaisesRegex(ValueError, "NOT_AVAILABLE_UNMATCHED_NODE"):
            self.ns["run_v2"](self.structure, self.library, "example")
        self.assertEqual(before, self.library_contents())
        self.assertFalse(self.calls[0][1].parent.exists())
        self.assertTrue(self.sentinel.is_file())

    def test_v2_unseen_formula_never_adds_type_one(self):
        self.ns["xyz2fomula"] = lambda path: "Cu2"
        before = self.library_contents()
        with self.assertRaisesRegex(ValueError, "NOT_AVAILABLE_UNMATCHED_NODE"):
            self.ns["run_v2"](self.structure, self.library, "example")
        self.assertEqual(before, self.library_contents())
        self.matcher.fit.assert_not_called()

    def test_v2_ambiguous_checks_all_candidates_and_returns_no_identifier(self):
        (self.library / "Zn1_Type-2.xyz").write_text("second node\n")
        before = self.library_contents()
        with self.assertRaisesRegex(ValueError, "NOT_AVAILABLE_AMBIGUOUS_NODE"):
            self.ns["run_v2"](self.structure, self.library, "example")
        self.assertEqual(self.matcher.fit.call_count, 2)
        self.assertEqual(before, self.library_contents())

    def test_extraction_failure_is_explicit_and_cleans_only_owned_output(self):
        self.ns["split_nodes_from_cif"] = lambda *args: 1
        with self.assertRaisesRegex(RuntimeError, "extraction is unavailable"):
            self.ns["run_v2"](self.structure, self.library, "example")
        self.assertFalse(self.calls[0][1].parent.exists())
        self.assertTrue(self.sentinel.is_file())

    def test_upstream_exception_cleans_only_owned_output(self):
        def fail(structure, output_path):
            self.calls.append((structure, Path(output_path)))
            raise RuntimeError("upstream failed")
        self.ns["cif2mofid"] = fail
        with self.assertRaisesRegex(RuntimeError, "upstream failed"):
            self.ns["run_v2"](self.structure, self.library, "example")
        self.assertFalse(self.calls[0][1].parent.exists())
        self.assertTrue(self.sentinel.is_file())

    def test_concurrent_calls_have_independent_outputs_without_chdir(self):
        barrier = threading.Barrier(2)
        original = self.ns["cif2mofid"]

        def parallel(structure, output_path):
            result = original(structure, output_path)
            barrier.wait(timeout=5)
            return result

        self.ns["cif2mofid"] = parallel
        with ThreadPoolExecutor(max_workers=2) as executor:
            futures = [executor.submit(self.ns["run_v2"], self.structure, self.library, f"case{i}")
                       for i in range(2)]
            results = [future.result(timeout=10) for future in futures]
        self.assertTrue(results[0].endswith(";case0"))
        self.assertTrue(results[1].endswith(";case1"))
        self.assertNotEqual(self.calls[0][1], self.calls[1][1])
        self.assertTrue(all(not path.parent.exists() for _, path in self.calls))
        self.assertEqual(Path.cwd(), self.root)
        self.assertTrue(self.sentinel.is_file())


if __name__ == "__main__":
    unittest.main()
