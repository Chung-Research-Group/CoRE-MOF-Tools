import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "examples" / "CoREMOF_dataset_splitting_quickstart.ipynb"


class DatasetSplittingNotebookTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))

    def test_notebook_is_clean_nbformat4(self):
        self.assertEqual(self.notebook["nbformat"], 4)
        self.assertGreaterEqual(self.notebook["nbformat_minor"], 5)
        cell_ids = [cell.get("id") for cell in self.notebook["cells"]]
        self.assertNotIn(None, cell_ids)
        self.assertEqual(len(cell_ids), len(set(cell_ids)))
        for cell_id in cell_ids:
            self.assertRegex(cell_id, r"^[A-Za-z0-9_-]{1,64}$")
        code_cells = [
            cell for cell in self.notebook["cells"] if cell["cell_type"] == "code"
        ]
        self.assertGreaterEqual(len(code_cells), 10)
        for cell in code_cells:
            self.assertIsNone(cell["execution_count"])
            self.assertEqual(cell["outputs"], [])

    def test_every_code_cell_compiles(self):
        for index, cell in enumerate(self.notebook["cells"]):
            if cell["cell_type"] != "code":
                continue
            source = "".join(cell["source"])
            compile(source, "{}#cell-{}".format(NOTEBOOK, index), "exec")

    def test_notebook_covers_safety_contract(self):
        text = "\n".join(
            "".join(cell["source"]) for cell in self.notebook["cells"]
        )
        for required in (
            "priority_main",
            "main_union",
            "verify_cif_files=RUN_STRICT_CIF_CHECK",
            "AMBIGUOUS",
            "UNCHECKED",
            "cross_split_block_count",
            "assignment_sha256",
            "provisional_input",
            "official_split",
            "licence-sanitization",
        ):
            self.assertIn(required, text)

    def test_package_documentation_links_notebook(self):
        relative = "examples/CoREMOF_dataset_splitting_quickstart.ipynb"
        public_paths = (
            ROOT / "README.md",
            ROOT / "README_DATASET_SPLITTING.md",
            ROOT / "examples" / "README.md",
            NOTEBOOK,
        )
        for path in public_paths[:2]:
            self.assertIn(relative, path.read_text(encoding="utf-8"))
        self.assertIn(
            "CoREMOF_dataset_splitting_quickstart.ipynb",
            (ROOT / "examples" / "README.md").read_text(encoding="utf-8"),
        )
        for path in public_paths:
            text = path.read_text(encoding="utf-8")
            self.assertNotIn("/home/yuc", text)
            self.assertNotIn("CoREMOF-COD", text)


if __name__ == "__main__":
    unittest.main()
