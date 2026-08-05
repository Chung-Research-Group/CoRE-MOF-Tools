"""Static regression checks for the dataset-splitting README handbook."""

import re
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
HANDBOOK = ROOT / "README_DATASET_SPLITTING.md"


def _heading_slug(value):
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9 _-]", "", value)
    return re.sub(r"[ _]+", "-", value)


class HandbookTests(unittest.TestCase):
    def test_main_readme_links_the_handbook(self):
        self.assertTrue(HANDBOOK.is_file())
        main_readme = (ROOT / "README.md").read_text(encoding="utf-8")
        self.assertIn(
            "https://github.com/Chung-Research-Group/CoRE-MOF-Tools/"
            "blob/main/README_DATASET_SPLITTING.md",
            main_readme,
        )

    def test_python_examples_are_syntactically_valid(self):
        text = HANDBOOK.read_text(encoding="utf-8")
        self.assertEqual(text.count("```") % 2, 0)
        blocks = re.findall(r"```python\n(.*?)\n```", text, flags=re.DOTALL)
        self.assertGreaterEqual(len(blocks), 40)
        for index, block in enumerate(blocks, start=1):
            with self.subTest(block=index):
                compile(block, "handbook-python-block-{}".format(index), "exec")

    def test_table_of_contents_anchors_resolve(self):
        text = HANDBOOK.read_text(encoding="utf-8")
        headings = {
            _heading_slug(match.group(1))
            for match in re.finditer(r"^#{1,6}\s+(.+?)\s*$", text, re.MULTILINE)
        }
        anchors = re.findall(r"\]\(#([^)]+)\)", text)
        self.assertGreaterEqual(len(anchors), 20)
        self.assertEqual(sorted(set(anchors).difference(headings)), [])

    def test_project_defined_split_identifiers_are_defined_before_examples(self):
        documents = {
            ROOT / "README.md": "Four project-defined API identifiers",
            HANDBOOK: "Project-defined identifiers used in this handbook",
            ROOT / "docs" / "source" / "splitting.rst": (
                "Project-defined split identifiers"
            ),
            ROOT / "examples" / "README.md": (
                "project-defined CoREMOF-tools identifiers"
            ),
        }
        for path, marker in documents.items():
            with self.subTest(path=path.name):
                text = path.read_text(encoding="utf-8")
                definition_index = text.index(marker)
                first_example = min(
                    index
                    for token in (
                        'parent_method="priority_main"',
                        "--parent-method priority_main",
                    )
                    if (index := text.find(token)) >= 0
                )
                self.assertLess(definition_index, first_example)
                normalized = " ".join(
                    text.replace("source-sibling", "source sibling").split()
                )
                for required in (
                    "PARENT_METHOD_CONFLICT",
                    "rac_status/rac_group/rac_size",
                    "mofid2_status/mofid2_group/mofid2_size",
                    "mofid1_status/mofid1_group/mofid1_size",
                    "RAC5",
                    "MOFid-v2",
                    "MOFid-v1",
                    "full CIF SHA-256",
                    "source sibling",
                    "parent_only",
                    "StructureMatcher",
                ):
                    self.assertIn(required, normalized)

    def test_identifier_and_optional_group_terms_are_defined_at_first_use(self):
        documents = (
            HANDBOOK,
            ROOT / "docs" / "source" / "splitting.rst",
        )
        for path in documents:
            with self.subTest(path=path.name):
                text = path.read_text(encoding="utf-8")

                canonical_start = text.index("canonicalized identifier text")
                canonical_context = text[canonical_start : canonical_start + 4200]
                for required in (
                    "whitespace",
                    "Unicode NFKC",
                    "case-fold",
                    "source_database",
                    "source_id",
                    "complete",
                    "fuzzy",
                    ".cif.gz",
                    "_ASR_pacman",
                    "MOFidv2.",
                    ".no_ref",
                    "coordinate",
                ):
                    self.assertIn(required, canonical_context)

                identity_start = text.index("identity_union")
                identity_context = " ".join(
                    text[identity_start : identity_start + 2400].split()
                )
                for required in (
                    "provisional source-ID/MOFid transitive groups",
                    "identity_size",
                    "transitive connected component",
                    "not a count of edges or identifiers",
                    "v26.0.1",
                    "MOFid-v2",
                    "MOFid-v1",
                    "no priority",
                    "missing",
                    "RAC5",
                    "Zeo++",
                    "StructureMatcher",
                    "not proof",
                    "excluded from both",
                ):
                    self.assertIn(required, identity_context)

                for prefix, meaning in (
                    ("RT-", "RAC5"),
                    ("M2T-", "MOFid-v2"),
                    ("SM-", "connected"),
                ):
                    start = text.index(prefix)
                    context = text[start : start + 1700]
                    self.assertIn(meaning, context)
                    self.assertIn("group", context)

                for required in (
                    "264 finite",
                    "no scaling",
                    "density_g_cm3",
                    "current CrystalNets scientific fingerprint",
                    "network dimension",
                    "catenation degree",
                    "topological genome",
                    "Runtime",
                    "ltol=stol=0.001",
                    "angle_tol=0.01",
                    "supercell_size=num_sites",
                    "no ignored species",
                ):
                    self.assertIn(required, text)

    def test_main_readme_defines_optional_method_keys_where_it_exposes_them(self):
        text = (ROOT / "README.md").read_text(encoding="utf-8")
        start = text.index("`rac5_topology`")
        context = " ".join(text[start : start + 2500].split())
        for required in (
            "(`RT-`)",
            "all 264 finite RAC5 values",
            "network/count/net/agreement",
            "`mofid_v2_topology` (`M2T-`)",
            "collapse each Unicode-whitespace run to one ASCII space, trim",
            "does not modify a CIF",
            "`structure_matcher_strict` (`SM-`)",
            "forward and reverse",
            "ltol=stol=0.001",
            "angle_tol=0.01",
            "supercell_size=num_sites",
            "no ignored species",
            "Missing or failed input adds no optional edge",
            "not a topology, MOFid, RMSD, or score",
            "do not change",
        ):
            self.assertIn(required, context)


if __name__ == "__main__":
    unittest.main()
