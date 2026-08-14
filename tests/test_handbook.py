"""Static regression checks for the dataset-splitting README handbook."""

import re
from pathlib import Path
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
HANDBOOK = ROOT / "README_DATASET_SPLITTING.md"


def _heading_slug(value):
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9 _-]", "", value)
    return re.sub(r"[ _]+", "-", value)


class HandbookTests(unittest.TestCase):
    def test_package_source_and_user_surfaces_pass_first_use_audit_normal_and_S(self):
        auditor = (
            ROOT.parent
            / "CoREMOF-COD"
            / "dataset_split"
            / "scripts"
            / "audit_first_use_terminology.py"
        )
        if not auditor.is_file():
            self.skipTest("CoREMOF-COD terminology auditor is not available")
        paths = (
            ROOT / "CoREMOF" / "parents.py",
            ROOT / "CoREMOF" / "splitters.py",
            ROOT / "CoREMOF" / "cli.py",
            ROOT / "README.md",
            HANDBOOK,
            ROOT / "docs" / "source" / "splitting.rst",
            ROOT / "examples" / "README.md",
            ROOT / "examples" / "screen_candidates.py",
            ROOT / "examples" / "CoREMOF_dataset_splitting_quickstart.ipynb",
        )
        for no_site in (False, True):
            command = [sys.executable]
            if no_site:
                command.append("-S")
            command.extend(["-B", str(auditor)])
            command.extend(str(path) for path in paths)
            completed = subprocess.run(
                command,
                cwd=str(ROOT),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )
            self.assertEqual(
                completed.returncode,
                0,
                msg="{}\n{}".format(" ".join(command), completed.stdout),
            )

    def test_docs_extra_matches_furo_configuration_and_scratch_log_is_absent(self):
        setup_text = (ROOT / "setup.py").read_text(encoding="utf-8")
        requirements = (ROOT / "docs" / "requirements.txt").read_text(
            encoding="utf-8"
        )
        configuration = (ROOT / "docs" / "source" / "conf.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("'furo'", setup_text)
        self.assertNotIn("pydata-sphinx-theme", setup_text)
        self.assertIn("furo==", requirements)
        self.assertIn('html_theme = "furo"', configuration)
        self.assertFalse((ROOT / "docs" / "api.log").exists())

    def test_distributable_artifacts_exclude_supporting_information_archives(self):
        setup_text = (ROOT / "setup.py").read_text(encoding="utf-8")
        manifest_text = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
        self.assertNotIn("'data/SI/*.zip',", setup_text)
        self.assertIn(
            "exclude_package_data={'CoREMOF': ['data/SI/*.zip']}", setup_text
        )
        self.assertEqual(
            manifest_text.strip(), "exclude CoREMOF/data/SI/*.zip"
        )

    def test_public_python_container_contract_is_documented(self):
        handbook = HANDBOOK.read_text(encoding="utf-8")
        sphinx = (ROOT / "docs" / "source" / "splitting.rst").read_text(
            encoding="utf-8"
        )
        for text in (handbook, sphinx):
            self.assertIn("ordered", text)
            self.assertIn("list", text)
            self.assertIn("tuple", text)
            self.assertIn("Numeric text", text)
            numeric_text = text[text.index("Numeric text"):][:100]
            self.assertRegex(numeric_text, r"(?:rejected|not coerced)")

    def test_public_docs_use_only_canonical_crystalnets_method_names(self):
        paths = (
            ROOT / "README.md",
            HANDBOOK,
            ROOT / "docs" / "source" / "splitting.rst",
            ROOT / "examples" / "README.md",
        )
        for path in paths:
            with self.subTest(path=path.name):
                text = path.read_text(encoding="utf-8")
                self.assertIn("rac5_crystalnets", text)
                self.assertIn("mofid_v2_crystalnets", text)
                self.assertNotIn("rac5_topology", text)
                self.assertNotIn("mofid_v2_topology", text)

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

    def test_stage_only_mofid_is_not_described_as_published(self):
        handbook = HANDBOOK.read_text(encoding="utf-8")
        sphinx = (ROOT / "docs" / "source" / "splitting.rst").read_text(
            encoding="utf-8"
        )
        main_readme = (ROOT / "README.md").read_text(encoding="utf-8")
        examples_readme = (ROOT / "examples" / "README.md").read_text(
            encoding="utf-8"
        )
        documents = {
            "handbook": handbook,
            "sphinx": sphinx,
            "main_readme": main_readme,
            "examples_readme": examples_readme,
        }
        for document, text in documents.items():
            with self.subTest(document=document):
                normalized = " ".join(text.split())
                self.assertIn("STAGE_ONLY", normalized)
                self.assertIn("non-published", normalized)
                self.assertIn("release-authorized current", normalized)
                self.assertIn(
                    "does not rank, schedule, or recalculate", normalized
                )
        for text in (handbook, sphinx):
            normalized = " ".join(text.split())
            self.assertIn("publication-eligible MOFid evidence", normalized)
            self.assertIn("separate curation-operation queue", normalized)

    def test_identifier_and_optional_group_terms_are_defined_at_first_use(self):
        documents = (
            HANDBOOK,
            ROOT / "docs" / "source" / "splitting.rst",
        )
        for path in documents:
            with self.subTest(path=path.name):
                text = path.read_text(encoding="utf-8")

                canonical_start = text.index("canonicalized identifier text")
                canonical_context = " ".join(
                    text[canonical_start : canonical_start + 4200].split()
                )
                for required in (
                    "whitespace",
                    "Unicode NFKC",
                    "case-fold",
                    "source_database",
                    "source_id",
                    "complete",
                    "fuzzy",
                    "SUCCESS_TOPOLOGY_UNKNOWN",
                    "freshly",
                    "does not seed",
                    "unresolved-reconciliation",
                    "ambiguous-node",
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
                    "independently",
                    "no earlier component",
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

                m2t_start = text.index("M2T-")
                m2t_context = " ".join(
                    text[m2t_start : m2t_start + 2600].split()
                ).casefold()
                for required in (
                    "SUCCESS",
                    "SUCCESS_TOPOLOGY_UNKNOWN",
                    "SUCCESS_TOPOLOGY_ERROR",
                    "SUCCESS_TOPOLOGY_TIMEOUT",
                    "every other",
                    "adds no edge",
                    "release-authorized mofid-v2 values change",
                    "rebuild",
                ):
                    self.assertIn(required.casefold(), m2t_context)

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
        start = text.index("`rac5_crystalnets`")
        context = " ".join(text[start : start + 3500].split())
        for required in (
            "(`RT-`)",
            "all 264 finite RAC5 values",
            "network/count/net/agreement",
            "`mofid_v2_crystalnets` (`M2T-`)",
            "collapse each Unicode-whitespace run to one ASCII space, trim",
            "does not modify a CIF",
            "`SUCCESS_TOPOLOGY_UNKNOWN`",
            "`SUCCESS_TOPOLOGY_ERROR`",
            "`SUCCESS_TOPOLOGY_TIMEOUT`",
            "every other MOFid-v2 status",
            "release-authorized MOFid-v2 values change",
            "rebuild the M2T groups before use",
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

    def test_examples_readme_states_the_m2t_rebuild_trigger(self):
        text = " ".join(
            (ROOT / "examples" / "README.md")
            .read_text(encoding="utf-8")
            .split()
        )
        self.assertIn("release-authorized MOFid-v2 values change", text)
        self.assertIn("rebuild the M2T groups before use", text)


if __name__ == "__main__":
    unittest.main()
