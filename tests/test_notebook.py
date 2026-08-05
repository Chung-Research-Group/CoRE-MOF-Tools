import contextlib
import csv
import hashlib
import io
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tests.test_dataset_labels import _make_release


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
        self.assertGreaterEqual(len(code_cells), 16)
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
            "TargetSource",
            "target_names",
            "required_targets",
            "leakage_blocks_use_full_release_universe",
            "doctor()",
            "merge_targets_from_config",
            "AliasRegistry",
            "CURRENT_FEATURE_TABLES",
            '"alias_columns": ["earlier_id"]',
            '"feature_tables": list(available_feature_tables)',
            "filter_precedes_assignment",
            "COREMOF_ROSENBLUTH_TARGET_CONFIG",
            "co2_rosenbluth_weight",
            "n2_rosenbluth_weight",
            "dimensionless Rosenbluth weights",
            "rac5_topology",
            "mofid_v2_topology",
            "structure_matcher_strict",
            "direct symmetric pair edges",
            "historical relaxed matcher",
            "parent_only",
            "screen_candidates.py",
            "filter_precedes_ranking",
            "cell_volume_A3",
            '"--required-target-mode", "all"',
            '"--limit", "100"',
            '"--parent-method", "priority_main"',
            '"--leakage-guard", "auto"',
            "descending Rosenbluth weight is a workflow demonstration",
            "not uptake/selectivity performance",
            "2026-08-05",
            "42,574 release rows",
            "17,245 rows complete",
            "137,960 target observations",
            "308 feature columns",
            "3,828 train / 479 validation / 479 test",
            "zero crossed leakage blocks",
            '"release_structure_count"',
            '"feature_table_count"',
            '"feature_column_count"',
            '"alias_registry_sha256"',
            '"config_sha256"',
            '"target_values_sha256"',
            '"split_counts"',
        ):
            self.assertIn(required, text)

    def test_project_defined_split_terms_are_explained_before_code_use(self):
        cells = self.notebook["cells"]
        definitions = [
            (index, cell)
            for index, cell in enumerate(cells)
            if cell.get("id") == "project-defined-split-identifiers"
        ]
        self.assertEqual(len(definitions), 1)
        definition_index, definition_cell = definitions[0]
        definition = "".join(definition_cell["source"])
        for required in (
            "CoREMOF-tools API identifiers",
            "not community-standard crystallographic terms",
            "RAC5 group anchors a component",
            "MOFid-v2 groups",
            "MOFid-v1 groups",
            "PARENT_METHOD_CONFLICT",
            "not row-wise first-nonmissing",
            "rac_status/rac_group",
            "mofid2_status/mofid2_group",
            "mofid1_status/mofid1_group",
            "five exact relations",
            "full CIF SHA-256",
            "database-namespaced source siblings",
            'leakage_guard=\"auto\"',
            "chooses `main_union`",
            "missing full hash fails closed",
            "adds no cross-method edge",
        ):
            self.assertIn(required, definition)

        first_code_use = min(
            index
            for index, cell in enumerate(cells)
            if cell["cell_type"] == "code"
            and any(
                token in "".join(cell["source"])
                for token in ("priority_main", "main_union")
            )
        )
        self.assertLess(definition_index, first_code_use)

    def test_notebook_defines_identity_canonicalization_and_prefixes_locally(self):
        cells = {cell["id"]: "".join(cell["source"]) for cell in self.notebook["cells"]}
        definition = cells["parent-hierarchy"]
        for required in (
            "canonicalized identifier text",
            "collapse each Unicode-whitespace run",
            "Unicode NFKC",
            "case-fold",
            "(source_database, source_id)",
            "never a prefix or fuzzy match",
            "_ASR_pacman",
            "MOFidv2.",
            ".no_ref",
            "never alter atoms",
            "`identity_union`",
            "provisional source-ID/MOFid transitive groups",
            "identity_status",
            "identity_size",
            "transitive connected component",
            "not a count of edges or identifiers",
            "v26.0.1 components",
            "no precedence or conflict rule",
            "Missing identifiers add no edge",
            "uses no RAC5",
            "not proof of structural identity",
            "`RT-`",
            "`M2T-`",
            "`SM-`",
            "Direct edges are authoritative",
            "not an all-pairs or duplicate-identity claim",
            "all 264 finite ordered descriptor columns",
            "selected Zeo++ fingerprint",
            "current CrystalNets scientific fingerprint",
            "catenation degree",
            "topology key/name",
            "ltol=stol=0.001",
            "angle_tol=0.01",
            "supercell_size=num_sites",
            "no ignored species",
        ):
            self.assertIn(required, definition)

    def test_receipt_example_exposes_machine_readable_split_definitions(self):
        cells = {cell["id"]: "".join(cell["source"]) for cell in self.notebook["cells"]}
        receipt_cell = cells["receipt-summary"]
        for required in (
            'receipt["parent_method_definition"]',
            'receipt["requested_leakage_guard"]',
            'receipt["requested_leakage_guard_definition"]',
            'receipt["leakage_guard_definition"]',
            '== "auto"',
            '== "main_union"',
        ):
            self.assertIn(required, receipt_cell)

    def test_notebook_screening_commands_follow_the_frozen_example_contract(self):
        cells = {cell["id"]: "".join(cell["source"]) for cell in self.notebook["cells"]}
        self.assertIn("screening-cli-contract", cells)
        self.assertIn("screen-release-metadata-command", cells)
        self.assertIn("screen-target-command", cells)
        self.assertIn("screen-rosenbluth-command", cells)

        metadata_command = cells["screen-release-metadata-command"]
        for value in (
            "cell_volume_A3",
            '"--checkers", "5checker"',
            '"--label", "CR"',
            '"--source", "COD"',
        ):
            self.assertIn(value, metadata_command)

        target_command = cells["screen-target-command"]
        for value in (
            "target_config_path",
            '"--rank-by", "xe_uptake"',
            '"--require-target", "xe_uptake"',
            '"--require-target", "xe_kr_selectivity"',
            '"--variant", "ASR"',
            '"--metal", "Cu"',
            '"--split"',
            '"--parent-method", "priority_main"',
            '"--leakage-guard", "auto"',
        ):
            self.assertIn(value, target_command)

        rosenbluth_command = cells["screen-rosenbluth-command"]
        self.assertIn("rosenbluth_config_path.is_file()", rosenbluth_command)
        self.assertIn('"co2_rosenbluth_weight"', rosenbluth_command)
        self.assertIn('"n2_rosenbluth_weight"', rosenbluth_command)

    def test_notebook_base_workflow_uses_no_third_party_notebook_helpers(self):
        code = "\n".join(
            "".join(cell["source"])
            for cell in self.notebook["cells"]
            if cell["cell_type"] == "code"
        )
        for forbidden in ("import pandas", "import numpy", "import nbformat", "import jupyter"):
            self.assertNotIn(forbidden, code)

    def test_notebook_executes_portable_workflow_on_a_small_release(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            release = root / "coremof_vtest"
            output = root / "outputs"
            release.mkdir()
            _make_release(release)
            environment = {
                "COREMOF_RELEASE": str(release),
                "COREMOF_V2602_RELEASE": str(release),
                "COREMOF_NOTEBOOK_OUTPUT": str(output),
            }
            namespace = {"__name__": "__main__"}
            with patch.dict(os.environ, environment, clear=False):
                os.environ.pop("COREMOF_ROSENBLUTH_TARGET_CONFIG", None)
                with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                    for index, cell in enumerate(self.notebook["cells"]):
                        if cell["cell_type"] != "code":
                            continue
                        source = "".join(cell["source"])
                        exec(
                            compile(source, "{}#cell-{}".format(NOTEBOOK, index), "exec"),
                            namespace,
                        )
            self.assertEqual(len(namespace["dataset"]), 4)
            self.assertEqual(
                set(namespace["target_dataset"].target_columns),
                {"xe_uptake", "xe_kr_selectivity"},
            )
            self.assertTrue(namespace["target_split"].leakage_audit["passed"])
            self.assertTrue(namespace["split"].leakage_audit["passed"])
            self.assertIsNone(namespace["real_rosenbluth_dataset"])
            self.assertIsNone(namespace["real_rosenbluth_summary"])
            self.assertIsNone(namespace["rosenbluth_screen_command"])
            self.assertIn("cell_volume_A3", namespace["metadata_screen_command"])
            self.assertIn(str(namespace["target_config_path"]), namespace["target_screen_command"])
            self.assertIn("xe_uptake", namespace["target_screen_command"])
            self.assertIn("xe_kr_selectivity", namespace["target_screen_command"])
            self.assertIn("--split", namespace["target_screen_command"])
            self.assertTrue((output / "example_target_config.json").is_file())
            self.assertTrue((output / "portable_target_features.csv").is_file())
            self.assertTrue((output / "cod_si_5checker_seed42.json").is_file())

    def test_guarded_rosenbluth_workflow_prints_a_path_sanitized_receipt_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            release = root / "coremof_vtest"
            output = root / "outputs"
            release.mkdir()
            _make_release(release)

            with (release / "metadata" / "metadata.csv").open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                structure_ids = [row["structure_id"] for row in csv.DictReader(handle)]

            targets_path = root / "rosenbluth_targets.csv"
            with targets_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=(
                        "structure_id",
                        "co2_rosenbluth_weight",
                        "n2_rosenbluth_weight",
                    ),
                )
                writer.writeheader()
                for index, structure_id in enumerate(structure_ids, start=1):
                    writer.writerow(
                        {
                            "structure_id": structure_id,
                            "co2_rosenbluth_weight": index * 10.0,
                            "n2_rosenbluth_weight": index * 2.0,
                        }
                    )

            aliases_path = root / "structure_aliases.csv"
            with aliases_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=("structure_id", "earlier_id")
                )
                writer.writeheader()
                for structure_id in structure_ids:
                    writer.writerow({"structure_id": structure_id, "earlier_id": ""})

            config_path = root / "rosenbluth_config.json"
            config_path.write_text(
                json.dumps(
                    {
                        "sources": [
                            {
                                "path": targets_path.name,
                                "name": "synthetic_rosenbluth",
                                "target_columns": [
                                    "co2_rosenbluth_weight",
                                    "n2_rosenbluth_weight",
                                ],
                                "value_types": {
                                    "co2_rosenbluth_weight": "float",
                                    "n2_rosenbluth_weight": "float",
                                },
                                "units": {
                                    "co2_rosenbluth_weight": "dimensionless",
                                    "n2_rosenbluth_weight": "dimensionless",
                                },
                                "conditions": {
                                    "co2_rosenbluth_weight": {"temperature_K": 298},
                                    "n2_rosenbluth_weight": {"temperature_K": 298},
                                },
                            }
                        ],
                        "alias_registry": {
                            "path": aliases_path.name,
                            "current_id_column": "structure_id",
                            "alias_columns": ["earlier_id"],
                        },
                        "feature_tables": [],
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n",
                encoding="utf-8",
            )

            environment = {
                "COREMOF_RELEASE": str(release),
                "COREMOF_V2602_RELEASE": str(release),
                "COREMOF_NOTEBOOK_OUTPUT": str(output),
                "COREMOF_ROSENBLUTH_TARGET_CONFIG": str(config_path),
            }
            namespace = {"__name__": "__main__"}
            with patch.dict(os.environ, environment, clear=False):
                with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                    for index, cell in enumerate(self.notebook["cells"]):
                        if cell["cell_type"] != "code":
                            continue
                        source = "".join(cell["source"])
                        exec(
                            compile(source, "{}#cell-{}".format(NOTEBOOK, index), "exec"),
                            namespace,
                        )

            summary = namespace["real_rosenbluth_summary"]
            self.assertEqual(summary["release_structure_count"], 4)
            self.assertEqual(summary["counts"]["observations"], 8)
            self.assertEqual(
                set(summary["target_columns"]),
                {"co2_rosenbluth_weight", "n2_rosenbluth_weight"},
            )
            self.assertEqual(summary["feature_table_count"], 0)
            self.assertEqual(summary["feature_column_count"], 0)
            self.assertEqual(
                summary["sources"],
                [
                    {
                        "file_name": targets_path.name,
                        "sha256": hashlib.sha256(targets_path.read_bytes()).hexdigest(),
                    }
                ],
            )
            self.assertEqual(
                summary["alias_registry_sha256"],
                hashlib.sha256(aliases_path.read_bytes()).hexdigest(),
            )
            self.assertEqual(
                summary["config_sha256"],
                hashlib.sha256(config_path.read_bytes()).hexdigest(),
            )
            self.assertRegex(summary["target_values_sha256"], r"^[0-9a-f]{64}$")
            self.assertEqual(
                summary["split_counts"],
                dict(namespace["real_rosenbluth_split"].counts),
            )
            self.assertEqual(summary["cross_split_block_count"], 0)
            self.assertNotIn(str(root), json.dumps(summary, sort_keys=True))

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
