import csv
import copy
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import pickle
import sys
import subprocess
import tempfile
from types import ModuleType
import unittest
from unittest.mock import patch

import CoREMOF.dataset as dataset_module
from CoREMOF.dataset import CoREMOFDataset, ReleaseValidationError
from CoREMOF.labels import (
    CHECKER_PRESETS,
    UnknownCheckerStatusError,
    consensus_label,
    resolve_checker_preset,
    resolve_checker_view,
)


CHECKER_COLUMNS = (
    "mofclassifier_status",
    "mofchecker_status",
    "chen_manz_status",
    "mosaec_status",
    "setc_gat_status",
)


def _write_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _metadata_rows():
    base = [
        (
            "ASR-COD-2026-0001",
            "COD",
            "ASR",
            "Cu",
            ("PASS", "PASS", "PASS", "PASS", "PASS"),
            ("CR", "CR", "CR"),
        ),
        (
            "FSR-COD-2026-0001",
            "COD",
            "FSR",
            "Fe;Zn",
            ("FAIL", "FAIL", "FAIL", "FAIL", "FAIL"),
            ("NCR", "NCR", "NCR"),
        ),
        (
            "ASR-SI-2025-0001",
            "SI",
            "ASR",
            "Zn",
            ("PASS", "FAIL", "PASS", "FAIL", "PASS"),
            ("AMBIGUOUS", "AMBIGUOUS", "AMBIGUOUS"),
        ),
        (
            "ION-CSD-2024-0001",
            "CSD",
            "ION",
            "Co",
            ("PASS", "PASS", "NOT_AVAILABLE", "PASS", "PASS"),
            ("UNCHECKED", "UNCHECKED", "UNCHECKED"),
        ),
    ]
    rows = []
    for structure_id, source, variant, metals, statuses, labels in base:
        row = {
            "structure_id": structure_id,
            "cif_file": "cifs/{}.cif".format(structure_id),
            "source_database": source,
            "source_id": "SOURCE-{}".format(structure_id),
            "structure_variant": variant,
            "metal_elements": metals,
            "mofid_v2": "MOF-SHARED",
            "mofid_v2_status": "SUCCESS",
            "label_3checker": labels[0],
            "label_4checker": labels[1],
            "label_5checker": labels[2],
        }
        row.update(dict(zip(CHECKER_COLUMNS, statuses)))
        rows.append(row)
    return rows


def _parent_rows():
    ids = [row["structure_id"] for row in _metadata_rows()]
    rows = []
    for index, structure_id in enumerate(ids):
        if index < 2:
            rac = ("MATCHED", "R-AAAA0001", "2")
        elif index == 2:
            rac = ("UNMATCHED", "R-AAAA0002", "1")
        else:
            rac = ("NOT_AVAILABLE", "R-AAAA0003", "1")
        row = {"structure_id": structure_id}
        for prefix, values in (
            ("rac_zeo", ("UNMATCHED", "RZ-{:08d}".format(index), "1")),
            ("rac", rac),
            ("zeo", ("UNMATCHED", "Z-{:08d}".format(index), "1")),
            ("source", ("UNMATCHED", "S-{:08d}".format(index), "1")),
            ("mofid2", ("UNMATCHED", "M2-{:08d}".format(index), "1")),
            ("mofid1", ("UNMATCHED", "M1-{:08d}".format(index), "1")),
            ("name", ("UNMATCHED", "N-{:08d}".format(index), "1")),
            ("identity", ("UNMATCHED", "I-{:08d}".format(index), "1")),
        ):
            row["{}_status".format(prefix)] = values[0]
            row["{}_group".format(prefix)] = values[1]
            row["{}_size".format(prefix)] = values[2]
        rows.append(row)
    return rows


def _make_release(root):
    metadata = _metadata_rows()
    parent_rows = _parent_rows()
    metadata_fields = tuple(metadata[0])
    parent_fields = tuple(parent_rows[0])
    _write_csv(root / "metadata" / "metadata.csv", metadata_fields, metadata)
    _write_csv(
        root / "parent_groups" / "parent_groups.csv", parent_fields, parent_rows
    )

    info = {
        "dataset_version": "vtest",
        "structure_count": len(metadata),
        "classification_definitions": {
            name: list(checkers) for name, checkers in CHECKER_PRESETS.items()
        },
    }
    (root / "dataset_info.json").write_text(
        json.dumps(info), encoding="utf-8"
    )
    methods = {
        "dataset_version": "vtest",
        "csv_column_prefixes": {
            "rac5_zeo": "rac_zeo",
            "rac5": "rac",
            "zeo": "zeo",
            "source_id": "source",
            "mofid_v2": "mofid2",
            "mofid_v1": "mofid1",
            "common_name": "name",
            "identity_union": "identity",
        },
    }
    (root / "parent_groups" / "parent_group_methods.json").write_text(
        json.dumps(methods), encoding="utf-8"
    )

    manifest = [
        {
            "structure_id": row["structure_id"],
            "cif_file": row["cif_file"],
            "size_bytes": "10",
            "sha256": hashlib.sha256(
                row["structure_id"].encode("utf-8")
            ).hexdigest(),
        }
        for row in metadata
    ]
    _write_csv(
        root / "manifests" / "cif_manifest.csv",
        ("structure_id", "cif_file", "size_bytes", "sha256"),
        manifest,
    )


def _add_structure_matcher_parent_columns(root):
    parents_path = root / "parent_groups" / "parent_groups.csv"
    with parents_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    for index, row in enumerate(rows):
        row["sm_status"] = "MATCHED" if index < 2 else "UNMATCHED"
        row["sm_group"] = (
            "SM-ABCD0001" if index < 2 else "SM-{:08X}".format(index)
        )
        row["sm_size"] = "2" if index < 2 else "1"
    _write_csv(parents_path, tuple(rows[0]), rows)
    return rows


def _add_crystalnets_parent_columns(root):
    parents_path = root / "parent_groups" / "parent_groups.csv"
    with parents_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    for index, row in enumerate(rows):
        for prefix, group_prefix in (
            ("rac_crystalnets", "RT-"),
            ("mofid2_crystalnets", "M2T-"),
        ):
            if index < 2:
                row["{}_status".format(prefix)] = "MATCHED"
                row["{}_group".format(prefix)] = group_prefix + "ABCD0001"
                row["{}_size".format(prefix)] = "2"
            else:
                row["{}_status".format(prefix)] = "NOT_AVAILABLE"
                row["{}_group".format(prefix)] = group_prefix + "{:08X}".format(index)
                row["{}_size".format(prefix)] = "1"
    _write_csv(parents_path, tuple(rows[0]), rows)
    return rows


def _declare_structure_matcher_contract(root):
    parents_path = root / "parent_groups" / "parent_groups.csv"
    with parents_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    receipt = {
        "schema_version": "coremof-structure-matcher-release-adapter-receipt/1.0",
        "status": "PASS",
        "dataset_version": "vtest",
        "method_id": "pymatgen_structure_matcher_strict_v2",
        "method_schema_version": "coremof-structure-matcher-method/2.0",
        "parent_groups_sha256": hashlib.sha256(parents_path.read_bytes()).hexdigest(),
        "strict_pair_ledger_sha256": hashlib.sha256(b"strict-pair-ledger").hexdigest(),
        "structure_count": len(rows),
        "candidate_pair_count": 1,
        "successful_pair_count": 1,
        "unresolved_pair_count": 0,
        "strict_direct_match_edge_count": 1,
        "not_available_structure_count": sum(
            row["sm_status"] == "NOT_AVAILABLE" for row in rows
        ),
        "historical_relaxed_executed": False,
        "historical_relaxed_exposed": False,
    }
    receipt_path = (
        root
        / "parent_groups"
        / "structure_matcher_strict_evidence_receipt.json"
    )
    receipt_path.write_text(json.dumps(receipt, sort_keys=True), encoding="utf-8")

    methods_path = root / "parent_groups" / "parent_group_methods.json"
    methods = json.loads(methods_path.read_text(encoding="utf-8"))
    methods["csv_column_prefixes"]["structure_matcher_strict"] = "sm"
    methods.setdefault("criteria", {})["structure_matcher_strict"] = {
        "role": "OPTIONAL_REFERENCE",
        "method_id": "pymatgen_structure_matcher_strict_v2",
        "method_schema_version": "coremof-structure-matcher-method/2.0",
        "authoritative_evidence": "DIRECT_SYMMETRIC_STRICT_MATCH_EDGES",
        "fit_policy": "FIT_SYMMETRIC_TRUE_REQUIRED",
        "component_semantics": (
            "CONNECTED_COMPONENT_CONVENIENCE_DIRECT_EDGES_AUTHORITATIVE"
        ),
        "component_completeness_policy": (
            "INCOMPLETE_COMPONENTS_NOT_AVAILABLE_UNIQUE_SINGLETON"
        ),
        "public_status_policy": ["MATCHED", "UNMATCHED", "NOT_AVAILABLE"],
        "included_in_priority_main": False,
        "included_in_main_union": False,
        "historical_relaxed_executed": False,
        "historical_relaxed_exposed": False,
        "evidence_receipt": {
            "file": "parent_groups/structure_matcher_strict_evidence_receipt.json",
            "sha256": hashlib.sha256(receipt_path.read_bytes()).hexdigest(),
        },
    }
    methods_path.write_text(json.dumps(methods, sort_keys=True), encoding="utf-8")
    return methods_path, receipt_path


def _declare_crystalnets_reference_contracts(root):
    methods_path = root / "parent_groups" / "parent_group_methods.json"
    methods = json.loads(methods_path.read_text(encoding="utf-8"))
    methods["csv_column_prefixes"].update(
        {
            "rac5_crystalnets": "rac_crystalnets",
            "mofid_v2_crystalnets": "mofid2_crystalnets",
        }
    )
    criteria = methods.setdefault("criteria", {})
    definitions = methods.setdefault("definitions", {})
    for method in ("rac5_crystalnets", "mofid_v2_crystalnets"):
        contract = json.loads(
            json.dumps(dataset_module._CRYSTALNETS_REFERENCE_CONTRACTS[method])
        )
        criteria[method] = contract
        definitions[method] = json.loads(json.dumps(contract))
    methods["crystalnets_reference_integration"] = json.loads(
        json.dumps(dict(dataset_module._CRYSTALNETS_REFERENCE_INTEGRATION))
    )
    methods_path.write_text(json.dumps(methods, sort_keys=True), encoding="utf-8")
    info_path = root / "dataset_info.json"
    info = json.loads(info_path.read_text(encoding="utf-8"))
    info_definitions = info.setdefault("definitions", {})
    parent_grouping = info.setdefault("parent_grouping", {})
    info_contracts = parent_grouping.setdefault("criterion_contracts", {})
    for method in ("rac5_crystalnets", "mofid_v2_crystalnets"):
        contract = json.loads(
            json.dumps(dataset_module._CRYSTALNETS_REFERENCE_CONTRACTS[method])
        )
        info_definitions[method] = contract
        info_contracts[method] = json.loads(json.dumps(contract))
    parent_grouping["crystalnets_reference_integration"] = json.loads(
        json.dumps(dict(dataset_module._CRYSTALNETS_REFERENCE_INTEGRATION))
    )
    info_path.write_text(json.dumps(info, sort_keys=True), encoding="utf-8")
    return methods_path


class LabelTests(unittest.TestCase):
    def test_strict_consensus_contract(self):
        self.assertEqual(consensus_label(["PASS", "PASS"]), "CR")
        self.assertEqual(consensus_label(["FAIL", "FAIL"]), "NCR")
        self.assertEqual(consensus_label(["PASS", "FAIL"]), "AMBIGUOUS")
        for status in ("NOT_AVAILABLE", "ERROR", "TIMEOUT", "PENDING"):
            self.assertEqual(consensus_label(["PASS", status]), "UNCHECKED")

    def test_unknown_status_fails_closed(self):
        for status in ("pass", "BROKEN", ""):
            with self.subTest(status=status):
                with self.assertRaises(UnknownCheckerStatusError):
                    consensus_label([status])
        with self.assertRaises(UnknownCheckerStatusError):
            consensus_label([None])

    def test_only_official_presets_are_accepted_as_presets(self):
        self.assertEqual(resolve_checker_preset(3)[0], "3checker")
        self.assertEqual(resolve_checker_preset("5checker")[0], "5checker")
        with self.assertRaises(ValueError):
            resolve_checker_preset("2checker")

    def test_custom_view_is_explicit_and_ordered(self):
        name, checkers, official = resolve_checker_view(
            ["MOSAEC", "MOFClassifier"]
        )
        self.assertEqual(name, "custom:MOSAEC+MOFClassifier")
        self.assertEqual(checkers, ("MOSAEC", "MOFClassifier"))
        self.assertFalse(official)
        with self.assertRaises(ValueError):
            resolve_checker_view(["MOSAEC", "MOSAEC"])
        with self.assertRaises(TypeError):
            resolve_checker_view({"MOSAEC", "MOFClassifier"})


class DatasetTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name, "release")
        self.root.mkdir()
        _make_release(self.root)

    def tearDown(self):
        self.temporary_directory.cleanup()

    def test_load_join_and_recompute_published_labels(self):
        dataset = CoREMOFDataset.from_release(self.root)
        self.assertEqual(dataset.dataset_version, "vtest")
        self.assertEqual(len(dataset), 4)
        self.assertEqual(dataset.structure_ids[0], "ASR-COD-2026-0001")
        self.assertEqual(len(dataset.metadata_rows), 4)
        self.assertEqual(dataset.parent_by_id[dataset.structure_ids[0]]["rac_size"], "2")
        self.assertEqual(
            dataset.cif_hashes[dataset.structure_ids[0]],
            hashlib.sha256(dataset.structure_ids[0].encode("utf-8")).hexdigest(),
        )
        self.assertIn("metadata/metadata.csv", dataset.input_hashes)
        with self.assertRaises(TypeError):
            dataset.dataset_info["classification_definitions"]["3checker"] = []
        first = dataset["ASR-COD-2026-0001"]
        self.assertEqual(first.parent_group("rac5").group_id, "R-AAAA0001")
        self.assertEqual(first.cif_manifest["size_bytes"], "10")

        classified = dataset.classify(3)
        self.assertEqual(
            classified.labels, ("CR", "NCR", "AMBIGUOUS", "UNCHECKED")
        )
        self.assertEqual(classified.label_counts()["CR"], 1)
        self.assertEqual(classified.cr_ids, ("ASR-COD-2026-0001",))
        self.assertEqual(classified.ncr_ids, ("FSR-COD-2026-0001",))

    def test_authenticated_generations_reject_copy_and_mutation_attacks(self):
        dataset = CoREMOFDataset.from_release(self.root)
        classified = dataset.classify("5checker")

        for value in (dataset, classified):
            for operation in (copy.copy, copy.deepcopy, pickle.dumps):
                with self.subTest(value=type(value).__name__, operation=operation.__name__):
                    with self.assertRaisesRegex(TypeError, "cannot be"):
                        operation(value)

        with self.assertRaisesRegex(AttributeError, "immutable"):
            dataset.records = tuple(reversed(dataset.records))
        with self.assertRaisesRegex(AttributeError, "immutable"):
            classified.label_by_id = {}

        object.__setattr__(dataset, "records", tuple(reversed(dataset.records)))
        with self.assertRaisesRegex(
            ValueError, "generation changed|fingerprint changed"
        ):
            dataset.classify("5checker")

        dataset = CoREMOFDataset.from_release(self.root)
        classified = dataset.classify("5checker")
        forged_record = replace(classified.records[0], label="NCR")
        object.__setattr__(
            classified,
            "records",
            (forged_record,) + classified.records[1:],
        )
        with self.assertRaisesRegex(ValueError, "checker-view generation changed"):
            classified.filter(labels=("CR",))

    def test_authenticated_generation_detects_nested_dataclass_bypass(self):
        dataset = CoREMOFDataset.from_release(self.root)
        original = dataset.records[0].structure_id
        object.__setattr__(dataset.records[0], "structure_id", original + "-FORGED")
        with self.assertRaisesRegex(
            ValueError, "generation changed|fingerprint changed"
        ):
            dataset.classify("5checker")

        dataset = CoREMOFDataset.from_release(self.root)
        classified = dataset.classify("5checker")
        object.__setattr__(classified.records[0], "label", "NCR")
        with self.assertRaisesRegex(
            ValueError, "checker-view generation changed|fingerprint changed"
        ):
            classified.filter(labels=("CR",))
        self.assertEqual(classified.ids_for_label("ambiguous"), ("ASR-SI-2025-0001",))
        self.assertEqual(classified.unchecked_ids, ("ION-CSD-2024-0001",))
        with self.assertRaises(ValueError):
            classified.ids_for_label("MAYBE")

    def test_optional_crystalnets_combined_parent_columns_are_validated_and_loaded(self):
        _add_crystalnets_parent_columns(self.root)
        _declare_crystalnets_reference_contracts(self.root)

        dataset = CoREMOFDataset.from_release(self.root)
        first = dataset[dataset.structure_ids[0]]
        second = dataset[dataset.structure_ids[1]]
        self.assertEqual(
            first.parent_group("rac5_crystalnets").group_id,
            second.parent_group("rac_crystalnets").group_id,
        )
        self.assertEqual(
            first.parent_group("mofid_v2_crystalnets").group_id,
            "M2T-ABCD0001",
        )

    def test_m2t_available_rows_require_explicit_eligible_mofid_v2_status(self):
        metadata_path = self.root / "metadata" / "metadata.csv"
        for value in (
            "ERROR",
            "TIMEOUT",
            "NOT_AVAILABLE",
            "success",
            " SUCCESS",
            "",
        ):
            with self.subTest(ineligible=value):
                _make_release(self.root)
                _add_crystalnets_parent_columns(self.root)
                _declare_crystalnets_reference_contracts(self.root)
                rows = _metadata_rows()
                rows[0]["mofid_v2_status"] = value
                _write_csv(metadata_path, tuple(rows[0]), rows)
                with self.assertRaisesRegex(
                    ReleaseValidationError, "ineligible mofid_v2_status"
                ):
                    CoREMOFDataset.from_release(self.root)

        for value in sorted(dataset_module.M2T_ELIGIBLE_MOFID_V2_STATUSES):
            with self.subTest(eligible=value):
                _make_release(self.root)
                _add_crystalnets_parent_columns(self.root)
                _declare_crystalnets_reference_contracts(self.root)
                rows = _metadata_rows()
                rows[0]["mofid_v2_status"] = value
                _write_csv(metadata_path, tuple(rows[0]), rows)
                CoREMOFDataset.from_release(self.root)

    def test_m2t_available_rows_require_one_complete_canonical_mofid_v2_per_group(self):
        metadata_path = self.root / "metadata" / "metadata.csv"
        for value in (None, "", "   ", "NOT_AVAILABLE", "error"):
            with self.subTest(invalid=value):
                _make_release(self.root)
                _add_crystalnets_parent_columns(self.root)
                _declare_crystalnets_reference_contracts(self.root)
                rows = _metadata_rows()
                if value is None:
                    for row in rows:
                        row.pop("mofid_v2", None)
                else:
                    rows[0]["mofid_v2"] = value
                _write_csv(metadata_path, tuple(rows[0]), rows)
                with self.assertRaisesRegex(
                    ReleaseValidationError,
                    "exact string|complete nonplaceholder mofid_v2",
                ):
                    CoREMOFDataset.from_release(self.root)

        _make_release(self.root)
        _add_crystalnets_parent_columns(self.root)
        _declare_crystalnets_reference_contracts(self.root)
        rows = _metadata_rows()
        rows[0]["mofid_v2"] = "MOF-A"
        rows[1]["mofid_v2"] = "MOF-B"
        _write_csv(metadata_path, tuple(rows[0]), rows)
        with self.assertRaisesRegex(
            ReleaseValidationError, "conflicting canonical mofid_v2"
        ):
            CoREMOFDataset.from_release(self.root)

        _make_release(self.root)
        _add_crystalnets_parent_columns(self.root)
        _declare_crystalnets_reference_contracts(self.root)
        rows = _metadata_rows()
        rows[0]["mofid_v2"] = "  ＭＯＦ-A  "
        rows[1]["mofid_v2"] = "mof-a"
        _write_csv(metadata_path, tuple(rows[0]), rows)
        CoREMOFDataset.from_release(self.root)

    def test_crystalnets_columns_require_machine_readable_reference_contracts(self):
        _add_crystalnets_parent_columns(self.root)
        methods_path = self.root / "parent_groups" / "parent_group_methods.json"
        methods = json.loads(methods_path.read_text(encoding="utf-8"))
        methods["csv_column_prefixes"].update(
            {
                "rac5_crystalnets": "rac_crystalnets",
                "mofid_v2_crystalnets": "mofid2_crystalnets",
            }
        )
        methods_path.write_text(json.dumps(methods), encoding="utf-8")
        with self.assertRaisesRegex(ReleaseValidationError, "criteria.rac5_crystalnets"):
            CoREMOFDataset.from_release(self.root)

    def test_crystalnets_reference_contract_fails_closed_on_semantic_drift(self):
        _add_crystalnets_parent_columns(self.root)
        methods_path = _declare_crystalnets_reference_contracts(self.root)
        original = json.loads(methods_path.read_text(encoding="utf-8"))

        def mutate_key(value, *path):
            document = json.loads(json.dumps(original))
            target = document
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = value
            methods_path.write_text(json.dumps(document), encoding="utf-8")

        cases = (
            (False, ("criteria", "rac5_crystalnets", "key")),
            ("R-WRONG", ("criteria", "rac5_crystalnets", "group_prefix")),
            (True, ("criteria", "rac5_crystalnets", "decision_contract", "publication_authorized")),
            (True, ("criteria", "rac5_crystalnets", "decision_contract", "priority_main", "consumed")),
            (True, ("criteria", "rac5_crystalnets", "decision_contract", "main_union", "changed")),
            (["mofid_v2", "rac5", "mofid_v1"], ("criteria", "rac5_crystalnets", "decision_contract", "priority_main", "authorized_criterion_order")),
            (["rac5_crystalnets"], ("criteria", "rac5_crystalnets", "decision_contract", "main_union", "edge_relations")),
            ("ERROR is eligible", ("criteria", "rac5_crystalnets", "inputs", "crystalnets", "availability_gate")),
            ("rounded to 3 decimals", ("criteria", "rac5_crystalnets", "inputs", "rac5")),
            ("missing values match each other", ("criteria", "rac5_crystalnets", "missing_error_behavior")),
            ("FINAL_PUBLISHED", ("criteria", "rac5_crystalnets", "release_state")),
            ("one-structure singleton", ("criteria", "rac5_crystalnets", "result_mapping", "MATCHED")),
            ("multi-member", ("criteria", "rac5_crystalnets", "result_mapping", "UNMATCHED")),
            ("eligible match", ("criteria", "rac5_crystalnets", "result_mapping", "NOT_AVAILABLE")),
            (["SUCCESS"], ("criteria", "mofid_v2_crystalnets", "inputs", "mofid_v2", "required_statuses")),
            (["SUCCESS", ["SUCCESS_TOPOLOGY_UNKNOWN"], "SUCCESS_TOPOLOGY_ERROR", "SUCCESS_TOPOLOGY_TIMEOUT"], ("criteria", "mofid_v2_crystalnets", "inputs", "mofid_v2", "required_statuses")),
            ("coordinates are changed", ("criteria", "mofid_v2_crystalnets", "inputs", "mofid_v2", "canonical_text", "scientific_effect")),
            ("ERROR is also eligible", ("criteria", "mofid_v2_crystalnets", "inputs", "mofid_v2", "status_gate")),
            (False, ("criteria", "mofid_v2_crystalnets", "decision_contract", "mofid_v2_crystalnets_rebuild_trigger")),
        )
        for value, path in cases:
            with self.subTest(path=path):
                mutate_key(value, *path)
                with self.assertRaisesRegex(
                    ReleaseValidationError, "exact closed terminal reference contract"
                ):
                    CoREMOFDataset.from_release(self.root)

        document = json.loads(json.dumps(original))
        document["criteria"]["rac5_crystalnets"]["publication_authorized"] = True
        methods_path.write_text(json.dumps(document), encoding="utf-8")
        with self.assertRaisesRegex(
            ReleaseValidationError, "exact closed terminal reference contract"
        ):
            CoREMOFDataset.from_release(self.root)

    def test_dataset_info_crystalnets_contract_copies_fail_closed_on_drift(self):
        _add_crystalnets_parent_columns(self.root)
        _declare_crystalnets_reference_contracts(self.root)
        info_path = self.root / "dataset_info.json"
        original = json.loads(info_path.read_text(encoding="utf-8"))

        cases = (
            (
                "DECISIVE PUBLISHED PARENT",
                ("definitions", "rac5_crystalnets", "purpose"),
                "definitions.rac5_crystalnets",
            ),
            (
                True,
                (
                    "parent_grouping",
                    "criterion_contracts",
                    "mofid_v2_crystalnets",
                    "publication_authorized",
                ),
                "criterion_contracts.mofid_v2_crystalnets",
            ),
            (
                True,
                (
                    "parent_grouping",
                    "crystalnets_reference_integration",
                    "priority_main_changed",
                ),
                "crystalnets_reference_integration",
            ),
            (
                True,
                (
                    "parent_grouping",
                    "crystalnets_reference_integration",
                    "main_union_changed",
                ),
                "crystalnets_reference_integration",
            ),
            (
                True,
                (
                    "parent_grouping",
                    "crystalnets_reference_integration",
                    "publication_authorized",
                ),
                "crystalnets_reference_integration",
            ),
            (
                0,
                (
                    "parent_grouping",
                    "crystalnets_reference_integration",
                    "publication_authorized",
                ),
                "crystalnets_reference_integration",
            ),
            (
                "FINAL_PUBLISHED",
                (
                    "parent_grouping",
                    "crystalnets_reference_integration",
                    "candidate_state",
                ),
                "crystalnets_reference_integration",
            ),
            (
                "LOCAL PATHS ARE PUBLIC",
                (
                    "parent_grouping",
                    "crystalnets_reference_integration",
                    "evidence_state",
                ),
                "crystalnets_reference_integration",
            ),
        )
        for value, path, message in cases:
            with self.subTest(path=path):
                document = json.loads(json.dumps(original))
                target = document
                for key in path[:-1]:
                    target = target[key]
                target[path[-1]] = value
                info_path.write_text(json.dumps(document), encoding="utf-8")
                with self.assertRaisesRegex(ReleaseValidationError, message):
                    CoREMOFDataset.from_release(self.root)

        for path, message in (
            (
                ("definitions", "rac5_crystalnets"),
                "definitions.rac5_crystalnets",
            ),
            (
                (
                    "parent_grouping",
                    "criterion_contracts",
                    "mofid_v2_crystalnets",
                ),
                "criterion_contracts.mofid_v2_crystalnets",
            ),
            (
                ("parent_grouping", "crystalnets_reference_integration"),
                "crystalnets_reference_integration",
            ),
        ):
            with self.subTest(missing=path):
                document = json.loads(json.dumps(original))
                target = document
                for key in path[:-1]:
                    target = target[key]
                del target[path[-1]]
                info_path.write_text(json.dumps(document), encoding="utf-8")
                with self.assertRaisesRegex(ReleaseValidationError, message):
                    CoREMOFDataset.from_release(self.root)

        document = json.loads(json.dumps(original))
        document["parent_grouping"]["crystalnets_reference_integration"][
            "extra_authority"
        ] = True
        info_path.write_text(json.dumps(document), encoding="utf-8")
        with self.assertRaisesRegex(
            ReleaseValidationError, "crystalnets_reference_integration"
        ):
            CoREMOFDataset.from_release(self.root)

    def test_parent_method_crystalnets_contract_copies_fail_closed_on_drift(self):
        _add_crystalnets_parent_columns(self.root)
        methods_path = _declare_crystalnets_reference_contracts(self.root)
        original = json.loads(methods_path.read_text(encoding="utf-8"))

        cases = (
            (
                "DECISIVE PUBLISHED PARENT",
                ("definitions", "rac5_crystalnets", "purpose"),
                "definitions.rac5_crystalnets",
            ),
            (
                True,
                (
                    "definitions",
                    "mofid_v2_crystalnets",
                    "publication_authorized",
                ),
                "definitions.mofid_v2_crystalnets",
            ),
            (
                True,
                ("crystalnets_reference_integration", "priority_main_changed"),
                "crystalnets_reference_integration",
            ),
            (
                True,
                ("crystalnets_reference_integration", "main_union_changed"),
                "crystalnets_reference_integration",
            ),
            (
                True,
                ("crystalnets_reference_integration", "publication_authorized"),
                "crystalnets_reference_integration",
            ),
            (
                0,
                ("crystalnets_reference_integration", "publication_authorized"),
                "crystalnets_reference_integration",
            ),
        )
        for value, path, message in cases:
            with self.subTest(path=path):
                document = json.loads(json.dumps(original))
                target = document
                for key in path[:-1]:
                    target = target[key]
                target[path[-1]] = value
                methods_path.write_text(json.dumps(document), encoding="utf-8")
                with self.assertRaisesRegex(ReleaseValidationError, message):
                    CoREMOFDataset.from_release(self.root)

        for key, message in (
            ("definitions", "definitions"),
            ("crystalnets_reference_integration", "crystalnets_reference_integration"),
        ):
            with self.subTest(missing=key):
                document = json.loads(json.dumps(original))
                del document[key]
                methods_path.write_text(json.dumps(document), encoding="utf-8")
                with self.assertRaisesRegex(ReleaseValidationError, message):
                    CoREMOFDataset.from_release(self.root)

        document = json.loads(json.dumps(original))
        document["crystalnets_reference_integration"]["extra_authority"] = True
        methods_path.write_text(json.dumps(document), encoding="utf-8")
        with self.assertRaisesRegex(
            ReleaseValidationError, "crystalnets_reference_integration"
        ):
            CoREMOFDataset.from_release(self.root)

    def test_old_topology_method_names_are_not_public_release_aliases(self):
        methods_path = self.root / "parent_groups" / "parent_group_methods.json"
        methods = json.loads(methods_path.read_text(encoding="utf-8"))
        methods["csv_column_prefixes"].update(
            {
                "rac5_topology": "rac_topology",
                "mofid_v2_topology": "mofid2_topology",
            }
        )
        methods_path.write_text(json.dumps(methods), encoding="utf-8")
        with self.assertRaisesRegex(
            ReleaseValidationError, "retired or reserved declaration key"
        ):
            CoREMOFDataset.from_release(self.root)

    def test_retired_and_reserved_authority_keys_are_rejected_recursively(self):
        for target, key in (
            ("methods", "rac5_topology"),
            ("methods", "ＲＡＣ５＿Ｔｏｐｏｌｏｇｙ"),
            ("methods", "_authority_generation_marker"),
            ("info", "mofid_v2_topology_status"),
            ("info", "_validated_release_token"),
        ):
            with self.subTest(target=target, key=key):
                _make_release(self.root)
                path = (
                    self.root / "parent_groups" / "parent_group_methods.json"
                    if target == "methods"
                    else self.root / "dataset_info.json"
                )
                document = json.loads(path.read_text(encoding="utf-8"))
                document.setdefault("nested", {}).setdefault("deeper", {})[key] = True
                path.write_text(json.dumps(document), encoding="utf-8")
                with self.assertRaisesRegex(
                    ReleaseValidationError, "retired or reserved declaration key"
                ):
                    CoREMOFDataset.from_release(self.root)

        authority_keys = (
            "publication_authorized",
            "PUBLICATION.AUTHORIZED",
            "Publication-Authorized",
            "Ｐｕｂｌｉｃａｔｉｏｎ＿Ａｕｔｈｏｒｉｚｅｄ",
            "publication_status",
            "Publication-Status",
            "publication_authorization",
            "Publication Authorization",
            "publication_ready",
            "Publication-Ready",
            "publishable",
            "authoritative",
            "official_split",
            "Official-Split",
            "release_authority",
            "Release Authority",
            "release_status",
            "Release-Status",
            "release_state",
            "Release-State",
            "candidate_status",
            "candidate_state",
            "evidence_state",
            "priority_main_changed",
            "main_union_changed",
            "checker_view_official",
            "Checker.View-Official",
        )
        authority_claims = tuple(
            (key, value)
            for key in authority_keys
            for value in (
                True,
                False,
                "FINAL",
                "PROVISIONAL_LATEST_AUDITED_SNAPSHOT",
            )
        )
        for target in ("methods", "info"):
            for key, value in authority_claims:
                with self.subTest(target=target, authority_key=key, value=value):
                    _make_release(self.root)
                    path = (
                        self.root / "parent_groups" / "parent_group_methods.json"
                        if target == "methods"
                        else self.root / "dataset_info.json"
                    )
                    document = json.loads(path.read_text(encoding="utf-8"))
                    document.setdefault("nested", {}).setdefault("deeper", {})[
                        key
                    ] = value
                    path.write_text(json.dumps(document), encoding="utf-8")
                    with self.assertRaisesRegex(
                        ReleaseValidationError,
                        "retired or reserved declaration key",
                    ):
                        CoREMOFDataset.from_release(self.root)

        for key, value in (
            ("rac5_topology", "forged"),
            ("ＲＡＣ５＿Ｔｏｐｏｌｏｇｙ", "forged"),
            ("_official_checker_view_token", "forged"),
            *authority_claims,
        ):
            with self.subTest(target="metadata", key=key, value=value):
                _make_release(self.root)
                metadata_path = self.root / "metadata" / "metadata.csv"
                with metadata_path.open("r", encoding="utf-8", newline="") as handle:
                    reader = csv.DictReader(handle)
                    fields = list(reader.fieldnames or ()) + [key]
                    rows = list(reader)
                rows[0][key] = value
                _write_csv(metadata_path, fields, rows)
                with self.assertRaisesRegex(
                    ReleaseValidationError, "retired or reserved declaration key"
                ):
                    CoREMOFDataset.from_release(self.root)

    def test_only_exact_staged_release_state_paths_are_accepted(self):
        staged_status = "PROVISIONAL_LATEST_AUDITED_SNAPSHOT"
        _add_crystalnets_parent_columns(self.root)
        _declare_crystalnets_reference_contracts(self.root)

        info_path = self.root / "dataset_info.json"
        methods_path = (
            self.root / "parent_groups" / "parent_group_methods.json"
        )
        info = json.loads(info_path.read_text(encoding="utf-8"))
        methods = json.loads(methods_path.read_text(encoding="utf-8"))
        info["release_status"] = staged_status
        info_path.write_text(json.dumps(info, sort_keys=True), encoding="utf-8")
        methods_path.write_text(
            json.dumps(methods, sort_keys=True), encoding="utf-8"
        )
        loaded = CoREMOFDataset.from_release(self.root)
        self.assertEqual(loaded.dataset_info["release_status"], staged_status)

        methods = json.loads(methods_path.read_text(encoding="utf-8"))
        integration = methods["crystalnets_reference_integration"]
        integration["Candidate-State"] = integration.pop("candidate_state")
        methods_path.write_text(
            json.dumps(methods, sort_keys=True), encoding="utf-8"
        )
        with self.assertRaisesRegex(
            ReleaseValidationError, "retired or reserved declaration key"
        ):
            CoREMOFDataset.from_release(self.root)

        for target in ("info", "methods", "both"):
            with self.subTest(target=target):
                _make_release(self.root)
                _add_crystalnets_parent_columns(self.root)
                _declare_crystalnets_reference_contracts(self.root)
                info = json.loads(info_path.read_text(encoding="utf-8"))
                methods = json.loads(methods_path.read_text(encoding="utf-8"))
                if target in ("info", "both"):
                    info["release_status"] = "FINAL"
                if target in ("methods", "both"):
                    methods["release_status"] = "FINAL"
                info_path.write_text(
                    json.dumps(info, sort_keys=True), encoding="utf-8"
                )
                methods_path.write_text(
                    json.dumps(methods, sort_keys=True), encoding="utf-8"
                )
                with self.assertRaisesRegex(
                    ReleaseValidationError,
                    "retired or reserved declaration key",
                ):
                    CoREMOFDataset.from_release(self.root)

    def test_optional_structure_matcher_parent_columns_are_validated_and_loaded(self):
        _add_structure_matcher_parent_columns(self.root)
        _declare_structure_matcher_contract(self.root)

        dataset = CoREMOFDataset.from_release(self.root)
        first = dataset[dataset.structure_ids[0]]
        self.assertEqual(
            first.parent_group("structure_matcher_strict").group_id,
            "SM-ABCD0001",
        )
        self.assertEqual(first.parent_group("sm").group_id, "SM-ABCD0001")
        self.assertIn(
            "parent_groups/structure_matcher_strict_evidence_receipt.json",
            dataset.input_hashes,
        )

    def test_structure_matcher_columns_require_the_exact_method_contract(self):
        _add_structure_matcher_parent_columns(self.root)
        methods_path = self.root / "parent_groups" / "parent_group_methods.json"
        methods = json.loads(methods_path.read_text(encoding="utf-8"))
        methods["csv_column_prefixes"]["structure_matcher_strict"] = "sm"
        methods_path.write_text(json.dumps(methods), encoding="utf-8")
        with self.assertRaisesRegex(
            ReleaseValidationError, "criteria.structure_matcher_strict"
        ):
            CoREMOFDataset.from_release(self.root)

    def test_structure_matcher_scientific_contract_fields_are_exact(self):
        _add_structure_matcher_parent_columns(self.root)
        cases = (
            ("method_id", "different_method"),
            ("method_schema_version", "different-schema/1.0"),
            ("authoritative_evidence", "ONE_WAY_EDGES"),
            ("fit_policy", "ONE_WAY_ALLOWED"),
            ("component_semantics", "ALL_PAIRS_MATCH"),
            ("component_completeness_policy", "ALLOW_PARTIAL"),
            ("public_status_policy", ["MATCHED", "PARTIAL"]),
            ("included_in_priority_main", True),
            ("included_in_priority_main", 0),
            ("included_in_main_union", True),
            ("included_in_main_union", 0),
            ("historical_relaxed_executed", True),
            ("historical_relaxed_executed", 0),
            ("historical_relaxed_exposed", 0),
        )
        for key, invalid in cases:
            with self.subTest(key=key):
                methods_path, _ = _declare_structure_matcher_contract(self.root)
                methods = json.loads(methods_path.read_text(encoding="utf-8"))
                methods["criteria"]["structure_matcher_strict"][key] = invalid
                methods_path.write_text(json.dumps(methods), encoding="utf-8")
                with self.assertRaisesRegex(ReleaseValidationError, key):
                    CoREMOFDataset.from_release(self.root)

    def test_structure_matcher_receipt_hash_and_relaxed_policy_fail_closed(self):
        _add_structure_matcher_parent_columns(self.root)
        methods_path, receipt_path = _declare_structure_matcher_contract(self.root)
        receipt_path.write_text("{}", encoding="utf-8")
        with self.assertRaisesRegex(ReleaseValidationError, "receipt SHA-256"):
            CoREMOFDataset.from_release(self.root)

        _declare_structure_matcher_contract(self.root)
        methods = json.loads(methods_path.read_text(encoding="utf-8"))
        methods["criteria"]["structure_matcher_strict"][
            "historical_relaxed_exposed"
        ] = True
        methods_path.write_text(json.dumps(methods), encoding="utf-8")
        with self.assertRaisesRegex(
            ReleaseValidationError, "historical_relaxed_exposed"
        ):
            CoREMOFDataset.from_release(self.root)

    def test_structure_matcher_receipt_schema_and_counts_fail_closed(self):
        _add_structure_matcher_parent_columns(self.root)
        methods_path, receipt_path = _declare_structure_matcher_contract(self.root)

        def rewrite_receipt(**changes):
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            receipt.update(changes)
            receipt_path.write_text(json.dumps(receipt, sort_keys=True), encoding="utf-8")
            methods = json.loads(methods_path.read_text(encoding="utf-8"))
            methods["criteria"]["structure_matcher_strict"]["evidence_receipt"][
                "sha256"
            ] = hashlib.sha256(receipt_path.read_bytes()).hexdigest()
            methods_path.write_text(json.dumps(methods, sort_keys=True), encoding="utf-8")

        rewrite_receipt(strict_pair_ledger_sha256="not-a-sha256")
        with self.assertRaisesRegex(ReleaseValidationError, "strict_pair_ledger_sha256"):
            CoREMOFDataset.from_release(self.root)

        _declare_structure_matcher_contract(self.root)
        rewrite_receipt(unresolved_pair_count=1)
        with self.assertRaisesRegex(ReleaseValidationError, "candidate_pair_count"):
            CoREMOFDataset.from_release(self.root)

        for field in (
            "historical_relaxed_executed",
            "historical_relaxed_exposed",
        ):
            with self.subTest(field=field):
                _declare_structure_matcher_contract(self.root)
                rewrite_receipt(**{field: 0})
                with self.assertRaisesRegex(ReleaseValidationError, field):
                    CoREMOFDataset.from_release(self.root)

    def test_structure_matcher_receipt_binds_the_parent_table(self):
        rows = _add_structure_matcher_parent_columns(self.root)
        _declare_structure_matcher_contract(self.root)
        parents_path = self.root / "parent_groups" / "parent_groups.csv"
        for row in rows:
            row["unbound_note"] = "changed-after-receipt"
        _write_csv(parents_path, tuple(rows[0]), rows)
        with self.assertRaisesRegex(ReleaseValidationError, "parent_groups.csv"):
            CoREMOFDataset.from_release(self.root)

    def test_manifest_is_optional(self):
        (self.root / "manifests" / "cif_manifest.csv").unlink()
        dataset = CoREMOFDataset.from_release(self.root)
        self.assertEqual(dict(dataset.cif_hashes), {})
        self.assertNotIn("manifests/cif_manifest.csv", dataset.input_hashes)
        with self.assertRaisesRegex(
            ReleaseValidationError, "requires manifests/cif_manifest.csv"
        ):
            CoREMOFDataset.from_release(self.root, verify_cif_files=True)

    def test_optional_cif_byte_verification_is_strict(self):
        content = b"0123456789"
        cif_directory = self.root / "cifs"
        cif_directory.mkdir()
        for row in _metadata_rows():
            (cif_directory / (row["structure_id"] + ".cif")).write_bytes(content)

        manifest_path = self.root / "manifests" / "cif_manifest.csv"
        with manifest_path.open(encoding="utf-8", newline="") as handle:
            manifest_rows = list(csv.DictReader(handle))
        expected_hash = hashlib.sha256(content).hexdigest()
        for row in manifest_rows:
            row["sha256"] = expected_hash
        _write_csv(manifest_path, tuple(manifest_rows[0]), manifest_rows)

        verified = CoREMOFDataset.from_release(
            self.root, verify_cif_files=True
        )
        self.assertEqual(len(verified), 4)
        from CoREMOF.splitters import split_release

        strict_split = split_release(
            verified.classify("5checker"),
            parent_method="none",
            verify_cif_files=True,
        )
        self.assertTrue(strict_split.cif_files_verified)
        self.assertTrue(strict_split.receipt()["integrity"]["cif_files_verified"])
        first = cif_directory / (_metadata_rows()[0]["structure_id"] + ".cif")
        first.write_bytes(b"changed!!!")
        with self.assertRaisesRegex(ReleaseValidationError, "SHA-256"):
            CoREMOFDataset.from_release(self.root, verify_cif_files=True)

    def test_cif_verification_and_official_view_flags_require_exact_booleans(self):
        for value in (0, 1, None, "false", [], {}):
            with self.subTest(api="from_release", value=value):
                with self.assertRaisesRegex(TypeError, "verify_cif_files"):
                    CoREMOFDataset.from_release(
                        self.root, verify_cif_files=value
                    )

        dataset = CoREMOFDataset.from_release(self.root)
        for value in (0, 1, None, "false", [], {}):
            with self.subTest(api="dataset_constructor", value=value):
                with self.assertRaisesRegex(TypeError, "cif_files_verified"):
                    CoREMOFDataset(
                        dataset.release_root,
                        dataset.records,
                        dataset.dataset_info,
                        dataset.parent_group_methods,
                        dataset.parent_by_id,
                        dataset.input_hashes,
                        value,
                    )
            with self.subTest(api="classified_constructor", value=value):
                with self.assertRaisesRegex(TypeError, "checker_view_official"):
                    dataset_module.ClassifiedDataset(
                        dataset, "5checker", (), checker_view_official=value
                    )

        custom = dataset_module.ClassifiedDataset(dataset, "custom:test", ())
        self.assertFalse(custom.checker_view_official)
        with self.assertRaisesRegex(
            ValueError, "requires a canonical official checker preset"
        ):
            dataset_module.ClassifiedDataset(
                dataset, "custom:test", (), checker_view_official=True
            )
        with self.assertRaisesRegex(
            ValueError, "internal validated-release checker-recomputation path"
        ):
            dataset_module.ClassifiedDataset(
                dataset, "5checker", (), checker_view_official=True
            )
        official = dataset.classify("5checker")
        self.assertTrue(official.checker_view_official)

        forged_label = replace(official.records[0], label="NCR")
        with self.assertRaisesRegex(ValueError, "internal validated-release"):
            dataset_module.ClassifiedDataset(
                dataset,
                "5checker",
                (forged_label,),
                checker_view_official=True,
                _official_checker_view_token=dataset_module._OFFICIAL_CHECKER_VIEW_TOKEN,
            )
        missing_status = replace(official.records[0], checker_statuses={})
        with self.assertRaisesRegex(ValueError, "internal validated-release"):
            dataset_module.ClassifiedDataset(
                dataset,
                "5checker",
                (missing_status,),
                checker_view_official=True,
                _official_checker_view_token=dataset_module._OFFICIAL_CHECKER_VIEW_TOKEN,
            )

        generic = CoREMOFDataset(
            dataset.release_root,
            dataset.records,
            dataset.dataset_info,
            dataset.parent_group_methods,
            dataset.parent_by_id,
            dataset.input_hashes,
            False,
        )
        self.assertFalse(generic.classify("5checker").checker_view_official)

        for value in (True, 1, 1.5, [], {}, None, "", " vtest", "vtest "):
            with self.subTest(dataset_version=value):
                malformed_info = dict(dataset.dataset_info)
                malformed_info["dataset_version"] = value
                malformed = CoREMOFDataset(
                    dataset.release_root,
                    dataset.records,
                    malformed_info,
                    dataset.parent_group_methods,
                    dataset.parent_by_id,
                    dataset.input_hashes,
                    False,
                )
                with self.assertRaisesRegex(
                    (TypeError, ValueError), "dataset_version"
                ):
                    _ = malformed.dataset_version

        for value in (" vtest", "vtest "):
            with self.subTest(loader_dataset_version=value):
                _make_release(self.root)
                info_path = self.root / "dataset_info.json"
                methods_path = (
                    self.root / "parent_groups" / "parent_group_methods.json"
                )
                info = json.loads(info_path.read_text(encoding="utf-8"))
                methods = json.loads(methods_path.read_text(encoding="utf-8"))
                info["dataset_version"] = value
                methods["dataset_version"] = value
                info_path.write_text(json.dumps(info), encoding="utf-8")
                methods_path.write_text(json.dumps(methods), encoding="utf-8")
                with self.assertRaisesRegex(
                    ReleaseValidationError,
                    "dataset_version must be an exact nonblank string",
                ):
                    CoREMOFDataset.from_release(self.root)

    def test_public_id_and_identity_fields_must_agree(self):
        path = self.root / "metadata" / "metadata.csv"
        rows = _metadata_rows()
        rows[0]["structure_variant"] = "FSR"
        _write_csv(path, tuple(rows[0]), rows)
        with self.assertRaisesRegex(
            ReleaseValidationError, "disagrees with its public ID"
        ):
            CoREMOFDataset.from_release(self.root)

    def test_import_does_not_pull_heavy_scientific_packages(self):
        code = (
            "import sys; import CoREMOF.dataset; "
            "blocked={'numpy','pandas','pymatgen','sklearn'}; "
            "loaded=blocked.intersection(sys.modules); "
            "assert not loaded, sorted(loaded)"
        )
        completed = subprocess.run(
            [sys.executable, "-c", code],
            cwd=str(Path(__file__).resolve().parents[1]),
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_filter_helpers_compose(self):
        classified = CoREMOFDataset.from_release(self.root).classify("5checker")
        self.assertFalse(classified.selection_filters["applied"])
        selected = classified.filter(
            labels=("NCR",), sources="cod", variants="fsr", metals="fe"
        )
        self.assertEqual(selected.structure_ids, ("FSR-COD-2026-0001",))
        self.assertTrue(selected.selection_filters["applied"])
        self.assertEqual(selected.selection_filters["selected_count"], 1)
        self.assertEqual(len(selected.selection_filters["steps"]), 1)
        sequential = classified.filter(sources=("COD", "SI")).filter(
            labels=("CR", "AMBIGUOUS")
        )
        self.assertEqual(
            sequential.structure_ids,
            ("ASR-COD-2026-0001", "ASR-SI-2025-0001"),
        )
        self.assertEqual(len(sequential.selection_filters["steps"]), 2)
        self.assertEqual(classified.select(metals="Cu").labels, ("CR",))
        self.assertEqual(classified.filter(labels="cr").structure_ids, classified.cr_ids)
        self.assertEqual(
            classified.filter(structure_ids=("ASR-SI-2025-0001",)).structure_ids,
            ("ASR-SI-2025-0001",),
        )
        with self.assertRaises(ValueError):
            classified.filter(labels="maybe")
        with self.assertRaisesRegex(KeyError, "unknown structure_id"):
            classified.filter(structure_ids=("absent",))
        for value in (True, 1, 1.5, [], {}, " ASR-COD-2026-0001"):
            with self.subTest(structure_id=value):
                with self.assertRaises((TypeError, ValueError)):
                    classified.filter(structure_ids=(value,))
        for value in (
            {"COD", "SI"},
            frozenset(("COD", "SI")),
            {"COD": True},
            (item for item in ("COD", "SI")),
            b"COD",
        ):
            with self.subTest(unordered_or_one_shot=type(value).__name__):
                with self.assertRaisesRegex(TypeError, "ordered list/tuple"):
                    classified.filter(sources=value)

    def test_custom_view_does_not_read_a_published_label_column(self):
        dataset = CoREMOFDataset.from_release(self.root)
        custom = dataset.classify(["MOFClassifier", "MOFChecker"])
        self.assertEqual(
            custom.checker_preset, "custom:MOFClassifier+MOFChecker"
        )
        # The final fixture row is UNCHECKED in every published preset because
        # Chen-Manz is unavailable, but it is complete in this two-checker view.
        self.assertEqual(custom["ION-CSD-2024-0001"].label, "CR")
        self.assertEqual(
            custom.label_by_id["ION-CSD-2024-0001"], "CR"
        )
        self.assertEqual(custom.checker_view, custom.checker_preset)
        self.assertFalse(custom.checker_view_official)

    def test_lazy_splitter_receives_filtered_dataset_and_safe_defaults(self):
        calls = {}
        fake_module = ModuleType("CoREMOF.splitters")

        class FakeSplitter:
            def __init__(self, dataset, **kwargs):
                calls["dataset"] = dataset
                calls["kwargs"] = kwargs

            def train_valid_test_split(self, fractions):
                calls["fractions"] = fractions
                return "split-result"

        fake_module.ParentGroupSplitter = FakeSplitter
        classified = CoREMOFDataset.from_release(self.root).classify(5)
        with patch.dict(sys.modules, {"CoREMOF.splitters": fake_module}):
            result = classified.train_valid_test_split()
        self.assertEqual(result, "split-result")
        self.assertIs(calls["dataset"], classified)
        self.assertEqual(calls["kwargs"]["labels"], ("CR", "NCR"))
        self.assertEqual(calls["kwargs"]["missing_parent"], "singleton")
        self.assertEqual(calls["kwargs"]["parent_method"], "priority_main")
        self.assertEqual(calls["fractions"], (0.8, 0.1, 0.1))

    def test_mismatched_published_label_is_rejected(self):
        path = self.root / "metadata" / "metadata.csv"
        rows = _metadata_rows()
        rows[0]["label_3checker"] = "NCR"
        _write_csv(path, tuple(rows[0]), rows)
        with self.assertRaisesRegex(ReleaseValidationError, "recompute to 'CR'"):
            CoREMOFDataset.from_release(self.root)

    def test_release_rejects_internal_checker_execution_tokens(self):
        path = self.root / "metadata" / "metadata.csv"
        rows = _metadata_rows()
        rows[0]["mofclassifier_status"] = "TIMEOUT"
        rows[0]["label_3checker"] = "UNCHECKED"
        rows[0]["label_4checker"] = "UNCHECKED"
        rows[0]["label_5checker"] = "UNCHECKED"
        _write_csv(path, tuple(rows[0]), rows)
        with self.assertRaisesRegex(ReleaseValidationError, "non-public MOFClassifier"):
            CoREMOFDataset.from_release(self.root)

    def test_parent_group_size_is_recomputed(self):
        path = self.root / "parent_groups" / "parent_groups.csv"
        rows = _parent_rows()
        rows[0]["rac_size"] = "3"
        _write_csv(path, tuple(rows[0]), rows)
        with self.assertRaisesRegex(ReleaseValidationError, "declares size 3"):
            CoREMOFDataset.from_release(self.root)

    def test_parent_group_size_rejects_signed_and_non_ascii_digits(self):
        path = self.root / "parent_groups" / "parent_groups.csv"
        for value in ("+1", "-1", "١", " 1 ", "01"):
            with self.subTest(value=value):
                rows = _parent_rows()
                rows[2]["rac_size"] = value
                _write_csv(path, tuple(rows[0]), rows)
                with self.assertRaisesRegex(ReleaseValidationError, "invalid rac size"):
                    CoREMOFDataset.from_release(self.root)

    def test_parent_group_id_prefix_and_hash_are_validated(self):
        path = self.root / "parent_groups" / "parent_groups.csv"
        rows = _parent_rows()
        rows[0]["rac_group"] = "R-not-a-public-code"
        _write_csv(path, tuple(rows[0]), rows)
        with self.assertRaisesRegex(ReleaseValidationError, "invalid rac group ID"):
            CoREMOFDataset.from_release(self.root)

    def test_parent_method_declaration_is_exact(self):
        path = self.root / "parent_groups" / "parent_group_methods.json"
        methods = json.loads(path.read_text(encoding="utf-8"))
        methods["csv_column_prefixes"]["not_a_method"] = methods[
            "csv_column_prefixes"
        ].pop("rac5")
        path.write_text(json.dumps(methods), encoding="utf-8")
        with self.assertRaisesRegex(ReleaseValidationError, "public method contract"):
            CoREMOFDataset.from_release(self.root)

    def test_parent_method_prefix_values_must_be_strings(self):
        path = self.root / "parent_groups" / "parent_group_methods.json"
        for value in ([], {}, True, 1):
            with self.subTest(value=value):
                methods = json.loads(path.read_text(encoding="utf-8"))
                methods["csv_column_prefixes"]["rac5"] = value
                path.write_text(json.dumps(methods), encoding="utf-8")
                with self.assertRaisesRegex(ReleaseValidationError, "must all be strings"):
                    CoREMOFDataset.from_release(self.root)
                _make_release(self.root)

    def test_orphan_and_unknown_parent_columns_fail_closed(self):
        path = self.root / "parent_groups" / "parent_groups.csv"
        with path.open("r", encoding="utf-8", newline="") as handle:
            original_rows = list(csv.DictReader(handle))

        orphan_rows = [dict(row) for row in original_rows]
        for row in orphan_rows:
            row["unregistered_group"] = "U-ABCD0001"
        _write_csv(path, tuple(orphan_rows[0]), orphan_rows)
        with self.assertRaisesRegex(
            ReleaseValidationError, "complete status/group/size"
        ):
            CoREMOFDataset.from_release(self.root)

        unknown_rows = [dict(row) for row in original_rows]
        for index, row in enumerate(unknown_rows):
            row["unregistered_status"] = "UNMATCHED"
            row["unregistered_group"] = "U-{:08X}".format(index)
            row["unregistered_size"] = "1"
        _write_csv(path, tuple(unknown_rows[0]), unknown_rows)
        with self.assertRaisesRegex(ReleaseValidationError, "CSV prefixes"):
            CoREMOFDataset.from_release(self.root)

        extra_rows = [dict(row) for row in original_rows]
        for row in extra_rows:
            row["audit_note"] = "not-part-of-the-public-schema"
        _write_csv(path, tuple(extra_rows[0]), extra_rows)
        with self.assertRaisesRegex(ReleaseValidationError, "must contain exactly"):
            CoREMOFDataset.from_release(self.root)

    def test_relaxed_parent_field_is_always_rejected(self):
        path = self.root / "parent_groups" / "parent_groups.csv"
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        for row in rows:
            row["historical_relaxed_note"] = "not-public"
        _write_csv(path, tuple(rows[0]), rows)
        with self.assertRaisesRegex(ReleaseValidationError, "historical relaxed"):
            CoREMOFDataset.from_release(self.root)

    def test_duplicate_and_cross_table_identifier_errors_are_rejected(self):
        metadata_path = self.root / "metadata" / "metadata.csv"
        rows = _metadata_rows()
        rows.append(dict(rows[0]))
        _write_csv(metadata_path, tuple(rows[0]), rows)
        with self.assertRaisesRegex(ReleaseValidationError, "duplicate structure_id"):
            CoREMOFDataset.from_release(self.root)

        _make_release(self.root)
        manifest_path = self.root / "manifests" / "cif_manifest.csv"
        with manifest_path.open(encoding="utf-8", newline="") as handle:
            manifest_rows = list(csv.DictReader(handle))
        _write_csv(manifest_path, tuple(manifest_rows[0]), manifest_rows[:-1])
        with self.assertRaisesRegex(ReleaseValidationError, "structure_id sets differ"):
            CoREMOFDataset.from_release(self.root)

    def test_duplicate_csv_header_is_rejected(self):
        path = self.root / "metadata" / "metadata.csv"
        text = path.read_text(encoding="utf-8")
        first, rest = text.split("\n", 1)
        path.write_text(first + ",structure_id\n" + rest, encoding="utf-8")
        with self.assertRaisesRegex(ReleaseValidationError, "duplicate header"):
            CoREMOFDataset.from_release(self.root)

    def test_release_snapshot_rejects_same_byte_path_replacement(self):
        metadata_path = self.root / "metadata" / "metadata.csv"
        original_reader = dataset_module._read_release_snapshot_bytes
        replaced = False

        def replace_during_capture(handle):
            nonlocal replaced
            data = original_reader(handle)
            if Path(handle.name) == metadata_path and not replaced:
                replacement = metadata_path.with_name("metadata.replacement")
                replacement.write_bytes(data)
                replacement.replace(metadata_path)
                replaced = True
            return data

        with patch.object(
            dataset_module,
            "_read_release_snapshot_bytes",
            side_effect=replace_during_capture,
        ):
            with self.assertRaisesRegex(
                ReleaseValidationError, "changed or was replaced"
            ):
                CoREMOFDataset.from_release(self.root)
        self.assertTrue(replaced)


if __name__ == "__main__":
    unittest.main()
