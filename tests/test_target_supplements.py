"""Synthetic regression tests for target-free loading and optional attachment."""

import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import CoREMOF.splitters as splitters_module
import CoREMOF.target_supplements as supplements_module
import CoREMOF.targets as targets_module
from CoREMOF.dataset import CoREMOFDataset, ReleaseValidationError
from CoREMOF.target_supplements import SUPPLEMENT_KIND
from CoREMOF.targets import TargetDataError, TargetSource, merge_targets
from CoREMOF.projections import export_source_projection
from test_dataset_labels import _make_release, _metadata_rows, _write_csv


TARGETS = (
    "ch4_loading_298K_65bar_mol_kg_framework",
    "co2_n2_henry_selectivity_298K",
    "h2_loading_77K_100bar_mol_kg_framework",
)
CORE_FILES = (
    "dataset_info.json", "metadata/metadata.csv", "manifests/cif_manifest.csv"
)


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_csv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)


class TargetSupplementTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.core = self.root / "core"
        _make_release(self.core)
        info_path = self.core / "dataset_info.json"
        info = json.loads(info_path.read_text())
        info["metadata_layout"] = "target_free_core/1.0"
        info_path.write_text(json.dumps(info))
        self.dataset = CoREMOFDataset.from_release(self.core)
        self.supplement = self.root / "targets"
        self.supplement.mkdir()
        self.definitions = {
            TARGETS[0]: {
                "unit": "mol/kg-framework",
                "conditions": {"temperature_K": 298.0, "pressure_Pa": 6500000.0},
            },
            TARGETS[1]: {
                "unit": "dimensionless",
                "conditions": {"temperature_K": 298.0, "pressure_Pa": None},
            },
            TARGETS[2]: {
                "unit": "mol/kg-framework",
                "conditions": {"temperature_K": 77.0, "pressure_Pa": 10000000.0},
            },
        }
        self.payloads = []
        csv_rows = []
        for index, identity in enumerate(self.dataset.structure_ids):
            value = (0.0, None, 1.25, 5.0)[index]
            status = "SUCCESS" if value is not None else "NOT_AVAILABLE"
            row = {"structure_id": identity}
            payloads = {}
            for name in TARGETS:
                row[name] = "" if value is None else str(value)
                row[name + "__status"] = status
                payloads[name] = dict(
                    self.definitions[name], value=value,
                    available=value is not None, execution_status=status,
                    diagnostic=None, calculation_method={"raspa_version": "RASPA3"},
                )
            csv_rows.append(row)
            self.payloads.append({"structure_id": identity, "targets": payloads})
        _write_csv(self.supplement / "targets.csv", tuple(csv_rows[0]), csv_rows)
        self._write_payloads()
        (self.supplement / "targets_schema.json").write_text(json.dumps({"schema_version": "targets/1.0"}))
        (self.supplement / "target_details.json").write_text(json.dumps({"simulation_software": "RASPA3"}))
        self.manifest = {
            "kind": SUPPLEMENT_KIND,
            "schema_version": SUPPLEMENT_KIND,
            "dataset_version": self.dataset.dataset_version,
            "structure_count": len(self.dataset),
            "compatible_core_sha256": {
                name: self.dataset.input_hashes[name] for name in CORE_FILES
            },
            "target_definitions": self.definitions,
            "files": {},
        }
        self._refresh()

    def _write_payloads(self):
        (self.supplement / "targets.jsonl").write_text(
            "".join(json.dumps(row, allow_nan=False) + "\n" for row in self.payloads)
        )

    def _refresh(self):
        for name in ("targets.csv", "targets.jsonl", "targets_schema.json", "target_details.json"):
            path = self.supplement / name
            self.manifest["files"][name] = {
                "sha256": _sha(path), "bytes": path.stat().st_size
            }
        path = self.supplement / "manifest.json"
        path.write_text(json.dumps(self.manifest, allow_nan=False))
        self.expected = _sha(path)

    def _attach(self):
        return self.dataset.attach_target_supplement(
            self.supplement, expected_sha256=self.expected
        )

    def _reject(self):
        with self.assertRaises(TargetDataError):
            self._attach()
        self.assertTrue(all(not set(TARGETS).intersection(row) for row in self.dataset.metadata_rows))

    def test_explicit_attachment_preserves_zero_null_ids_groups_and_hashes(self):
        before = {name: _sha(self.supplement / name) for name in self.manifest["files"]}
        attached = self._attach()
        self.assertEqual(attached.structure_ids, self.dataset.structure_ids)
        self.assertEqual(attached.parent_by_id, self.dataset.parent_by_id)
        self.assertEqual(attached.target_columns, TARGETS)
        first, second = self.dataset.structure_ids[:2]
        self.assertEqual(attached.target_values_by_id[first][TARGETS[0]], 0.0)
        self.assertIsNone(attached.target_values_by_id[second][TARGETS[0]])
        self.assertEqual(attached.target_definitions[TARGETS[1]]["unit"], "dimensionless")
        for row in attached.metadata_rows:
            self.assertTrue(all(name + "__status" not in row for name in TARGETS))
        receipt = attached.receipt()
        self.assertEqual(receipt["sources"][0]["sha256"], self.manifest["files"]["targets.csv"]["sha256"])
        self.assertEqual(receipt["sources"][0]["name"], "target-supplement:" + self.expected)
        self.assertEqual(before, {name: _sha(self.supplement / name) for name in before})
        for name, digest in self.dataset.input_hashes.items():
            self.assertEqual(attached.input_hashes[name], digest)

    def test_default_loading_never_discovers_corrupt_or_absent_supplements(self):
        for state in ("missing", "corrupt"):
            with self.subTest(state=state):
                if state == "corrupt":
                    (self.core / "target_supplement").mkdir()
                    (self.core / "target_supplement" / "manifest.json").write_text("not JSON")
                loaded = CoREMOFDataset.from_release(self.core)
                self.assertEqual(loaded.structure_ids, self.dataset.structure_ids)
                self.assertTrue(all(not set(TARGETS).intersection(row) for row in loaded.metadata_rows))

    def test_supplement_receipt_binds_imported_validator_and_writes_it(self):
        attached = self._attach()
        expected = dict(targets_module._IMPORTED_IMPLEMENTATION_HASHES)
        expected["target_supplements.py"] = _sha(Path(supplements_module.__file__))
        self.assertEqual(
            expected["target_supplements.py"],
            supplements_module._IMPORTED_IMPLEMENTATION_HASHES["target_supplements.py"],
        )
        self.assertEqual(attached.receipt()["implementation"]["source_sha256"], expected)
        receipt_path = attached.write(self.root / "attached")[2]
        self.assertEqual(
            json.loads(receipt_path.read_text())["implementation"]["source_sha256"],
            expected,
        )

    def test_generic_merge_and_split_keep_legacy_implementation_closure(self):
        changed = {"target_supplements.py": "0" * 64}
        with patch.object(supplements_module, "_current_implementation_hashes", return_value=changed):
            generic = self.dataset.merge_targets((TargetSource(
                self.supplement / "targets.csv", target_columns=TARGETS,
                value_types={name: "float" for name in TARGETS},
            ),))
            receipt = generic.receipt()
            self.assertEqual(
                receipt["implementation"]["source_sha256"],
                dict(targets_module._IMPORTED_IMPLEMENTATION_HASHES),
            )
            split = generic.classify("5checker").train_valid_test_split(
                parent_method="none", leakage_guard="parent_only",
                fractions=(1.0, 0.0, 0.0), random_state=7,
            )
            self.assertNotIn(
                "target_supplements.py", split.receipt()["implementation"]["source_sha256"]
            )

    def test_validator_drift_before_attachment_is_rejected(self):
        changed = {"target_supplements.py": "0" * 64}
        with patch.object(supplements_module, "_current_implementation_hashes", return_value=changed):
            with self.assertRaisesRegex(TargetDataError, "source changed after module import: target_supplements.py"):
                self._attach()

    def test_validator_drift_during_validation_is_rejected(self):
        original_parse = supplements_module._parse_csv
        changed = {"target_supplements.py": "0" * 64}

        def drift_after_parse(*args):
            rows = original_parse(*args)
            drift = patch.object(supplements_module, "_current_implementation_hashes", return_value=changed)
            drift.start()
            self.addCleanup(drift.stop)
            return rows

        with patch.object(supplements_module, "_parse_csv", side_effect=drift_after_parse):
            with self.assertRaisesRegex(TargetDataError, "source changed after module import: target_supplements.py"):
                self._attach()

    def test_validator_drift_after_attachment_blocks_receipt_and_split(self):
        attached = self._attach()
        changed = {"target_supplements.py": "0" * 64}
        with patch.object(supplements_module, "_current_implementation_hashes", return_value=changed):
            with self.assertRaisesRegex(TargetDataError, "source changed after module import: target_supplements.py"):
                attached.receipt()
            with self.assertRaisesRegex(TargetDataError, "source changed after module import: target_supplements.py"):
                attached.classify("5checker").train_valid_test_split(
                    parent_method="none", leakage_guard="parent_only",
                    fractions=(1.0, 0.0, 0.0), random_state=7,
                )

    def test_supplement_split_propagates_and_checks_validator_identity(self):
        attached = self._attach()
        split = attached.classify("5checker").train_valid_test_split(
            parent_method="none", leakage_guard="parent_only",
            fractions=(1.0, 0.0, 0.0), random_state=7,
        )
        receipt = split.receipt()
        digest = _sha(Path(supplements_module.__file__))
        self.assertEqual(receipt["implementation"]["source_sha256"]["target_supplements.py"], digest)
        self.assertEqual(receipt["target_data"]["implementation"]["source_sha256"]["target_supplements.py"], digest)
        self.assertFalse(receipt["official_split"])
        tampered = attached.receipt()
        tampered["implementation"]["source_sha256"]["target_supplements.py"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "changed between target merge and split import"):
            splitters_module._implementation_hashes(include_targets=True, target_receipt=tampered)
        changed = {"target_supplements.py": "0" * 64}
        with patch.object(supplements_module, "_current_implementation_hashes", return_value=changed):
            with self.assertRaisesRegex(ValueError, "source changed after module import: target_supplements.py"):
                split.receipt()

    def test_reimported_validator_cannot_rebind_an_existing_receipt(self):
        attached = self._attach()
        changed = {"target_supplements.py": "0" * 64}
        with patch.object(supplements_module, "_IMPORTED_IMPLEMENTATION_HASHES", changed):
            with patch.object(supplements_module, "_current_implementation_hashes", return_value=changed):
                with self.assertRaisesRegex(TargetDataError, "does not match executing sources"):
                    attached.receipt()

    def test_supplement_identity_requires_authenticated_factory(self):
        source = TargetSource(self.supplement / "targets.csv", target_columns=TARGETS)
        with self.assertRaisesRegex(TargetDataError, "only by attach_target_supplement"):
            merge_targets(self.dataset, (source,), _supplement_factory_token=object())

    def test_legacy_columns_hidden_by_default_and_opt_in_preserves_original_hash(self):
        path = self.core / "metadata" / "metadata.csv"
        fields, rows = _read_csv(path)
        for row in rows:
            for name in TARGETS:
                row[name] = "0"
                row[name + "__status"] = "SUCCESS"
        _write_csv(path, tuple(rows[0]), rows)
        info_path = self.core / "dataset_info.json"
        info = json.loads(info_path.read_text())
        info.pop("metadata_layout")
        info["target_metadata"] = {"coverage": 4}
        info["metadata_revision"] = {
            "revision_id": "target_metadata_20260923_v3",
            "targets_modified": True,
            "frozen_splits_modified": False,
            "publication_gate": "unchanged",
        }
        info["tabular_files"] = {
            "metadata/targets.jsonl": {"row_count": 4, "size_bytes": 2000},
            "metadata/targets_schema.json": {"size_bytes": 1000},
            "metadata/metadata.csv": {
                "columns": list(rows[0]), "row_count": 4, "size_bytes": path.stat().st_size,
            },
            "features/zeo_features.csv": {
                "columns": ["structure_id", "density", "target_columns_consumed"],
                "row_count": 4, "size_bytes": 900,
            },
        }
        info_path.write_text(json.dumps(info))
        loaded = CoREMOFDataset.from_release(self.core)
        legacy = CoREMOFDataset.from_release(self.core, include_legacy_targets=True)
        self.assertTrue(all(not set(TARGETS).intersection(row) for row in loaded.metadata_rows))
        self.assertNotIn("target_metadata", loaded.dataset_info)
        self.assertNotIn("targets_modified", loaded.dataset_info["metadata_revision"])
        self.assertEqual(loaded.dataset_info["metadata_revision"]["revision_id"], "legacy_core_view")
        self.assertEqual(loaded.dataset_info["metadata_revision"]["publication_gate"], "unchanged")
        self.assertEqual(legacy.metadata_rows[0][TARGETS[0]], "0")
        self.assertIn("target_metadata", legacy.dataset_info)
        self.assertNotIn("metadata/targets.jsonl", loaded.dataset_info["tabular_files"])
        self.assertNotIn("metadata/targets_schema.json", loaded.dataset_info["tabular_files"])
        self.assertFalse(set(TARGETS).intersection(loaded.dataset_info["tabular_files"]["metadata/metadata.csv"]["columns"]))
        self.assertEqual(loaded.dataset_info["tabular_files"]["metadata/metadata.csv"]["size_bytes"], path.stat().st_size)
        self.assertEqual(loaded.dataset_info["tabular_files"]["features/zeo_features.csv"], legacy.dataset_info["tabular_files"]["features/zeo_features.csv"])
        self.assertIn("metadata/targets.jsonl", legacy.dataset_info["tabular_files"])
        self.assertEqual(loaded.input_hashes, legacy.input_hashes)
        self.assertEqual(loaded.input_hashes["metadata/metadata.csv"], _sha(path))

    def test_legacy_opt_in_requires_boolean(self):
        for value in (None, 0, "yes"):
            with self.subTest(value=value), self.assertRaises(TypeError):
                CoREMOFDataset.from_release(self.core, include_legacy_targets=value)

    def test_strict_core_rejects_target_columns_even_with_legacy_opt_in(self):
        path = self.core / "metadata" / "metadata.csv"
        fields, rows = _read_csv(path)
        for row in rows:
            row[TARGETS[0]] = "0"
        _write_csv(path, tuple(rows[0]), rows)
        for flag in (False, True):
            with self.subTest(flag=flag), self.assertRaises(ReleaseValidationError):
                CoREMOFDataset.from_release(self.core, include_legacy_targets=flag)

    def test_strict_core_rejects_target_metadata_without_value_columns(self):
        path = self.core / "dataset_info.json"
        info = json.loads(path.read_text())
        info["target_metadata"] = {}
        path.write_text(json.dumps(info))
        for flag in (False, True):
            with self.subTest(flag=flag), self.assertRaises(ReleaseValidationError):
                CoREMOFDataset.from_release(self.core, include_legacy_targets=flag)

    def test_strict_core_rejects_target_table_column_and_revision_declarations(self):
        path = self.core / "dataset_info.json"
        original = json.loads(path.read_text())
        for addition in (
            {"tabular_files": {"metadata/targets.jsonl": {"row_count": 4}}},
            {"tabular_files": {"metadata/targets_schema.json": {"size_bytes": 100}}},
            {"tabular_files": {"metadata/metadata.csv": {"columns": ["structure_id", TARGETS[0]]}}},
            {"metadata_revision": {"targets_modified": False}},
            {"metadata_revision": {"revision_id": "target_metadata_legacy"}},
        ):
            with self.subTest(addition=addition):
                info = dict(original, **addition)
                path.write_text(json.dumps(info))
                for flag in (False, True):
                    with self.assertRaises(ReleaseValidationError):
                        CoREMOFDataset.from_release(self.core, include_legacy_targets=flag)

    def test_target_free_analysis_rule_does_not_consume_targets(self):
        # An assertion of independence is not a stored response or an implicit
        # join. This supplemental analysis-policy file is not discovered at load.
        (self.core / "metadata" / "analysis_metadata.json").write_text(
            json.dumps({"target_columns_consumed": []})
        )
        loaded = CoREMOFDataset.from_release(self.core)
        self.assertEqual(loaded.structure_ids, self.dataset.structure_ids)
        self.assertTrue(all(not set(TARGETS).intersection(row) for row in loaded.metadata_rows))

    def test_manifest_expected_hash_is_required_and_checked(self):
        for digest in ("0" * 64, "not-a-hash", "A" * 64):
            with self.subTest(digest=digest), self.assertRaises(TargetDataError):
                self.dataset.attach_target_supplement(self.supplement, expected_sha256=digest)

    def test_wrong_kind_schema_version_dataset_version_or_count(self):
        original = dict(self.manifest)
        for key, value in (
            ("kind", "coremof-core/1.0"), ("schema_version", "unknown/2.0"),
            ("dataset_version", "other"), ("structure_count", len(self.dataset) - 1),
            ("structure_count", True),
        ):
            with self.subTest(key=key, value=value):
                self.manifest = dict(original)
                self.manifest[key] = value
                self._refresh()
                self._reject()

    def test_incompatible_or_incomplete_core_bindings(self):
        self.manifest["compatible_core_sha256"][CORE_FILES[0]] = "0" * 64
        self._refresh()
        self._reject()
        self.manifest["compatible_core_sha256"].pop(CORE_FILES[0])
        self._refresh()
        self._reject()

    def test_corrupt_bound_file_and_wrong_byte_count(self):
        path = self.supplement / "targets.csv"
        path.write_bytes(path.read_bytes() + b"\n")
        self._reject()
        self._refresh()
        self.manifest["files"]["targets.csv"]["bytes"] += 1
        (self.supplement / "manifest.json").write_text(json.dumps(self.manifest))
        self.expected = _sha(self.supplement / "manifest.json")
        self._reject()

    def test_missing_required_file_binding(self):
        self.manifest["files"].pop("targets_schema.json")
        path = self.supplement / "manifest.json"
        path.write_text(json.dumps(self.manifest))
        self.expected = _sha(path)
        self._reject()

    def test_malformed_json_companions_are_rejected_even_when_hash_bound(self):
        path = self.supplement / "target_details.json"
        for data in ("not JSON", "[]"):
            with self.subTest(data=data):
                path.write_text(data)
                self._refresh()
                self._reject()

    def test_missing_or_extra_target_definitions_are_rejected(self):
        original = dict(self.manifest["target_definitions"])
        for definitions in (
            {name: definition for name, definition in original.items() if name != TARGETS[0]},
            dict(original, invented={"unit": "dimensionless", "conditions": {}}),
        ):
            with self.subTest(definitions=tuple(definitions)):
                self.manifest["target_definitions"] = definitions
                self._refresh()
                self._reject()

    def test_missing_and_duplicate_csv_headers_are_rejected(self):
        path = self.supplement / "targets.csv"
        original = path.read_bytes()
        lines = original.splitlines(keepends=True)
        for header in (
            lines[0].replace((TARGETS[0] + ",").encode(), b"", 1),
            lines[0].replace((TARGETS[1] + ",").encode(), (TARGETS[0] + ",").encode(), 1),
        ):
            with self.subTest(header=header):
                path.write_bytes(header + b"".join(lines[1:]))
                self._refresh()
                self._reject()

    def test_path_traversal_is_rejected(self):
        outside = self.root / "outside.json"
        outside.write_text("{}")
        self.manifest["files"]["../outside.json"] = {"sha256": _sha(outside), "bytes": 2}
        self._refresh()
        self._reject()

    def test_symlink_escape_is_rejected(self):
        outside = self.root / "outside.json"
        outside.write_text("{}")
        (self.supplement / "escape.json").symlink_to(outside)
        self.manifest["files"]["escape.json"] = {"sha256": _sha(outside), "bytes": 2}
        self._refresh()
        self._reject()

    def test_duplicate_csv_or_jsonl_ids(self):
        for name in ("targets.csv", "targets.jsonl"):
            with self.subTest(name=name):
                path = self.supplement / name
                original = path.read_bytes()
                lines = original.splitlines(keepends=True)
                path.write_bytes(original + lines[1 if name.endswith("csv") else 0])
                self._refresh()
                self._reject()
                path.write_bytes(original)
                self._refresh()

    def test_unknown_or_missing_ids_even_when_both_formats_agree(self):
        self.payloads[0]["structure_id"] = "2026[Cu][nan]3[ASR]999999"
        self._write_payloads()
        path = self.supplement / "targets.csv"
        fields, rows = _read_csv(path)
        rows[0]["structure_id"] = self.payloads[0]["structure_id"]
        _write_csv(path, fields, rows)
        self._refresh()
        self._reject()

    def test_csv_jsonl_value_or_status_disagreement(self):
        for column, value in ((TARGETS[0], "9"), (TARGETS[0] + "__status", "ERROR")):
            with self.subTest(column=column):
                path = self.supplement / "targets.csv"
                original = path.read_bytes()
                fields, rows = _read_csv(path)
                rows[0][column] = value
                _write_csv(path, fields, rows)
                self._refresh()
                self._reject()
                path.write_bytes(original)
                self._refresh()

    def test_bool_response_and_inconsistent_availability_or_status(self):
        for key, value in (("value", True), ("available", False), ("execution_status", "ERROR"), ("execution_status", "EXISTING")):
            with self.subTest(key=key, value=value):
                payload = self.payloads[0]["targets"][TARGETS[0]]
                original = payload[key]
                payload[key] = value
                self._write_payloads()
                self._refresh()
                self._reject()
                payload[key] = original

    def test_nan_or_infinity_csv_values(self):
        path = self.supplement / "targets.csv"
        fields, rows = _read_csv(path)
        for value in ("nan", "inf", "-inf"):
            with self.subTest(value=value):
                rows[0][TARGETS[0]] = value
                _write_csv(path, fields, rows)
                self._refresh()
                self._reject()

    def test_json_nonfinite_and_duplicate_keys(self):
        path = self.supplement / "targets.jsonl"
        original = path.read_text()
        for altered in (
            original.replace('"value": 0.0', '"value": NaN', 1),
            original.replace('"value": 0.0', '"value": 0.0, "value": 0.0', 1),
        ):
            with self.subTest(altered=altered[:40]):
                path.write_text(altered)
                self._refresh()
                self._reject()
        path.write_text(original)

    def test_overflowed_json_numbers_are_rejected_globally(self):
        for name, altered in (
            ("target_details.json", '{"reported_error": 1e999}'),
            ("targets.jsonl", (self.supplement / "targets.jsonl").read_text().replace(
                '"diagnostic": null', '"diagnostic": {"reported_error": 1e999}', 1
            )),
        ):
            with self.subTest(name=name):
                path = self.supplement / name
                original = path.read_text()
                path.write_text(altered)
                self._refresh()
                self._reject()
                path.write_text(original)
        self.manifest["diagnostic_note"] = "PLACEHOLDER"
        self._refresh()
        path = self.supplement / "manifest.json"
        path.write_text(path.read_text().replace('"PLACEHOLDER"', '1e999'))
        self.expected = _sha(path)
        self._reject()

    def test_unit_and_condition_mismatch(self):
        for key, value in (("unit", "invented"), ("conditions", {"temperature_K": 1000})):
            with self.subTest(key=key):
                payload = self.payloads[0]["targets"][TARGETS[0]]
                original = payload[key]
                payload[key] = value
                self._write_payloads()
                self._refresh()
                self._reject()
                payload[key] = original

    def test_matching_but_wrong_endpoint_definitions_are_rejected(self):
        original = json.loads(json.dumps(self.definitions))
        for key, value in (
            ("unit", " mol/kg-framework "),
            ("conditions", {}),
            ("conditions", {"temperature_K": 299, "pressure_Pa": 6500000}),
            ("conditions", {"temperature_K": True, "pressure_Pa": 6500000}),
            ("conditions", {"temperature_K": "298", "pressure_Pa": 6500000}),
        ):
            with self.subTest(key=key, value=value):
                self.manifest["target_definitions"] = json.loads(json.dumps(original))
                self.manifest["target_definitions"][TARGETS[0]][key] = value
                for row in self.payloads:
                    row["targets"][TARGETS[0]].update(json.loads(json.dumps(original[TARGETS[0]])))
                    row["targets"][TARGETS[0]][key] = value
                self._write_payloads()
                self._refresh()
                self._reject()

    def _projection(self, *, legacy=False):
        if legacy:
            metadata_path = self.core / "metadata" / "metadata.csv"
            fields, rows = _read_csv(metadata_path)
            for row in rows:
                for name in TARGETS:
                    row[name] = "0"
                    row[name + "__status"] = "SUCCESS"
            _write_csv(metadata_path, tuple(rows[0]), rows)
            info_path = self.core / "dataset_info.json"
            info = json.loads(info_path.read_text())
            info.pop("metadata_layout")
            info_path.write_text(json.dumps(info))
            complete = CoREMOFDataset.from_release(self.core, include_legacy_targets=True)
        else:
            complete = self.dataset
        source = self.root / "cod"
        selected = {row["structure_id"] for row in _metadata_rows() if row["source_database"] == "COD"}
        for name in ("metadata/metadata.csv", "parent_groups/parent_groups.csv", "manifests/cif_manifest.csv"):
            fields, rows = _read_csv(self.core / name)
            _write_csv(source / name, fields, [row for row in rows if row["structure_id"] in selected])
        methods = "parent_groups/parent_group_methods.json"
        (source / methods).write_bytes((self.core / methods).read_bytes())
        contract_path = source / "source_projection.json"
        exported = export_source_projection(
            complete, source, contract_path, sources=("COD",),
            group_profiles=(("rac5",),), diversity="none",
        )
        return source, contract_path, exported["sha256"]

    def test_legacy_source_projection_hides_targets_and_preserves_all_bindings(self):
        source, path, digest = self._projection(legacy=True)
        contract = json.loads(path.read_text())
        contract["dataset_info"]["target_metadata"] = {"coverage": 2}
        contract["dataset_info"]["metadata_revision"] = {
            "revision_id": "target_metadata_old", "targets_modified": True,
        }
        path.write_text(json.dumps(contract))
        digest = _sha(path)
        default = CoREMOFDataset.from_projection(source, path, expected_sha256=digest)
        legacy = CoREMOFDataset.from_projection(source, path, expected_sha256=digest, include_legacy_targets=True)
        self.assertTrue(all(not set(TARGETS).intersection(row) for row in default.metadata_rows))
        self.assertEqual(legacy.metadata_rows[0][TARGETS[0]], "0")
        self.assertNotIn("target_metadata", default.dataset_info)
        self.assertEqual(default.input_hashes, legacy.input_hashes)
        self.assertEqual(default.input_hashes["metadata/metadata.csv"], _sha(source / "metadata" / "metadata.csv"))
        self.assertEqual(default.parent_by_id, legacy.parent_by_id)
        self.assertEqual(default._authority_extra_state, legacy._authority_extra_state)

    def test_target_free_source_projection_rejects_embedded_columns_after_hash_binding(self):
        source, path, digest = self._projection()
        contract = json.loads(path.read_text())
        self.assertEqual(contract["dataset_info"]["metadata_layout"], "target_free_core/1.0")
        metadata_path = source / "metadata" / "metadata.csv"
        fields, rows = _read_csv(metadata_path)
        for row in rows:
            row[TARGETS[0]] = "0"
        _write_csv(metadata_path, tuple(rows[0]), rows)
        contract["selected_files"]["metadata/metadata.csv"] = {
            "sha256": _sha(metadata_path), "size_bytes": metadata_path.stat().st_size,
        }
        path.write_text(json.dumps(contract))
        digest = _sha(path)
        for flag in (False, True):
            with self.subTest(flag=flag), self.assertRaises(ReleaseValidationError):
                CoREMOFDataset.from_projection(source, path, expected_sha256=digest, include_legacy_targets=flag)

    def test_projection_target_stripping_never_bypasses_original_input_hashes(self):
        source, path, digest = self._projection(legacy=True)
        metadata_path = source / "metadata" / "metadata.csv"
        fields, rows = _read_csv(metadata_path)
        rows[0][TARGETS[0]] = "9"
        _write_csv(metadata_path, fields, rows)
        with self.assertRaises(ReleaseValidationError):
            CoREMOFDataset.from_projection(source, path, expected_sha256=digest)

    def test_projection_legacy_opt_in_requires_boolean(self):
        source, path, digest = self._projection()
        for value in (None, 0, "yes"):
            with self.subTest(value=value), self.assertRaises(TypeError):
                CoREMOFDataset.from_projection(source, path, expected_sha256=digest, include_legacy_targets=value)

    def test_repeated_attachment_is_rejected(self):
        attached = self._attach()
        with self.assertRaises(TargetDataError):
            attached.attach_target_supplement(self.supplement, expected_sha256=self.expected)

    def test_target_column_conflict_is_rejected(self):
        path = self.core / "metadata" / "metadata.csv"
        fields, rows = _read_csv(path)
        for row in rows:
            row[TARGETS[0]] = "0"
        _write_csv(path, tuple(rows[0]), rows)
        info_path = self.core / "dataset_info.json"
        info = json.loads(info_path.read_text())
        info.pop("metadata_layout")
        info_path.write_text(json.dumps(info))
        self.dataset = CoREMOFDataset.from_release(self.core, include_legacy_targets=True)
        self.manifest["compatible_core_sha256"] = {
            name: self.dataset.input_hashes[name] for name in CORE_FILES
        }
        self._refresh()
        with self.assertRaises(TargetDataError):
            self._attach()
        self.assertEqual(self.dataset.metadata_rows[0][TARGETS[0]], "0")

    def test_csv_changed_between_verification_and_merge_is_rejected(self):
        original_merge = supplements_module.merge_targets

        def raced_merge(dataset, sources, **kwargs):
            path = self.supplement / "targets.csv"
            fields, rows = _read_csv(path)
            rows[0][TARGETS[0]] = "2"
            _write_csv(path, fields, rows)
            return original_merge(dataset, sources, **kwargs)

        with patch.object(supplements_module, "merge_targets", raced_merge):
            self._reject()


if __name__ == "__main__":
    unittest.main()
