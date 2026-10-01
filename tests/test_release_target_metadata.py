"""Private release target enrichment, using disposable loader-valid fixtures."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from CoREMOF.dataset import CoREMOFDataset, ReleaseValidationError
from test_dataset_labels import _make_release, _metadata_rows


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("extend_release_metadata", ROOT / "examples/extend_release_metadata.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
try:
    import jsonschema
except ImportError:
    jsonschema = None


def target_row(sid="2000[Cu][nan]3[FSR]1", endpoint="ch4", value="0.0", status="SUCCESS"):
    name, unit, temperature, pressure = mod.ENDPOINTS[endpoint]
    return {"structure_id": sid, "endpoint": endpoint, "cif_sha256": "",
        "target_name": name, "unit": unit, "temperature_K": str(temperature),
        "pressure_Pa": "" if pressure is None else str(pressure),
        "assignment_status": status, "available": "true" if status == "SUCCESS" else "false",
        "value": value, "reported_error": "", "source_record_id": "PRIVATE-OLD-ID",
        "source_logical_name": "/private/worker/path"}


class TargetProjectionTests(unittest.TestCase):
    def test_zero_and_native_null_are_distinct(self):
        self.assertEqual(mod.public_target(target_row())["value"], 0.0)
        null = mod.public_target(target_row(endpoint="widom", value="", status="SUCCESS_NULL"))
        self.assertIsNone(null["value"])
        self.assertFalse(null["available"])
        self.assertEqual(null["diagnostic"]["code"], "UNDEFINED_SCIENTIFIC_RESULT")
        self.assertIsNone(null["conditions"]["pressure_Pa"])

    def test_rejects_nonfinite_inconsistent_and_changed_conditions(self):
        for field, value in (("value", "NaN"), ("value", "Infinity"), ("value", "true"),
                ("value", ""), ("reported_error", "-1"), ("unit", "cm3/g"),
                ("temperature_K", "300"), ("pressure_Pa", "1"), ("available", "false")):
            with self.subTest(field=field, value=value):
                row = target_row()
                row[field] = value
                with self.assertRaises(ValueError):
                    mod.public_target(row)

    def test_projection_does_not_export_internal_evidence(self):
        value = json.dumps(mod.public_target(target_row()))
        self.assertNotIn("PRIVATE", value)
        self.assertNotIn("/private", value)
        self.assertNotIn("cif_sha256", value)

    def test_recorded_method_fields_survive_but_private_values_fail(self):
        row = target_row()
        row.update(protocol_charge_method="Ewald", raspa_version="3.1.0",
            cif_preprocessing_action="UNCHANGED_COMPLETE_CHARGES")
        self.assertEqual(mod.public_target(row)["calculation_method"]["protocol_charge_method"], "Ewald")
        self.assertEqual(mod.public_target(row)["value"], 0.0)
        row["raspa_version"] = "/home/private/raspa"
        with self.assertRaises(ValueError):
            mod.public_target(row)
        with self.assertRaises(ValueError):
            mod.number(True)

    def test_unsafe_ledger_paths_and_duplicate_json_fail(self):
        for value in ("../outside", "/etc/passwd", "a/../b", "a\\b", "a//b"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                mod.relative_path(value)
        with self.assertRaises(ValueError):
            mod.strict_json('{"value": 0, "value": 1}')


@unittest.skipUnless(jsonschema is not None, "jsonschema is required for full record validation")
class ReleaseTargetMetadataTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / "source"
        self.targets = self.root / "targets"
        self.source.mkdir()
        self.targets.mkdir()
        _make_release(self.source)
        metadata = _metadata_rows()
        self.ids = [row["structure_id"] for row in metadata]
        cif_data = {sid: ("data_" + sid + "\n_cell_length_a 10.0\n").encode() for sid in self.ids}
        info = json.loads((self.source / "dataset_info.json").read_text())
        info.update({"release_status": "PROVISIONAL_LATEST_AUDITED_SNAPSHOT",
            "metadata_revision": {"revision_id": "source", "targets_or_frozen_splits_modified": False},
            "per_structure_json": {"schema_version": "coremof-structure-record/1.0"},
            "tabular_files": {"metadata/metadata.csv": {"row_count": len(self.ids)}}})
        schema = {"type": "object", "additionalProperties": False,
            "required": ["schema_version", "structure_id", "identity"],
            "properties": {"schema_version": {"const": "coremof-structure-record/1.0"},
                "structure_id": {"type": "string"}, "identity": {"type": "object"}}}
        for row in metadata:
            row["same_feature"] = "1.23456789012345"
        self.records = {sid: {"schema_version": "coremof-structure-record/1.0", "structure_id": sid,
            "identity": {"source_database": "COD", "mofid": "SOURCE_UNCHANGED"}} for sid in self.ids}
        payloads = {
            "dataset_info.json": mod.canonical(info), "README.md": b"Source method documentation.\n",
            "metadata/structure_record_schema.json": mod.canonical(schema),
            "metadata/analysis_metadata_schema.json": mod.canonical({"formats": {
                "metadata/structure_annotations.jsonl": "joins structure-record/1.0 without changing it"},
                "target_columns_consumed": []}),
            "metadata/metadata.csv": mod.csv_bytes(list(metadata[0]), metadata),
            "metadata/metadata.jsonl": ("\n".join(json.dumps(row) for row in metadata) + "\n").encode(),
            "manifests/cif_manifest.csv": mod.csv_bytes(["structure_id", "cif_file", "size_bytes", "sha256"], [
                {"structure_id": sid, "cif_file": f"cifs/{sid}.cif", "size_bytes": len(data),
                 "sha256": mod.digest(data)} for sid, data in cif_data.items()]),
            "manifests/structure_json_manifest.csv": b"structure_id,json_file,size_bytes,sha256\n",
            "manifests/validation.json": b'{"validation_status":"PASS"}',
            "manifests/metadata_revision_validation.json": b'{"validation_status":"PASS"}',
            "parent_groups/parent_groups.csv": (self.source / "parent_groups/parent_groups.csv").read_bytes(),
            "parent_groups/parent_group_methods.json": (self.source / "parent_groups/parent_group_methods.json").read_bytes(),
        }
        payloads.update({f"cifs/{sid}.cif": data for sid, data in cif_data.items()})
        payloads.update({f"metadata/structures/{sid}.json": mod.canonical(record) for sid, record in self.records.items()})
        self.payloads = payloads
        for name, data in payloads.items():
            path = self.source / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        self.source_hash = self.seal(self.source, payloads)
        self.rows = []
        for sid in self.ids:
            for endpoint in mod.ENDPOINTS:
                row = target_row(sid, endpoint)
                row["cif_sha256"] = mod.digest(cif_data[sid])
                if sid == self.ids[1] and endpoint == "widom":
                    row.update(value="", assignment_status="SUCCESS_NULL", available="false")
                self.rows.append(row)
        self.reseal_targets()

    @staticmethod
    def seal(root, payloads):
        data = "".join(f"{mod.digest(data)}  {name}\n" for name, data in sorted(payloads.items())).encode()
        (root / "SHA256SUMS").write_bytes(data)
        return mod.digest(data)

    def reseal_targets(self):
        data = mod.csv_bytes(list(self.rows[0]), self.rows)
        (self.targets / "target_assignments.csv").write_bytes(data)
        self.target_hash = self.seal(self.targets, {"target_assignments.csv": data})

    def build(self):
        return mod.build(self.source, self.targets, self.root / "output", self.source_hash, self.target_hash)

    def test_complete_output_preserves_inputs_and_agrees_in_every_format(self):
        receipt = self.build()
        output = self.root / "output"
        self.assertFalse(receipt["publication_authorized"])
        self.assertEqual(receipt["per_structure_schema_validations"], len(self.ids))
        loaded = CoREMOFDataset.from_release(output, verify_cif_files=True)
        self.assertEqual(len(loaded), len(self.ids))
        self.assertTrue(loaded.cif_files_verified)
        info = json.loads((output / "dataset_info.json").read_text())
        self.assertEqual(info["release_status"], "PROVISIONAL_LATEST_AUDITED_SNAPSHOT")
        self.assertNotIn("official_split", info)
        self.assertNotIn("publication_authorized", info)
        binding = json.loads((output / "manifests/target_metadata_inputs.json").read_text())
        self.assertEqual(binding["accepted_targets_ledger_sha256"], self.target_hash)
        analysis = json.loads((output / "metadata/analysis_metadata_schema.json").read_text())
        self.assertIn("structure-record/1.1", analysis["formats"]["metadata/structure_annotations.jsonl"])
        self.assertEqual(analysis["target_columns_consumed"], [])
        records = {}
        for sid in self.ids:
            records[sid] = json.loads((output / f"metadata/structures/{sid}.json").read_text())
            original = {key: value for key, value in records[sid].items() if key != "targets"}
            original["schema_version"] = "coremof-structure-record/1.0"
            self.assertEqual(original, self.records[sid])
        for name, data in self.payloads.items():
            self.assertEqual((self.source / name).read_bytes(), data)
            if name.startswith(("cifs/", "parent_groups/")):
                self.assertEqual((output / name).read_bytes(), data)
        for name in ("metadata/targets.jsonl", "metadata/metadata.jsonl"):
            for line in (output / name).read_text().splitlines():
                row = json.loads(line)
                self.assertEqual(row["targets"], records[row["structure_id"]]["targets"])
        with (output / "metadata/metadata.csv").open() as stream:
            for row in csv.DictReader(stream):
                for name, target in records[row["structure_id"]]["targets"].items():
                    self.assertEqual(row[name], "" if target["value"] is None else repr(target["value"]))
        for name, expected in mod.ledger(output, receipt["output_ledger_sha256"]).items():
            mod.read_bound(output, name, expected)

    def test_loader_rejection_prevents_publishing_and_preserves_source(self):
        # False-valued authority declarations are reserved too. Never erase
        # or weaken the loader's evidence boundary to accept an export.
        for key, value in (("official_split", False), ("publication_authorized", False),
                ("release_status", "STAGED_CANDIDATE_NOT_PUBLISHED")):
            with self.subTest(key=key):
                altered = dict(self.payloads)
                info = json.loads(altered["dataset_info.json"])
                info[key] = value
                altered["dataset_info.json"] = mod.canonical(info)
                (self.source / "dataset_info.json").write_bytes(altered["dataset_info.json"])
                self.source_hash = self.seal(self.source, altered)
                with self.assertRaises(ReleaseValidationError):
                    self.build()
                self.assertFalse((self.root / "output").exists())
                self.assertEqual((self.source / "dataset_info.json").read_bytes(), altered["dataset_info.json"])

    def test_generated_invalid_metadata_is_rejected_before_output_publication(self):
        serialize = mod.canonical

        def corrupt_generated_info(value):
            if isinstance(value, dict) and "target_metadata" in value:
                value = {**value, "official_split": False}
            return serialize(value)

        original_info = (self.source / "dataset_info.json").read_bytes()
        with patch.object(mod, "canonical", side_effect=corrupt_generated_info):
            with self.assertRaises(ReleaseValidationError):
                self.build()
        self.assertFalse((self.root / "output").exists())
        self.assertEqual((self.source / "dataset_info.json").read_bytes(), original_info)

    def test_corrupt_source_and_existing_output_are_not_overwritten(self):
        (self.source / f"cifs/{self.ids[0]}.cif").write_bytes(b"corrupt")
        with self.assertRaises(ValueError):
            self.build()
        self.assertFalse((self.root / "output").exists())
        (self.root / "output").mkdir()
        sentinel = self.root / "output/keep.txt"
        sentinel.write_text("keep")
        with self.assertRaises(FileExistsError):
            self.build()
        self.assertEqual(sentinel.read_text(), "keep")

    def test_missing_duplicate_and_mismatched_structure_targets_fail(self):
        original = list(self.rows)
        for variant in (original[:-1], original + [dict(original[0])],
                [{**row, "cif_sha256": "0" * 64} if i == 0 else row for i, row in enumerate(original)]):
            with self.subTest(rows=len(variant)):
                self.rows = variant
                self.reseal_targets()
                with self.assertRaises(ValueError):
                    self.build()
                self.assertFalse((self.root / "output").exists())

    def test_missing_or_duplicated_metadata_table_rows_fail(self):
        variants = [
            ("metadata/metadata.csv", mod.csv_bytes(["structure_id"], [{"structure_id": self.ids[0]}] * 2)),
            ("metadata/metadata.jsonl", (json.dumps({"structure_id": self.ids[0]}) + "\n").encode()),
            ("metadata/metadata.jsonl", ((json.dumps({"structure_id": self.ids[0]}) + "\n") * 2).encode()),
        ]
        for name, bad in variants:
            with self.subTest(file=name, size=len(bad)):
                altered = dict(self.payloads)
                altered[name] = bad
                for path, data in altered.items():
                    (self.source / path).write_bytes(data)
                self.source_hash = self.seal(self.source, altered)
                with self.assertRaises(ValueError):
                    self.build()
                self.assertFalse((self.root / "output").exists())

    def test_symlinked_output_parent_cannot_modify_inputs(self):
        alias = self.root / "alias"
        alias.symlink_to(self.source, target_is_directory=True)
        with self.assertRaises(ValueError):
            mod.build(self.source, self.targets, alias / "new", self.source_hash, self.target_hash)
        self.assertFalse((self.source / "new").exists())


if __name__ == "__main__":
    unittest.main()
