import csv
import hashlib
import json
from pathlib import Path
import tempfile
from types import MappingProxyType
import unittest
from unittest.mock import patch

from CoREMOF import __version__
import CoREMOF.splitters as splitters_module
import CoREMOF.targets as targets_module
from CoREMOF.dataset import CoREMOFDataset, StructureRecord
from CoREMOF.targets import (
    AliasRegistry,
    TARGET_API_VERSION,
    TargetDataError,
    TargetSource,
    merge_targets,
    merge_targets_from_config,
)


IDS = (
    "ASR-COD-2026-0001",
    "ASR-CSD-2026-0001",
    "FSR-COD-2026-0001",
    "ION-SI-2026-0001",
)


def _write_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _dataset(root):
    records = []
    hashes = {}
    for index, structure_id in enumerate(IDS):
        variant, source, _, _ = structure_id.split("-")
        # The hidden CSD row bridges the first and third rows through exact
        # hash and source evidence in the full-universe main_union graph.
        source_id = "BRIDGE" if index in (1, 2) else "SRC-{}".format(index)
        sha = hashlib.sha256(structure_id.encode("utf-8")).hexdigest()
        if index in (0, 1):
            sha = "a" * 64
        hashes[structure_id] = sha
        metadata = {
            "structure_id": structure_id,
            "cif_file": "cifs/{}.cif".format(structure_id),
            "source_database": source,
            "source_id": source_id,
            "structure_variant": variant,
            "metal_elements": "Cu",
            "mofclassifier_status": "PASS",
            "mofchecker_status": "PASS",
            "chen_manz_status": "PASS",
            "mosaec_status": "PASS",
            "setc_gat_status": "PASS",
            "label_3checker": "CR",
            "label_4checker": "CR",
            "label_5checker": "CR",
        }
        records.append(
            StructureRecord(
                structure_id=structure_id,
                metadata=MappingProxyType(metadata),
                parent_groups=MappingProxyType({}),
                cif_manifest=MappingProxyType(
                    {
                        "structure_id": structure_id,
                        "cif_file": metadata["cif_file"],
                        "size_bytes": "1",
                        "sha256": sha,
                    }
                ),
            )
        )
    parent_by_id = {
        IDS[1]: {
            "mofid2_status": "MATCHED",
            "mofid2_group": "M2-ABCD0001",
            "mofid2_size": "2",
        },
        IDS[2]: {
            "mofid2_status": "MATCHED",
            "mofid2_group": "M2-ABCD0001",
            "mofid2_size": "2",
        },
    }
    dataset = CoREMOFDataset(
        release_root=root,
        records=records,
        dataset_info={"dataset_version": "vtest", "release_status": "FINAL"},
        parent_group_methods={"release_status": "FINAL"},
        parent_by_id=parent_by_id,
        input_hashes={"metadata/metadata.csv": "b" * 64},
        cif_files_verified=False,
    )
    feature_rows = [
        {"structure_id": structure_id, "rac5_available": "true", "rac_x": str(index)}
        for index, structure_id in enumerate(IDS)
    ]
    _write_csv(
        root / "features" / "rac5_features.csv",
        ("structure_id", "rac5_available", "rac_x"),
        feature_rows,
    )
    return dataset


class TargetMergeTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.dataset = _dataset(self.root)
        self.alias_path = self.root / "aliases.csv"
        _write_csv(
            self.alias_path,
            ("structure_id", "earlier_id"),
            [
                {"structure_id": structure_id, "earlier_id": "OLD-{}".format(index)}
                for index, structure_id in enumerate(IDS)
            ],
        )
        self.aliases = AliasRegistry(
            self.alias_path, alias_columns=("earlier_id",)
        )

    def tearDown(self):
        self.temporary_directory.cleanup()

    def _sources(self):
        uptake = self.root / "uptake.csv"
        _write_csv(
            uptake,
            ("name", "uptake"),
            [
                {"name": "OLD-0", "uptake": "1.5"},
                {"name": IDS[1], "uptake": ""},
                {"name": IDS[2], "uptake": "2.0"},
            ],
        )
        selectivity = self.root / "selectivity.jsonl"
        selectivity.write_text(
            "\n".join(
                (
                    json.dumps({"structure": IDS[0], "selectivity": 10.0}),
                    json.dumps({"structure": "OLD-2", "selectivity": 20.0}),
                    json.dumps({"structure": IDS[3], "selectivity": None}),
                )
            )
            + "\n",
            encoding="utf-8",
        )
        return (
            TargetSource(
                uptake,
                id_column="name",
                target_columns=("uptake",),
                value_types={"uptake": "float"},
                units={"uptake": "mol/kg"},
                conditions={"uptake": {"temperature_K": 298, "pressure_bar": 1}},
            ),
            TargetSource(
                selectivity,
                id_column="structure",
                target_columns=("selectivity",),
                units={"selectivity": "dimensionless"},
                conditions={
                    "selectivity": {
                        "temperature_K": 298,
                        "mixture": {"Xe": 0.2, "Kr": 0.8},
                    }
                },
            ),
        )

    def test_multiple_files_aliases_features_nulls_and_provenance(self):
        merged = merge_targets(
            self.dataset,
            self._sources(),
            alias_registry=self.aliases,
            feature_tables=("rac5",),
        )
        self.assertEqual(merged.target_columns, ("uptake", "selectivity"))
        self.assertIn("rac_x", merged.feature_columns)
        self.assertEqual(merged[IDS[0]].metadata["uptake"], 1.5)
        self.assertEqual(merged[IDS[0]].metadata["selectivity"], 10.0)
        self.assertIsNone(merged[IDS[1]].metadata["uptake"])
        self.assertIsNone(merged[IDS[1]].metadata["selectivity"])
        self.assertEqual(merged[IDS[2]].metadata["rac_x"], "2")
        self.assertIsNone(merged.target_values(IDS[3])["selectivity"])
        observation = merged.target_provenance_by_id[IDS[0]]["uptake"][0]
        self.assertEqual(observation["id_resolution"], "ALIAS")
        self.assertEqual(observation["input_structure_id"], "OLD-0")
        receipt = merged.receipt()
        self.assertEqual(receipt["counts"]["source_files"], 2)
        self.assertEqual(receipt["counts"]["alias_id_rows"], 2)
        self.assertEqual(receipt["target_definitions"]["uptake"]["unit"], "mol/kg")
        self.assertFalse(receipt["policies"]["units_inferred"])
        self.assertIn("targets/uptake", merged.input_hashes)
        with self.assertRaises(TypeError):
            merged.target_values(IDS[0])["uptake"] = 99

    def test_json_targets_are_canonical_immutable_and_order_independent(self):
        first = self.root / "profile.csv"
        second = self.root / "profile.jsonl"
        _write_csv(
            first,
            ("structure_id", "profile"),
            [
                {
                    "structure_id": IDS[0],
                    "profile": '{"b":[2,1],"a":{"x":1}}',
                }
            ],
        )
        second.write_text(
            json.dumps(
                {
                    "structure_id": IDS[0],
                    "profile": {"a": {"x": 1}, "b": [2, 1]},
                }
            )
            + "\n",
            encoding="utf-8",
        )
        common = {
            "target_columns": ("profile",),
            "value_types": {"profile": "json"},
            "units": {"profile": "user-defined"},
        }
        merged = merge_targets(
            self.dataset,
            (
                TargetSource(first, name="profile_csv", **common),
                TargetSource(second, name="profile_jsonl", **common),
            ),
        )
        profile = merged[IDS[0]].metadata["profile"]
        self.assertEqual(dict(profile["a"]), {"x": 1})
        self.assertEqual(profile["b"], (2, 1))
        self.assertEqual(
            len(merged.target_provenance_by_id[IDS[0]]["profile"]), 2
        )
        with self.assertRaises(TypeError):
            profile["a"]["x"] = 2
        with self.assertRaises(TypeError):
            profile["b"][0] = 3

        written = merged.write(self.root / "json-target-output")[0]
        with written.open("r", encoding="utf-8", newline="") as handle:
            row = next(csv.DictReader(handle))
        self.assertEqual(row["profile"], '{"a":{"x":1},"b":[2,1]}')

    def test_json_target_conflicts_and_nested_nonfinite_values_fail_closed(self):
        first = self.root / "first-profile.csv"
        second = self.root / "second-profile.csv"
        _write_csv(
            first,
            ("structure_id", "profile"),
            [{"structure_id": IDS[0], "profile": '{"score":1}'}],
        )
        _write_csv(
            second,
            ("structure_id", "profile"),
            [{"structure_id": IDS[0], "profile": '{"score":2}'}],
        )
        sources = tuple(
            TargetSource(
                path,
                name="profile_{}".format(index),
                target_columns=("profile",),
                value_types={"profile": "json"},
            )
            for index, path in enumerate((first, second), start=1)
        )
        with self.assertRaisesRegex(TargetDataError, "conflicting duplicate target"):
            merge_targets(self.dataset, sources)

        _write_csv(
            second,
            ("structure_id", "profile"),
            [{"structure_id": IDS[1], "profile": '{"score":NaN}'}],
        )
        with self.assertRaisesRegex(TargetDataError, "not finite JSON data"):
            merge_targets(self.dataset, sources)
        with self.assertRaisesRegex(TargetDataError, "not finite JSON data"):
            TargetSource(
                first,
                target_columns=("profile",),
                conditions={"profile": {"not_json": {1, 2}}},
            )

    def test_merge_receipt_binds_deterministic_implementation_sources(self):
        merged = merge_targets(
            self.dataset,
            self._sources(),
            alias_registry=self.aliases,
            feature_tables=("rac5",),
        )
        receipt = merged.receipt()
        implementation = receipt["implementation"]
        self.assertEqual(implementation["package"], "CoREMOF-tools")
        self.assertEqual(implementation["package_version"], __version__)
        self.assertEqual(
            implementation["target_api_version"], TARGET_API_VERSION
        )
        self.assertEqual(receipt["target_api_version"], TARGET_API_VERSION)

        package_root = Path(targets_module.__file__).resolve().parent
        expected = {
            filename: hashlib.sha256(
                (package_root / filename).read_bytes()
            ).hexdigest()
            for filename in ("dataset.py", "labels.py", "targets.py")
        }
        self.assertEqual(implementation["source_sha256"], expected)
        self.assertTrue(
            all(len(digest) == 64 for digest in expected.values())
        )

        written_receipt = merged.write(self.root / "implementation-receipt")[2]
        self.assertEqual(
            json.loads(written_receipt.read_text(encoding="utf-8"))[
                "implementation"
            ],
            implementation,
        )

    def test_merge_fails_if_imported_implementation_sources_drift(self):
        changed = dict(targets_module._IMPORTED_IMPLEMENTATION_HASHES)
        changed["targets.py"] = "0" * 64
        with patch.object(
            targets_module,
            "_current_implementation_hashes",
            return_value=changed,
        ):
            with self.assertRaisesRegex(
                TargetDataError, "implementation source changed after module import"
            ):
                targets_module._implementation_hashes()

    def test_file_snapshot_rejects_replacement_and_in_place_mutation(self):
        source = self.root / "snapshot.csv"
        source.write_text("structure_id,target\nA,1\n", encoding="utf-8")

        def replace_during_read(handle):
            data = handle.read()
            replacement = self.root / "snapshot.replacement"
            replacement.write_text(
                "structure_id,target\nB,2\n", encoding="utf-8"
            )
            replacement.replace(source)
            return data

        with patch.object(
            targets_module,
            "_read_snapshot_bytes",
            side_effect=replace_during_read,
        ):
            with self.assertRaisesRegex(TargetDataError, "changed or was replaced"):
                targets_module._capture_file(source)

        source.write_text("structure_id,target\nA,1\n", encoding="utf-8")

        def mutate_during_read(handle):
            data = handle.read()
            with source.open("r+b") as writer:
                writer.seek(0)
                writer.write(b"X")
                writer.flush()
                writer.seek(0)
            return data

        with patch.object(
            targets_module,
            "_read_snapshot_bytes",
            side_effect=mutate_during_read,
        ):
            with self.assertRaisesRegex(TargetDataError, "changed or was replaced"):
                targets_module._capture_file(source)

    def test_parsing_and_receipts_use_the_captured_byte_generation(self):
        source = self.root / "captured.csv"
        old_bytes = b"structure_id,target\nASR-COD-2026-0001,1\n"
        source.write_bytes(old_bytes)
        snapshot = targets_module._capture_file(source)
        source.write_text(
            "structure_id,target\nASR-COD-2026-0001,999\n", encoding="utf-8"
        )
        fields, rows = targets_module._read_csv_records(snapshot)
        self.assertEqual(fields, ("structure_id", "target"))
        self.assertEqual(rows[0][1]["target"], "1")
        self.assertEqual(snapshot.size_bytes, len(old_bytes))
        self.assertEqual(snapshot.sha256, hashlib.sha256(old_bytes).hexdigest())

    def test_all_current_feature_tables_join_without_column_collisions(self):
        feature_specs = {
            "zeo_features.csv": ("zeo_available", "n2_volume"),
            "zeo_zero_probe_features.csv": ("zero_probe_available", "zero_volume"),
            "topology_features.csv": ("topology_available", "network_dimension"),
        }
        for filename, columns in feature_specs.items():
            _write_csv(
                self.root / "features" / filename,
                ("structure_id",) + columns,
                [
                    {
                        "structure_id": structure_id,
                        columns[0]: "true",
                        columns[1]: str(index),
                    }
                    for index, structure_id in enumerate(IDS)
                ],
            )
        merged = merge_targets(
            self.dataset,
            (self._sources()[0],),
            alias_registry=self.aliases,
            feature_tables=("rac5", "zeo", "zeo_zero_probe", "topology"),
        )
        self.assertEqual(len(merged.feature_columns), 8)
        self.assertEqual(merged[IDS[2]].metadata["rac_x"], "2")
        self.assertEqual(merged[IDS[2]].metadata["n2_volume"], "2")
        self.assertEqual(merged[IDS[2]].metadata["zero_volume"], "2")
        self.assertEqual(merged[IDS[2]].metadata["network_dimension"], "2")

    def test_required_targets_filter_before_assignment_keeps_full_universe_bridge(self):
        merged = merge_targets(
            self.dataset,
            self._sources(),
            alias_registry=self.aliases,
        )
        split = merged.classify("5checker").train_valid_test_split(
            parent_method="priority_main",
            leakage_guard="main_union",
            labels=("CR",),
            sources=("COD",),
            required_targets=("uptake", "selectivity"),
            required_target_mode="all",
            fractions=(0.5, 0.0, 0.5),
            random_state=7,
        )
        self.assertEqual(set(split.assignments), {IDS[0], IDS[2]})
        self.assertEqual(split.assignments[IDS[0]], split.assignments[IDS[2]])
        self.assertEqual(split.leakage_groups[IDS[0]], split.leakage_groups[IDS[2]])
        self.assertEqual(split.exclusions[IDS[1]], "MISSING_REQUIRED_TARGET")
        self.assertTrue(split.leakage_audit["passed"])
        target_filter = split.receipt()["filters"]["targets"]
        self.assertTrue(target_filter["filter_precedes_assignment"])
        self.assertTrue(target_filter["leakage_blocks_use_full_release_universe"])
        self.assertEqual(target_filter["eligible_release_count"], 2)
        self.assertEqual(
            split.receipt()["target_data"]["target_values_sha256"],
            merged.receipt()["target_values_sha256"],
        )
        merge_source_hashes = merged.receipt()["implementation"]["source_sha256"]
        split_source_hashes = split.receipt()["implementation"]["source_sha256"]
        self.assertEqual(
            split_source_hashes["targets.py"], merge_source_hashes["targets.py"]
        )
        for filename in ("dataset.py", "labels.py"):
            self.assertEqual(split_source_hashes[filename], merge_source_hashes[filename])

        tampered = json.loads(json.dumps(merged.receipt()))
        tampered["implementation"]["source_sha256"]["dataset.py"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "changed between target merge"):
            splitters_module._implementation_hashes(
                include_targets=True, target_receipt=tampered
            )

    def test_any_target_mode_and_unknown_requirements(self):
        merged = merge_targets(
            self.dataset, self._sources(), alias_registry=self.aliases
        )
        split = merged.classify("5checker").train_valid_test_split(
            parent_method="none",
            leakage_guard="parent_only",
            labels=("CR",),
            required_targets=("uptake", "selectivity"),
            required_target_mode="any",
            fractions=(1.0, 0.0, 0.0),
        )
        self.assertEqual(set(split.assignments), {IDS[0], IDS[2]})
        with self.assertRaisesRegex(ValueError, "unknown required target"):
            merged.classify("5checker").train_valid_test_split(
                parent_method="none", required_targets=("not_a_target",)
            )
        with self.assertRaisesRegex(ValueError, "created by merge_targets"):
            self.dataset.classify("5checker").train_valid_test_split(
                parent_method="none", required_targets=("uptake",)
            )

    def test_unknown_and_ambiguous_aliases_fail_closed(self):
        source = self.root / "unknown.csv"
        _write_csv(source, ("structure_id", "target"), [{"structure_id": "OLD-0", "target": "1"}])
        with self.assertRaisesRegex(TargetDataError, "neither a current public ID"):
            merge_targets(self.dataset, (source,))

        _write_csv(
            self.alias_path,
            ("structure_id", "earlier_id"),
            [
                {"structure_id": IDS[0], "earlier_id": "SAME"},
                {"structure_id": IDS[1], "earlier_id": "SAME"},
            ],
        )
        with self.assertRaisesRegex(TargetDataError, "ambiguous alias"):
            merge_targets(
                self.dataset,
                (source,),
                alias_registry=AliasRegistry(
                    self.alias_path, alias_columns=("earlier_id",)
                ),
            )

    def test_conflicting_duplicates_and_definitions_fail_closed(self):
        first = self.root / "first.json"
        second = self.root / "second.json"
        first.write_text(
            json.dumps([{"structure_id": IDS[0], "uptake": 1.0}]),
            encoding="utf-8",
        )
        second.write_text(
            json.dumps([{"structure_id": IDS[0], "uptake": 2.0}]),
            encoding="utf-8",
        )
        with self.assertRaisesRegex(TargetDataError, "conflicting duplicate target"):
            merge_targets(
                self.dataset,
                (
                    TargetSource(first, units={"uptake": "mol/kg"}),
                    TargetSource(second, units={"uptake": "mol/kg"}),
                ),
            )
        second.write_text(
            json.dumps([{"structure_id": IDS[1], "uptake": 1.0}]),
            encoding="utf-8",
        )
        with self.assertRaisesRegex(TargetDataError, "conflicting unit/condition"):
            merge_targets(
                self.dataset,
                (
                    TargetSource(first, units={"uptake": "mol/kg"}),
                    TargetSource(second, units={"uptake": "cm3/g"}),
                ),
            )

    def test_generic_input_columns_can_be_renamed_to_canonical_targets(self):
        xe = self.root / "xe.csv"
        kr = self.root / "kr.csv"
        _write_csv(
            xe,
            ("structure_id", "uptake"),
            [{"structure_id": IDS[0], "uptake": "1.25"}],
        )
        _write_csv(
            kr,
            ("structure_id", "uptake"),
            [{"structure_id": IDS[0], "uptake": "2.50"}],
        )
        merged = merge_targets(
            self.dataset,
            (
                TargetSource(
                    xe,
                    name="xe",
                    target_columns=("uptake",),
                    target_names={"uptake": "xe_uptake"},
                    value_types={"xe_uptake": "float"},
                    units={"xe_uptake": "mol/kg"},
                    conditions={"xe_uptake": {"temperature_K": 298}},
                ),
                TargetSource(
                    kr,
                    name="kr",
                    target_columns=("uptake",),
                    target_names={"uptake": "kr_uptake"},
                    value_types={"kr_uptake": "float"},
                    units={"kr_uptake": "mol/kg"},
                    conditions={"kr_uptake": {"temperature_K": 298}},
                ),
            ),
        )
        self.assertEqual(merged.target_columns, ("xe_uptake", "kr_uptake"))
        self.assertEqual(merged.target_values(IDS[0])["xe_uptake"], 1.25)
        self.assertEqual(merged.target_values(IDS[0])["kr_uptake"], 2.5)
        xe_provenance = merged.target_provenance_by_id[IDS[0]]["xe_uptake"][0]
        self.assertEqual(xe_provenance["input_column"], "uptake")
        source_receipt = merged.receipt()["sources"][0]
        self.assertEqual(source_receipt["target_names"], {"uptake": "xe_uptake"})

    def test_output_is_deterministic_and_keeps_current_ids(self):
        merged = merge_targets(
            self.dataset, self._sources(), alias_registry=self.aliases
        )
        first = self.root / "out1"
        second = self.root / "out2"
        paths1 = merged.write(first)
        paths2 = merged.write(second)
        self.assertEqual([path.read_bytes() for path in paths1], [path.read_bytes() for path in paths2])
        with paths1[0].open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual([row["structure_id"] for row in rows], list(IDS))
        self.assertEqual(len(rows), len(IDS))
        self.assertNotIn("earlier_id", rows[0])

    def test_nonoverwrite_rollback_preserves_a_concurrent_replacement(self):
        merged = merge_targets(
            self.dataset, self._sources(), alias_registry=self.aliases
        )
        output = self.root / "replacement-race"
        real_link = targets_module.os.link
        publication_count = 0
        concurrent_bytes = b"concurrent writer generation\n"

        def replace_then_fail(source, target):
            nonlocal publication_count
            destination = Path(target)
            if destination.parent == output:
                publication_count += 1
            if publication_count == 1 and destination.parent == output:
                real_link(source, target)
                replacement = output / ".concurrent-replacement"
                replacement.write_bytes(concurrent_bytes)
                targets_module.os.replace(str(replacement), str(target))
                return None
            if publication_count == 2 and destination.parent == output:
                raise OSError("deterministic second publication failure")
            return real_link(source, target)

        with patch.object(targets_module.os, "link", side_effect=replace_then_fail):
            with self.assertRaisesRegex(
                OSError, "deterministic second publication failure"
            ):
                merged.write(output, stem="race", overwrite=False)

        self.assertEqual((output / "race.csv").read_bytes(), concurrent_bytes)
        self.assertFalse((output / "race.provenance.jsonl").exists())
        self.assertFalse((output / "race.json").exists())
        self.assertFalse(any(path.name.startswith(".race-") for path in output.iterdir()))

    def test_nonoverwrite_rollback_cannot_unlink_a_post_identity_replacement(self):
        merged = merge_targets(
            self.dataset, self._sources(), alias_registry=self.aliases
        )
        output = self.root / "post-check-replacement-race"
        target = output / "race.csv"
        real_link = targets_module.os.link
        real_replace = targets_module.os.replace
        real_samefile = targets_module.os.path.samefile
        publication_count = 0
        identity_checked = False
        concurrent_bytes = b"replacement installed after identity check\n"

        def fail_second_publication(source, destination):
            nonlocal publication_count
            publication_count += 1
            if publication_count == 2:
                raise OSError("deterministic second publication failure")
            return real_link(source, destination)

        def replace_final_after_identity_check(left, right):
            nonlocal identity_checked
            result = real_samefile(left, right)
            if result and not identity_checked:
                identity_checked = True
                replacement = output / ".post-check-replacement"
                replacement.write_bytes(concurrent_bytes)
                real_replace(str(replacement), str(target))
            return result

        with patch.object(
            targets_module.os, "link", side_effect=fail_second_publication
        ), patch.object(
            targets_module.os.path,
            "samefile",
            side_effect=replace_final_after_identity_check,
        ):
            with self.assertRaisesRegex(
                OSError, "deterministic second publication failure"
            ):
                merged.write(output, stem="race", overwrite=False)

        self.assertTrue(identity_checked)
        self.assertEqual(target.read_bytes(), concurrent_bytes)
        self.assertFalse((output / "race.provenance.jsonl").exists())
        self.assertFalse((output / "race.json").exists())
        self.assertFalse(any(path.name.startswith(".race-") for path in output.iterdir()))

    def test_overwrite_rollback_removes_new_output_and_restores_existing_outputs(self):
        merged = merge_targets(
            self.dataset, self._sources(), alias_registry=self.aliases
        )
        output = self.root / "overwrite-publication-failure"
        output.mkdir()
        provenance = output / "race.provenance.jsonl"
        receipt = output / "race.json"
        old_provenance = b"previous provenance generation\n"
        old_receipt = b'{"previous": true}\n'
        provenance.write_bytes(old_provenance)
        receipt.write_bytes(old_receipt)
        real_replace = targets_module.os.replace
        replacement_count = 0

        def fail_second_replacement(source, target):
            nonlocal replacement_count
            replacement_count += 1
            if replacement_count == 2:
                raise OSError("deterministic second replacement failure")
            return real_replace(source, target)

        with patch.object(
            targets_module.os, "replace", side_effect=fail_second_replacement
        ):
            with self.assertRaisesRegex(
                OSError, "deterministic second replacement failure"
            ):
                merged.write(output, stem="race", overwrite=True)

        self.assertFalse((output / "race.csv").exists())
        self.assertEqual(provenance.read_bytes(), old_provenance)
        self.assertEqual(receipt.read_bytes(), old_receipt)
        self.assertFalse(any(path.name.startswith(".race-") for path in output.iterdir()))

    def test_overwrite_is_documented_and_supported_as_single_writer(self):
        merged = merge_targets(
            self.dataset, self._sources(), alias_registry=self.aliases
        )
        documentation = type(merged).write.__doc__ or ""
        self.assertIn("single-writer semantics", documentation)
        self.assertIn("external lock", documentation)

        output = self.root / "sequential-overwrite"
        first = merged.write(output, stem="targets")
        expected = tuple(path.read_bytes() for path in first)
        second = merged.write(output, stem="targets", overwrite=True)
        self.assertEqual(tuple(path.read_bytes() for path in second), expected)

    def test_json_config_paths_are_relative_and_hash_bound(self):
        sources = self._sources()
        config = self.root / "targets.json"
        config.write_text(
            json.dumps(
                {
                    "sources": [
                        {
                            "path": sources[0].path.name,
                            "name": "uptake",
                            "id_column": "name",
                            "target_columns": ["uptake"],
                            "value_types": {"uptake": "float"},
                            "units": {"uptake": "mol/kg"},
                            "conditions": {
                                "uptake": {"temperature_K": 298, "pressure_bar": 1}
                            },
                        }
                    ],
                    "alias_registry": {
                        "path": self.alias_path.name,
                        "alias_columns": ["earlier_id"],
                    },
                    "feature_tables": ["rac5"],
                }
            ),
            encoding="utf-8",
        )
        merged = merge_targets_from_config(self.dataset, config)
        self.assertIn("targets/config", merged.input_hashes)
        self.assertEqual(merged.receipt()["config"]["file_name"], "targets.json")
        self.assertEqual(merged[IDS[0]].metadata["uptake"], 1.5)


if __name__ == "__main__":
    unittest.main()
