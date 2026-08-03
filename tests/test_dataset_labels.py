import csv
import hashlib
import json
from pathlib import Path
import sys
import subprocess
import tempfile
from types import ModuleType
import unittest
from unittest.mock import patch

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
        self.assertEqual(classified.ids_for_label("ambiguous"), ("ASR-SI-2025-0001",))
        self.assertEqual(classified.unchecked_ids, ("ION-CSD-2024-0001",))
        with self.assertRaises(ValueError):
            classified.ids_for_label("MAYBE")

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


if __name__ == "__main__":
    unittest.main()
