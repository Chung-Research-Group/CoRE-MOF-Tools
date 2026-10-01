"""Public dataset naming keeps release identity exact and algorithms unchanged."""
import json
from pathlib import Path
import tempfile
import unittest

from CoREMOF.dataset import CoREMOFDataset, ReleaseValidationError
from CoREMOF.splitters import split_release
from tests.test_dataset_labels import _make_release


class PublicDatasetNameTests(unittest.TestCase):
    def test_named_and_legacy_like_inputs_keep_their_exact_receipt_identity(self):
        for name in ("CoREMOF-COD", "fixture-base", "fixture-expanded"):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                _make_release(root)
                for relative in ("dataset_info.json", "parent_groups/parent_group_methods.json"):
                    path = root / relative
                    document = json.loads(path.read_text())
                    document["dataset_version"] = name
                    path.write_text(json.dumps(document))
                dataset = CoREMOFDataset.from_release(root)
                result = split_release(dataset, parent_method="none", random_state=42)
                self.assertEqual(dataset.dataset_version, name)
                self.assertEqual(result.receipt()["dataset_version"], name)
                self.assertFalse(result.receipt()["official_split"])

    def test_public_name_does_not_relax_matching_registry_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _make_release(root)
            path = root / "dataset_info.json"
            document = json.loads(path.read_text())
            document["dataset_version"] = "CoREMOF-COD"
            path.write_text(json.dumps(document))
            with self.assertRaisesRegex(ReleaseValidationError, "dataset_version does not match"):
                CoREMOFDataset.from_release(root)


if __name__ == "__main__":
    unittest.main()
