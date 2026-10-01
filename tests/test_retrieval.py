"""Bounded release retrieval integration tests; no external service is used."""

import copy
from functools import partial
import hashlib
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch
from urllib.error import HTTPError

from CoREMOF.dataset import CoREMOFDataset, ReleaseValidationError
from CoREMOF.retrieval import CATALOG_SCHEMA, RECEIPT_PATH, RetrievalError, fetch_release
from test_dataset_labels import _make_release, _metadata_rows, _write_csv


class ReleaseRetrievalTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.source = self.root / "served" / "data" / "vtest"
        _make_release(self.source)
        manifest = []
        for row in _metadata_rows():
            data = row["structure_id"].encode("utf-8")
            cif = self.source / row["cif_file"]
            cif.parent.mkdir(parents=True, exist_ok=True)
            cif.write_bytes(data)
            manifest.append({"structure_id": row["structure_id"], "cif_file": row["cif_file"],
                             "size_bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()})
        _write_csv(self.source / "manifests" / "cif_manifest.csv", tuple(manifest[0]), manifest)
        files = []
        for path in sorted(self.source.rglob("*")):
            if path.is_file():
                logical = path.relative_to(self.source).as_posix()
                data = path.read_bytes()
                files.append({"path": logical, "url": "data/vtest/" + logical,
                              "size_bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()})
        self.catalog = {"schema_version": CATALOG_SCHEMA, "releases": {"vtest": {"files": files}}}
        self.catalog_path = self.root / "served" / "catalog.json"
        self.destination = self.root / "retrieved"
        self.pin = self.write_catalog()

    def write_catalog(self, catalog=None):
        data = json.dumps(self.catalog if catalog is None else catalog, sort_keys=True).encode("utf-8")
        self.catalog_path.write_bytes(data)
        return hashlib.sha256(data).hexdigest()

    def fetch(self, **kwargs):
        return fetch_release(self.catalog_path, "vtest", self.destination,
                             catalog_sha256=self.pin, **kwargs)

    def assert_unpublished(self):
        self.assertFalse(self.destination.exists())
        self.assertEqual(list(self.root.glob(".coremof-retrieval-*")), [])

    def start_server(self):
        requests = []

        class Handler(SimpleHTTPRequestHandler):
            def do_GET(self):
                requests.append(self.path)
                super().do_GET()

            def log_message(self, *args):
                pass

        handler = partial(Handler, directory=str(self.catalog_path.parent))
        server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
        worker = threading.Thread(target=server.serve_forever, daemon=True)
        worker.start()

        def finish():
            server.shutdown()
            worker.join(timeout=5)
            server.server_close()

        self.addCleanup(finish)
        return "http://127.0.0.1:{}/catalog.json".format(server.server_port), requests

    def test_local_catalog_strict_cifs_and_receipt_preserve_exact_bytes(self):
        result = self.fetch(verify_cif_files=True)
        self.assertEqual(result, self.destination)
        dataset = CoREMOFDataset.from_release(result, verify_cif_files=True)
        self.assertEqual((dataset.dataset_version, len(dataset)), ("vtest", 4))
        receipt = json.loads((result / RECEIPT_PATH).read_text())
        self.assertEqual(receipt["catalog_sha256"], self.pin)
        self.assertEqual(receipt["release_input_hashes"], dict(dataset.input_hashes))
        self.assertTrue(receipt["cif_files_verified"])
        self.assertEqual(receipt["downloaded_bytes"], sum(f["size_bytes"] for f in self.catalog["releases"]["vtest"]["files"]))
        self.assertEqual(receipt["file_count"], len(receipt["files"]))
        self.assertEqual(len(receipt["retrieval_source_sha256"]), 64)
        self.assertNotIn(str(self.root), json.dumps(receipt))
        self.assertNotIn("url", json.dumps(receipt))
        for item in receipt["files"]:
            self.assertEqual((result / item["path"]).read_bytes(), (self.source / item["path"]).read_bytes())

    def test_http_catalog_relative_resources_and_exact_version_selection(self):
        # An unrelated version is not requested or silently selected.
        self.catalog["releases"]["vother"] = {"files": [{"path": "unused", "url": "missing",
                                                          "sha256": "0" * 64, "size_bytes": 1}]}
        self.pin = self.write_catalog()
        url, requests = self.start_server()
        result = fetch_release(url + "?private_token=secret", "vtest", self.destination,
                               catalog_sha256=self.pin, verify_cif_files=True)
        self.assertEqual(result, self.destination)
        self.assertEqual(len(requests), 1 + len(self.catalog["releases"]["vtest"]["files"]))
        self.assertNotIn("/missing", requests)
        receipt = (result / RECEIPT_PATH).read_text()
        self.assertNotIn("secret", receipt)
        self.assertNotIn("127.0.0.1", receipt)

    def test_wrong_catalog_pin_stops_before_payload_fetch(self):
        with patch("CoREMOF.retrieval._download_file") as download:
            self.pin = "0" * 64
            with self.assertRaisesRegex(RetrievalError, "catalog checksum mismatch"):
                self.fetch()
            download.assert_not_called()
        self.assert_unpublished()

    def test_corrupt_truncated_or_oversized_resource_is_not_published(self):
        original = copy.deepcopy(self.catalog)
        for mutation in ("checksum", "larger_declared_size", "smaller_declared_size"):
            with self.subTest(mutation=mutation):
                self.catalog = copy.deepcopy(original)
                item = self.catalog["releases"]["vtest"]["files"][0]
                if mutation == "checksum":
                    item["sha256"] = "0" * 64
                else:
                    item["size_bytes"] += 1 if mutation == "larger_declared_size" else -1
                self.pin = self.write_catalog()
                with self.assertRaises(RetrievalError):
                    self.fetch()
                self.assert_unpublished()

    def test_unsafe_duplicate_reserved_and_file_directory_paths_are_rejected(self):
        original = copy.deepcopy(self.catalog)
        for unsafe in ("../outside", "/outside", "a/../outside", "a//b", "a\\b", "./a",
                       "manifests", RECEIPT_PATH, "DATASET_INFO.JSON", "dataset_info.json/child"):
            with self.subTest(path=unsafe):
                self.catalog = copy.deepcopy(original)
                item = copy.deepcopy(self.catalog["releases"]["vtest"]["files"][0])
                item["path"] = unsafe
                self.catalog["releases"]["vtest"]["files"].append(item)
                self.pin = self.write_catalog()
                with patch("CoREMOF.retrieval._download_file") as download:
                    with self.assertRaises(RetrievalError):
                        self.fetch()
                    download.assert_not_called()
                self.assert_unpublished()
        self.assertFalse((self.root / "outside").exists())

    def test_absent_version_and_mislabeled_dataset_are_rejected(self):
        with self.assertRaisesRegex(RetrievalError, "absent"):
            fetch_release(self.catalog_path, "latest", self.destination, catalog_sha256=self.pin)
        self.catalog["releases"]["vwrong"] = self.catalog["releases"].pop("vtest")
        self.pin = self.write_catalog()
        with self.assertRaisesRegex(RetrievalError, "dataset version differs"):
            fetch_release(self.catalog_path, "vwrong", self.destination, catalog_sha256=self.pin)
        self.assert_unpublished()

    def test_existing_destination_and_publication_failure_preserve_state(self):
        self.destination.mkdir()
        sentinel = self.destination / "user.txt"
        sentinel.write_text("preserve me")
        with self.assertRaises(FileExistsError):
            self.fetch()
        self.assertEqual(sentinel.read_text(), "preserve me")
        # Choose a distinct absent output; do not delete the protected directory.
        self.destination = self.root / "failed-publication"
        with patch("CoREMOF.retrieval.publish_directory", side_effect=OSError("synthetic publish failure")):
            with self.assertRaisesRegex(OSError, "synthetic publish failure"):
                self.fetch()
        self.assert_unpublished()
        self.assertEqual(sentinel.read_text(), "preserve me")

    def test_metadata_only_retrieval_and_opt_in_cif_requirement(self):
        self.catalog["releases"]["vtest"]["files"] = [
            item for item in self.catalog["releases"]["vtest"]["files"] if not item["path"].startswith("cifs/")]
        self.pin = self.write_catalog()
        with self.assertRaises((ReleaseValidationError, FileNotFoundError)):
            self.fetch(verify_cif_files=True)
        self.assert_unpublished()
        self.fetch()
        receipt = json.loads((self.destination / RECEIPT_PATH).read_text())
        self.assertFalse(receipt["cif_files_verified"])
        self.assertFalse((self.destination / "cifs").exists())

    def test_total_byte_limit_rejects_before_payload_fetch(self):
        total = sum(item["size_bytes"] for item in self.catalog["releases"]["vtest"]["files"])
        with patch("CoREMOF.retrieval._download_file") as download:
            with self.assertRaisesRegex(RetrievalError, "max_total_bytes"):
                self.fetch(max_total_bytes=total - 1)
            download.assert_not_called()
        self.assert_unpublished()
        self.fetch(max_total_bytes=total)

    def test_duplicate_json_keys_and_nonfinite_sizes_are_rejected(self):
        valid = self.catalog_path.read_bytes()
        for data in (valid.replace(b'"releases":', b'"schema_version":"duplicate", "releases":', 1),
                     valid.replace(b'"size_bytes": ', b'"size_bytes": NaN, "ignored": ', 1)):
            with self.subTest(data=data[:50]):
                self.catalog_path.write_bytes(data)
                self.pin = hashlib.sha256(data).hexdigest()
                with self.assertRaises(RetrievalError):
                    self.fetch()
                self.assert_unpublished()

    def test_missing_http_resource_omits_private_url_from_error(self):
        self.catalog["releases"]["vtest"]["files"][0]["url"] = "missing?private_token=secret"
        self.pin = self.write_catalog()
        url, _ = self.start_server()
        with self.assertRaisesRegex(RetrievalError, "HTTP 404") as caught:
            fetch_release(url, "vtest", self.destination, catalog_sha256=self.pin)
        self.assertNotIn("secret", str(caught.exception))
        self.assert_unpublished()
        body = io.BytesIO(b"private response body")
        failure = HTTPError(url + "?private_token=secret", 403, "forbidden", {}, body)
        with patch("CoREMOF.retrieval.urlopen", side_effect=failure):
            with self.assertRaisesRegex(RetrievalError, "HTTP 403"):
                self.fetch()
        self.assertTrue(body.closed)

    def test_executable_example_uses_same_verified_api(self):
        script = Path(__file__).resolve().parents[1] / "examples" / "fetch_release.py"
        result = subprocess.run([sys.executable, "-S", str(script), str(self.catalog_path),
                                 "vtest", str(self.destination), "--catalog-sha256", self.pin,
                                 "--verify-cifs"], capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), str(self.destination))
        self.assertEqual(len(CoREMOFDataset.from_release(self.destination, verify_cif_files=True)), 4)


if __name__ == "__main__":
    unittest.main()
