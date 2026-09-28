"""No-network regressions for model download failure and cache publication."""

import hashlib
import os
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

from CoREMOF._prediction_download import download_model_file


class Response:
    def __init__(self, chunks=(b"model",), failure=None):
        self.chunks = chunks
        self.failure = failure
        self.closed = False

    def raise_for_status(self):
        pass

    def iter_content(self, chunk_size):
        yield from self.chunks
        if self.failure:
            raise self.failure

    def close(self):
        self.closed = True


class PredictionDownloadsTests(unittest.TestCase):
    def request(self, response):
        module = types.ModuleType("requests")
        module.get = lambda *args, **kwargs: response
        return patch.dict(sys.modules, {"requests": module})

    def test_basename_download_and_complete_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            old = Path.cwd()
            try:
                os.chdir(directory)
                response = Response((b"mo", b"", b"del"))
                with self.request(response):
                    self.assertTrue(download_model_file("fixture", "model", expected_sha256=hashlib.sha256(b"model").hexdigest()))
                self.assertEqual(Path("model").read_bytes(), b"model")
                self.assertEqual(list(Path('.').iterdir()), [Path("model")])
                self.assertTrue(response.closed)
            finally:
                os.chdir(old)

    def test_failures_never_publish_partial_data(self):
        for response, digest in ((Response((b"partial",), RuntimeError("interrupted")), None),
                                 (Response(()), None), (Response(), "0" * 64)):
            with self.subTest(digest=digest), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                with self.request(response), self.assertRaises((ValueError, RuntimeError)):
                    download_model_file("fixture", root / "model", expected_sha256=digest)
                self.assertEqual(list(root.iterdir()), [])
                self.assertTrue(response.closed)

    def test_existing_cache_checked_without_requests(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model"
            path.write_bytes(b"existing")
            with patch.dict(sys.modules, {"requests": None}):
                self.assertFalse(download_model_file("fixture", path))
                self.assertFalse(download_model_file("fixture", path, expected_sha256=hashlib.sha256(b"existing").hexdigest().upper()))
                with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                    download_model_file("fixture", path, expected_sha256="0" * 64)
            self.assertEqual(path.read_bytes(), b"existing")

    def test_empty_cache_and_symlink_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "empty"
            path.touch()
            with self.assertRaises(ValueError):
                download_model_file("fixture", path)
            link = root / "linked"
            link.symlink_to(root / "missing")
            with self.assertRaises(ValueError):
                download_model_file("fixture", link)
            folder = root / "folder"
            folder.symlink_to(root, target_is_directory=True)
            with self.assertRaises(ValueError):
                download_model_file("fixture", folder / "model")

    def test_concurrent_winner_is_not_replaced(self):
        for correct in (True, False):
            with self.subTest(correct=correct), tempfile.TemporaryDirectory() as directory:
                destination = Path(directory) / "model"
                response = Response()

                def raced(source, target):
                    destination.write_bytes(b"model" if correct else b"other")
                    raise FileExistsError(target)

                with self.request(response), patch("CoREMOF._prediction_download.os.link", side_effect=raced):
                    if correct:
                        self.assertFalse(download_model_file("fixture", destination, expected_sha256=hashlib.sha256(b"model").hexdigest()))
                    else:
                        with self.assertRaises(ValueError):
                            download_model_file("fixture", destination, expected_sha256=hashlib.sha256(b"model").hexdigest())
                self.assertEqual(destination.read_bytes(), b"model" if correct else b"other")
                self.assertEqual(list(Path(directory).iterdir()), [destination])

    def test_malformed_checksum_rejected_before_network(self):
        for digest in ("abc", "g" * 64, 123):
            with self.subTest(digest=digest), self.assertRaises(ValueError):
                download_model_file("fixture", "unused", expected_sha256=digest)


if __name__ == "__main__":
    unittest.main()
