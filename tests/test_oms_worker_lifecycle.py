"""Exercise OMS process supervision with disposable workers, without chemistry."""

import ast
from contextlib import contextmanager, redirect_stdout
import io
import multiprocessing
import os
from pathlib import Path
import signal
import tempfile
import time
import unittest
from unittest.mock import Mock, patch


SOURCE = Path(__file__).resolve().parents[1] / "CoREMOF/calculation/mof_collection.py"
tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
collection = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == "MofCollection")
definitions = [node for node in collection.body if getattr(node, "name", None) in {"analyse_mofs", "_run_batch"}]
NAMESPACE = {"time": time}
exec(compile(ast.Module(body=definitions, type_ignores=[]), str(SOURCE), "exec"), NAMESPACE)


class StubCollection:
    separator = "OMS test"
    analyse_mofs = NAMESPACE["analyse_mofs"]

    def __init__(self, batches):
        self.batches = batches
        self.mof_coll = []
        self._make_batches = Mock()
        self._validate_properties = Mock()

    def _run_batch(self, index, batch, overwrite, status):
        if batch and batch[0]["action"] == "exit_in_status_update":
            if hasattr(status, "get_lock"):
                status.get_lock().acquire()
            os._exit(7)
        if batch and batch[0]["action"] == "false_complete":
            status[index] = -1
            os._exit(9)
        NAMESPACE["_run_batch"](self, index, batch, overwrite, status)

    def _analyse(self, item, overwrite):
        if item["action"] == "exit":
            os._exit(7)
        if item["action"] == "clean_exit":
            os._exit(0)
        if item["action"] == "write":
            # This worker must be allowed to finish after its sibling fails.
            time.sleep(0.05)
            Path(item["path"]).write_text("completed stub result", encoding="utf-8")


@contextmanager
def watchdog(seconds=3):
    def expired(signum, frame):
        raise AssertionError("OMS supervisor did not notice terminated workers")

    previous = signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


@unittest.skipUnless("fork" in multiprocessing.get_all_start_methods(), "requires disposable POSIX fork workers")
class OmsWorkerLifecycleTests(unittest.TestCase):
    def setUp(self):
        self.context = multiprocessing.get_context("fork")
        self.workers = []

        def process(**kwargs):
            worker = self.context.Process(**kwargs)
            worker.join = Mock(wraps=worker.join)
            self.workers.append(worker)
            return worker

        self.environment = patch.dict(NAMESPACE, Process=process, Array=self.context.Array)
        self.environment.start()
        self.addCleanup(self.environment.stop)
        self.addCleanup(self.cleanup_workers)

    def cleanup_workers(self):
        # These are only the exact disposable test workers created above.
        for worker in self.workers:
            if worker.is_alive():
                worker.terminate()
            worker.join(timeout=2)

    def run_collection(self, batches):
        collection = StubCollection(batches)
        with watchdog(), redirect_stdout(io.StringIO()):
            collection.analyse_mofs(num_batches=len(batches))
        return collection

    def assert_joined(self):
        for worker in self.workers:
            self.assertFalse(worker.is_alive())
            worker.join.assert_called()

    def test_abrupt_nonzero_exit_is_reported_instead_of_polling_forever(self):
        collection = StubCollection([[{"mof_name": "crash", "action": "exit"}]])
        with watchdog(), redirect_stdout(io.StringIO()), self.assertRaisesRegex(RuntimeError, "batch 1 exited with code 7"):
            collection.analyse_mofs()
        self.assertEqual(self.workers[0].exitcode, 7)
        collection._validate_properties.assert_not_called()
        self.assert_joined()

    def test_zero_exit_without_completion_marker_is_unavailable(self):
        with self.assertRaisesRegex(RuntimeError, "code 0 without a completion status"):
            self.run_collection([[{"mof_name": "incomplete", "action": "clean_exit"}]])
        self.assert_joined()

    def test_exit_during_status_update_cannot_poison_a_shared_lock(self):
        with self.assertRaisesRegex(RuntimeError, "batch 1 exited with code 7"):
            self.run_collection([[{"mof_name": "status update crash", "action": "exit_in_status_update"}]])
        self.assert_joined()

    def test_completion_marker_does_not_hide_nonzero_process_exit(self):
        with self.assertRaisesRegex(RuntimeError, "batch 1 exited with code 9"):
            self.run_collection([[{"mof_name": "false marker", "action": "false_complete"}]])
        self.assert_joined()

    def test_failed_worker_does_not_cancel_healthy_sibling_or_delete_results(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            existing = root / "existing.txt"
            existing.write_bytes(b"retained evidence")
            output = root / "healthy.txt"
            batches = [
                [{"mof_name": "crash", "action": "exit"}],
                [{"mof_name": "healthy", "action": "write", "path": str(output)}],
            ]
            with self.assertRaisesRegex(RuntimeError, "Existing results were preserved"):
                self.run_collection(batches)
            self.assertEqual(existing.read_bytes(), b"retained evidence")
            self.assertEqual(output.read_text(), "completed stub result")
            self.assertEqual([worker.exitcode for worker in self.workers], [7, 0])
            self.assert_joined()

    def test_successful_and_empty_batches_are_joined_before_final_validation(self):
        collection = self.run_collection([[{"mof_name": "healthy", "action": "success"}], []])
        self.assertEqual([worker.exitcode for worker in self.workers], [0, 0])
        collection._validate_properties.assert_called_once_with(["has_oms"])
        self.assert_joined()


if __name__ == "__main__":
    unittest.main()
