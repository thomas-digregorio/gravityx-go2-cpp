"""Tiny startup/teardown fixtures, never a full GO2 solve or evaluation."""
import concurrent.futures
import os
from pathlib import Path
import sys
import tempfile
import threading
import time
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_experiment import CompetitionTimeout, StreamingSerialEvaluation, PersistentEvaluatorProcess


class DelayedWorker:
    def __init__(self, worker_id, *args, wait_for_ready=True, **kwargs):
        self.worker_id = worker_id
        self.ready = False
        self.startup_seconds = None
        self.task_count = 0
        self.release = threading.Event()
        self.waiting = threading.Event()
        self.terminated = False
        self.fail = False
        if wait_for_ready:
            self._wait_until_ready()

    def _wait_until_ready(self):
        self.waiting.set()
        if not self.release.wait(5.0):
            raise RuntimeError("tiny fixture was not released")
        if self.terminated or self.fail:
            raise RuntimeError("startup failed")
        self.startup_seconds = 0.1
        self.ready = True

    def evaluate(self, record, deadline):
        assert self.ready
        self.task_count += 1
        return {"wall_seconds": 0.01}

    def terminate(self):
        self.terminated = True
        self.release.set()

    def stop(self, deadline=None):
        self.terminate()


class EvaluatorStartupTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        root = Path(self.temp.name)
        self.case = root / "case"
        self.output = root / "output"
        self.case.mkdir()
        self.output.mkdir()
        (self.case / "case.raw").write_text("raw\n")
        (self.case / "case.json").write_text("{}\n")
        (self.case / "case.con").write_text(
            "CONTINGENCY CTG_A\nREMOVE GENERATOR 1 FROM BUS 1\nEND\n"
            "CONTINGENCY CTG_B\nREMOVE GENERATOR 2 FROM BUS 2\nEND\nEND\n")
        for name in ("BASECASE", "CTG_A", "CTG_B"):
            (self.output / f"solution_{name}.txt").write_text("tiny\n")
        self.evaluator = root / "tiny_evaluator.py"
        self.evaluator.write_text("import sys\nsys.stdin.read()\n")

    def manager(self, deadline=None, early=False):
        m = StreamingSerialEvaluation(
            Path(sys.executable), self.evaluator, self.case, self.output,
            self.output / "internal", ["CTG_A", "CTG_B"], 2, 1,
            deadline if deadline is not None else time.perf_counter() + 10.0,
            post_screen_maximum_processes=2, completion_order_groups=True,
            persistent_evaluator_processes=True,
            prepare_persistent_pool_early=early)
        self.addCleanup(m.abort)
        return m

    @mock.patch("run_experiment.PersistentEvaluatorProcess", DelayedWorker)
    def test_delayed_ready_does_not_block_reporting_and_cannot_receive_task(self):
        m = self.manager()
        pool = concurrent.futures.ThreadPoolExecutor(max_workers=1)
        try:
            # Old synchronous construction would block until the fixture timeout.
            pool.submit(m.mark_completed, "CTG_A").result(timeout=1.0)
            worker = m.persistent_workers[0]
            self.assertTrue(worker.waiting.wait(1.0))
            m.mark_completed("CTG_B")
            self.assertEqual(m.completed_labels, {"CTG_A", "CTG_B"})
            self.assertEqual(m.running_records, {})
            self.assertEqual(worker.task_count, 0)
            worker.release.set()
            m.persistent_startup_futures[0].result(timeout=1.0)
            with m.lock:
                m._collect_persistent_startups_locked()
                self.assertEqual(m.available_persistent_workers, [worker])
                m._launch_ready_locked()
            self.assertEqual(len(m.running_records), 1)
            list(m.running_records.values())[0]["evaluation_future"].result(timeout=1.0)
            self.assertEqual(worker.task_count, 1)
        finally:
            m.abort()
            pool.shutdown(wait=True)

    @mock.patch("run_experiment.PersistentEvaluatorProcess", DelayedWorker)
    def test_startup_failure_is_not_accepted_or_hidden(self):
        m = self.manager()
        m.mark_completed("CTG_A")
        worker = m.persistent_workers[0]
        worker.fail = True
        worker.release.set()
        with self.assertRaises(RuntimeError):
            m.persistent_startup_futures[0].result(timeout=1.0)
        with m.lock, self.assertRaisesRegex(RuntimeError, "startup failed"):
            m._launch_ready_locked()
        self.assertTrue(m.aborted)
        self.assertTrue(worker.terminated)
        self.assertEqual(worker.task_count, 0)

    @mock.patch("run_experiment.PersistentEvaluatorProcess", DelayedWorker)
    def test_abort_stops_unready_process_and_releases_waiter(self):
        m = self.manager()
        m.mark_completed("CTG_A")
        worker = m.persistent_workers[0]
        self.assertTrue(worker.waiting.wait(1.0))
        future = m.persistent_startup_futures[0]
        m.abort()
        with self.assertRaises(RuntimeError):
            future.result(timeout=1.0)
        self.assertTrue(worker.terminated)
        self.assertFalse(worker.ready)
        self.assertEqual(m.running_records, {})

    @mock.patch("run_experiment.PersistentEvaluatorProcess", DelayedWorker)
    def test_promotion_does_not_import_unneeded_extra_workers(self):
        m = self.manager()
        m.promote_maximum_processes(2, remaining_screen_groups=0)
        self.assertEqual(m.maximum_processes, 2)
        self.assertEqual(len(m.persistent_workers), 1)

    @mock.patch("run_experiment.PersistentEvaluatorProcess")
    def test_expired_deadline_prevents_process_creation(self, constructor):
        m = self.manager(deadline=time.perf_counter() - 1.0)
        with self.assertRaises(CompetitionTimeout):
            m.mark_completed("CTG_A")
        constructor.assert_not_called()
        self.assertTrue(m.aborted)

    @mock.patch("run_experiment.PersistentEvaluatorProcess", DelayedWorker)
    @mock.patch("run_experiment.finalize_serial_evaluation_shard", side_effect=lambda record: record)
    def test_early_pool_preparation_does_not_increase_evaluation_concurrency(self, finalize):
        m = self.manager(early=True)
        self.assertEqual(len(m.persistent_workers), 2)
        self.assertEqual(m.maximum_processes, 1)
        self.assertFalse(m.completed_labels)
        self.assertFalse(m.running_records)
        for w in m.persistent_workers:
            w.release.set()
        for f in m.persistent_startup_futures.values():
            f.result(timeout=1.0)
        m.mark_completed("CTG_A")
        m.mark_completed("CTG_B")
        # A finished first tiny evaluation may be collected during the second
        # mark, but the configured in-flight evaluation cap remains one.
        self.assertLessEqual(len(m.running_records), 1)
        self.assertEqual(m.maximum_processes, 1)
        self.assertEqual(len(m.persistent_workers), 2)

    @mock.patch("run_experiment.PersistentEvaluatorProcess")
    def test_expired_early_pool_does_not_create_processes(self, constructor):
        with self.assertRaises(CompetitionTimeout):
            self.manager(deadline=time.perf_counter() - 1.0, early=True)
        constructor.assert_not_called()

    def test_early_pool_requires_persistent_protocol(self):
        with self.assertRaisesRegex(ValueError, "requires persistent"):
            StreamingSerialEvaluation(Path(sys.executable), self.evaluator,
                self.case, self.output, self.output / "internal", ["CTG_A"],
                1, 1, time.perf_counter() + 1,
                prepare_persistent_pool_early=True)

    def test_real_process_without_ready_is_terminated_at_deadline(self):
        m = self.manager(deadline=time.perf_counter() + 0.2)
        m.mark_completed("CTG_A")
        m.mark_completed("CTG_B")
        with self.assertRaises(CompetitionTimeout):
            m.finish()
        self.assertTrue(m.aborted)
        self.assertTrue(all(w.process.poll() is not None for w in m.persistent_workers))
        self.assertTrue(all(w.task_count == 0 for w in m.persistent_workers))

    def test_real_process_exit_before_ready_fails_closed(self):
        self.evaluator.write_text("raise SystemExit(7)\n")
        m = self.manager()
        m.mark_completed("CTG_A")
        future = m.persistent_startup_futures[0]
        with self.assertRaisesRegex(RuntimeError, "ready handshake"):
            future.result(timeout=3.0)
        with m.lock, self.assertRaisesRegex(RuntimeError, "startup failed"):
            m._launch_ready_locked()
        self.assertTrue(m.aborted)
        self.assertEqual(m.persistent_workers[0].process.returncode, 7)

    def test_idle_shutdown_is_parallel_and_uses_same_absolute_deadline(self):
        m = self.manager()
        barrier = threading.Barrier(3, timeout=2.0)
        deadlines = []
        class Worker:
            ready = True
            def terminate(self): pass
            def stop(self, deadline):
                deadlines.append(deadline)
                barrier.wait()
        m.persistent_workers = [Worker() for _ in range(3)]
        m._shutdown_evaluator_executor(True)
        self.assertEqual(deadlines, [m.deadline] * 3)
        self.assertGreaterEqual(m.evaluator_shutdown_seconds, 0.0)
        # The second cleanup is idempotent and cannot re-enter the barrier.
        m._shutdown_evaluator_executor(True)

    def test_one_shutdown_error_does_not_skip_other_worker_cleanup(self):
        m = self.manager()
        stopped = []
        class Worker:
            ready = True
            def terminate(self): pass
            def __init__(self, index): self.index = index
            def stop(self, deadline):
                stopped.append(self.index)
                if self.index == 1: raise RuntimeError("tiny shutdown failure")
        m.persistent_workers = [Worker(i) for i in range(3)]
        with self.assertRaisesRegex(RuntimeError, "tiny shutdown failure"):
            m._shutdown_evaluator_executor(True)
        self.assertEqual(sorted(stopped), [0, 1, 2])
        self.assertTrue(m.evaluator_executor_shutdown)

    def real_idle_worker(self, after_input):
        self.evaluator.write_text(
            "import sys, time\n"
            "print('GRAVITYX_EVALUATOR_READY', flush=True)\n"
            "sys.stdin.readline()\n" + after_input)
        worker = PersistentEvaluatorProcess(0, Path(sys.executable), self.evaluator,
            os.environ.copy(), self.output / "idle-worker.log")
        self.addCleanup(worker.terminate)
        return worker

    def test_real_shutdown_drains_output_and_preserves_success(self):
        worker = self.real_idle_worker("print('tiny shutdown message', flush=True)\n")
        worker.stop(time.perf_counter() + 2.0)
        self.assertEqual(worker.process.returncode, 0)
        self.assertIn("tiny shutdown message", worker.log_path.read_text())
        self.assertIn("GRAVITYX_EVALUATOR_SHUTDOWN", worker.log_path.read_text())
        worker.stop(time.perf_counter() + 2.0)

    def test_real_unresponsive_shutdown_is_bounded_and_terminated(self):
        worker = self.real_idle_worker("time.sleep(30)\n")
        with self.assertRaises(CompetitionTimeout):
            worker.stop(time.perf_counter() + 0.1)
        self.assertIsNotNone(worker.process.poll())
        self.assertTrue(worker.log_handle.closed)

    def test_real_failed_worker_exit_is_not_success(self):
        worker = self.real_idle_worker("raise SystemExit(7)\n")
        with self.assertRaisesRegex(RuntimeError, "status 7"):
            worker.stop(time.perf_counter() + 2.0)
        self.assertEqual(worker.process.returncode, 7)

    def test_exit_failure_observed_before_stop_is_still_failure(self):
        worker = self.real_idle_worker("raise SystemExit(7)\n")
        worker.process.stdin.write('{}\n')
        worker.process.stdin.flush()
        worker.process.wait(timeout=2.0)
        with self.assertRaisesRegex(RuntimeError, "status 7"):
            worker.stop(time.perf_counter() + 2.0)
        self.assertTrue(worker.log_handle.closed)


if __name__ == "__main__":
    unittest.main()
