"""Tiny protocol tests: never retry a write or accept an unacknowledged task."""
import errno
import json
from pathlib import Path
import subprocess
import sys
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_experiment import CompetitionTimeout, send_worker_message


class WorkerInputTests(unittest.TestCase):
    def setUp(self):
        self.process = mock.Mock()
        self.process.poll.return_value = None
        self.lines = []
        self.options = dict(worker_name="fast-screen worker 2", deadline=120.0,
                            deadline_name="work deadline", log_path=Path("tiny.log"),
                            output_lines=self.lines)

    def evidence(self):
        self.assertEqual(len(self.lines), 1)
        self.assertTrue(self.lines[0].startswith("GRAVITYX_WORKER_INPUT_FAILURE "))
        return json.loads(self.lines[0].partition(" ")[2])

    @mock.patch("run_experiment.time.perf_counter", return_value=100.0)
    def test_success_sends_exactly_once(self, clock):
        send_worker_message(self.process, {"label": "CTG_1"}, **self.options)
        self.process.stdin.write.assert_called_once_with('{"label":"CTG_1"}\n')
        self.process.stdin.flush.assert_called_once_with()
        self.assertEqual(self.lines, [])

    @mock.patch("run_experiment.time.perf_counter", return_value=121.0)
    def test_no_task_is_sent_after_deadline(self, clock):
        with self.assertRaises(CompetitionTimeout):
            send_worker_message(self.process, {"label": "late"}, **self.options)
        self.process.stdin.write.assert_not_called()
        event = self.evidence()
        self.assertIsNone(event["returncode"])
        self.assertLess(event["seconds_remaining"], 0)

    @mock.patch("run_experiment.time.perf_counter", return_value=100.0)
    def test_known_timeout_exit_before_write(self, clock):
        self.process.poll.return_value = 124
        with self.assertRaises(CompetitionTimeout):
            send_worker_message(self.process, {"label": "unstarted"}, **self.options)
        self.process.stdin.write.assert_not_called()
        self.assertEqual(self.evidence()["returncode"], 124)

    @mock.patch("run_experiment.time.perf_counter", return_value=100.0)
    def test_windows_einval_collects_timeout_exit_without_retry(self, clock):
        original = OSError(errno.EINVAL, "Invalid argument")
        self.process.stdin.write.side_effect = original
        self.process.poll.side_effect = [None, None]
        self.process.wait.return_value = 124
        with self.assertRaises(CompetitionTimeout) as caught:
            send_worker_message(self.process, {"label": "CTG_last"}, **self.options)
        self.assertIs(caught.exception.__cause__, original)
        self.process.stdin.write.assert_called_once()
        self.process.stdin.flush.assert_not_called()
        self.process.wait.assert_called_once_with(timeout=0.05)
        event = self.evidence()
        self.assertEqual(event["errno"], errno.EINVAL)
        self.assertEqual(event["returncode"], 124)
        self.assertEqual(event["label"], "CTG_last")

    @mock.patch("run_experiment.time.perf_counter", return_value=100.0)
    def test_flush_failure_does_not_resend_or_hide_crash(self, clock):
        self.process.stdin.flush.side_effect = BrokenPipeError(errno.EPIPE, "closed")
        self.process.poll.side_effect = [None, 9]
        with self.assertRaises(RuntimeError) as caught:
            send_worker_message(self.process, {"label": "crashed"}, **self.options)
        self.assertNotIsInstance(caught.exception, CompetitionTimeout)
        self.process.stdin.write.assert_called_once()
        self.process.stdin.flush.assert_called_once()
        self.process.wait.assert_not_called()
        self.assertEqual(self.evidence()["returncode"], 9)

    @mock.patch("run_experiment.time.perf_counter", return_value=100.0)
    def test_unresolved_early_failure_remains_unknown_not_timeout(self, clock):
        self.process.stdin.write.side_effect = OSError(errno.EINVAL, "Invalid argument")
        self.process.wait.side_effect = subprocess.TimeoutExpired("tiny", 0.05)
        with self.assertRaises(RuntimeError) as caught:
            send_worker_message(self.process, {"label": "unknown"}, **self.options)
        self.assertNotIsInstance(caught.exception, CompetitionTimeout)
        self.process.stdin.write.assert_called_once()
        self.process.wait.assert_called_once_with(timeout=0.05)
        self.assertIsNone(self.evidence()["returncode"])

    @mock.patch("run_experiment.time.perf_counter", side_effect=[119.999, 120.001, 120.002])
    def test_deadline_race_reports_elapsed_limit_without_inventing_exit(self, clock):
        self.process.stdin.write.side_effect = OSError(errno.EINVAL, "Invalid argument")
        with self.assertRaises(CompetitionTimeout):
            send_worker_message(self.process, {"label": "race"}, **self.options)
        self.process.stdin.write.assert_called_once()
        self.process.wait.assert_not_called()
        event = self.evidence()
        self.assertIsNone(event["returncode"])
        self.assertEqual(event["errno"], errno.EINVAL)
        self.assertLess(event["seconds_remaining"], 0)

    @mock.patch("run_experiment.time.perf_counter", return_value=121.0)
    def test_shutdown_can_use_finalization_reserve_but_sends_no_task(self, clock):
        send_worker_message(self.process, {"stop": True}, enforce_deadline=False,
                            **self.options)
        self.process.stdin.write.assert_called_once_with('{"stop":true}\n')
        self.process.stdin.flush.assert_called_once_with()
        self.assertEqual(self.lines, [])

    def test_task_cannot_opt_out_of_deadline(self):
        with self.assertRaisesRegex(ValueError, "only a stop message"):
            send_worker_message(self.process, {"label": "forbidden"},
                                enforce_deadline=False, **self.options)
        self.process.stdin.write.assert_not_called()


if __name__ == "__main__":
    unittest.main()
