from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest

from scripts import release_build_monitor as monitor


class ReleaseBuildMonitorTests(unittest.TestCase):
    def test_process_tree_counts_rustc_separately(self) -> None:
        table = """
          10 1 100 cargo
          11 10 200 rustc
          12 11 50 cc
          13 1 999 rustc
        """
        self.assertEqual(
            monitor.parse_process_table(table, 10),
            {
                "tree_rss_kib": 350,
                "rustc_rss_kib": 200,
                "rustc_processes": 1,
                "largest_rustc_rss_kib": 200,
            },
        )

    def test_linux_pressure_swap_and_cgroup_are_sampled(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "proc/pressure").mkdir(parents=True)
            (root / "sys/fs/cgroup").mkdir(parents=True)
            (root / "proc/meminfo").write_text(
                "MemAvailable: 1234 kB\nSwapTotal: 500 kB\nSwapFree: 400 kB\n"
            )
            (root / "proc/pressure/memory").write_text(
                "some avg10=1.25 avg60=0.50 avg300=0.10 total=12\n"
                "full avg10=0.75 avg60=0.25 avg300=0.05 total=5\n"
            )
            (root / "sys/fs/cgroup/memory.current").write_text("1000\n")
            (root / "sys/fs/cgroup/memory.max").write_text("2000\n")
            (root / "sys/fs/cgroup/memory.events").write_text("oom 2\noom_kill 1\n")
            sample = monitor.linux_memory_sample(root)
        self.assertEqual(sample["linux_memavailable_kib"], 1234)
        self.assertEqual(sample["linux_swapfree_kib"], 400)
        self.assertEqual(sample["linux_pressure_some_avg10"], "1.25")
        self.assertEqual(sample["linux_pressure_full_avg10"], "0.75")
        self.assertEqual(sample["cgroup_memory_events"], "oom:2,oom_kill:1")

    def test_runner_streams_output_and_preserves_exit_status(self) -> None:
        with tempfile.TemporaryDirectory(prefix="release monitor ") as directory:
            log = Path(directory) / "build.log"
            for code in (0, 101, 143):
                with self.subTest(code=code):
                    log.unlink(missing_ok=True)
                    completed = subprocess.run(
                        [
                            sys.executable,
                            str(Path(monitor.__file__)),
                            "--log",
                            str(log),
                            "--sample-seconds",
                            "0.01",
                            "--",
                            sys.executable,
                            "-c",
                            f"import time; print('compiler output'); time.sleep(.03); raise SystemExit({code})",
                        ],
                        capture_output=True,
                        text=True,
                        timeout=10,
                    )
                    self.assertEqual(completed.returncode, code)
                    self.assertIn("compiler output", completed.stdout)
                    self.assertIn("release_resource_sample", completed.stdout)
                    self.assertIn("release_resource_peak", completed.stdout)
                    self.assertEqual(completed.stdout, log.read_text())

    @unittest.skipUnless(hasattr(signal, "SIGTERM"), "requires POSIX-style signals")
    def test_runner_forwards_termination_and_fails(self) -> None:
        with tempfile.TemporaryDirectory(prefix="release monitor signal ") as directory:
            log = Path(directory) / "build.log"
            process = subprocess.Popen(
                [
                    sys.executable,
                    str(Path(monitor.__file__)),
                    "--log",
                    str(log),
                    "--sample-seconds",
                    "0.01",
                    "--",
                    sys.executable,
                    "-c",
                    "import time; time.sleep(10)",
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            time.sleep(0.1)
            process.terminate()
            output, _ = process.communicate(timeout=10)
        self.assertEqual(process.returncode, 128 + signal.SIGTERM)
        self.assertIn(f"signal={signal.SIGTERM}", output)
        self.assertIn("release_resource_peak", output)


if __name__ == "__main__":
    unittest.main()
