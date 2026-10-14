"""Synthetic subprocess checks for the Cargo default-stack test runner.

Inline Python children and mocked exit statuses reproduce the environment and
argv boundary on every platform. These tests never invoke Cargo or Rust.
"""

from contextlib import redirect_stdout
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tomllib
import unittest
from unittest.mock import patch

from scripts import default_stack_test_runner as runner


class DefaultStackRunnerTests(unittest.TestCase):
    def test_config_quotes_native_interpreter_and_wrapper_paths(self):
        interpreter = r"C:\Synthetic Tools\Python\python.exe"
        wrapper = Path(r"C:\Synthetic Project\default_stack_test_runner.py")
        with patch.object(runner.sys, "executable", interpreter), \
                patch.object(runner.Path, "resolve", return_value=wrapper):
            config = runner.cargo_runner_config()
        self.assertEqual(tomllib.loads(config)["target"]["cfg(all())"]["runner"],
                         [interpreter, str(wrapper)])

    def test_real_child_cannot_inherit_cargos_stack_override(self):
        command = [sys.executable, "-c",
                   "import os; assert 'RUST_MIN_STACK' not in os.environ"]
        output = io.StringIO()
        with patch.dict(os.environ, {"RUST_MIN_STACK": "16777216"}), redirect_stdout(output):
            self.assertEqual(runner.run_test(command), 0)
            self.assertEqual(os.environ["RUST_MIN_STACK"], "16777216")
        receipt = json.loads(output.getvalue())
        self.assertEqual(receipt["inherited_rust_min_stack"], "16777216")
        self.assertIsNone(receipt["child_rust_min_stack"])

    def test_runner_preserves_arguments_and_cargo_runtime_environment(self):
        command = [r"C:\Synthetic Project\gentle-tests.exe", "dbsnp_", "--test-threads=1"]
        with patch.dict(os.environ, {"RUST_MIN_STACK": "16777216",
                                     "LD_LIBRARY_PATH": "/synthetic/target/debug/deps"}), \
                patch.object(runner.subprocess, "run", return_value=subprocess.CompletedProcess(command, 7)) as run, \
                redirect_stdout(io.StringIO()):
            self.assertEqual(runner.run_test(command), 7)
        self.assertEqual(run.call_args.args, (command,))
        environment = run.call_args.kwargs["env"]
        self.assertNotIn("RUST_MIN_STACK", environment)
        self.assertEqual(environment["LD_LIBRARY_PATH"], "/synthetic/target/debug/deps")
        self.assertIs(run.call_args.kwargs["check"], False)
        self.assertNotIn("shell", run.call_args.kwargs)

    def test_signal_failure_remains_nonzero(self):
        with patch.object(runner.subprocess, "run", return_value=subprocess.CompletedProcess(["synthetic"], -6)), \
                redirect_stdout(io.StringIO()):
            self.assertEqual(runner.run_test(["synthetic"]), 134)

    def test_empty_command_is_refused_before_spawn(self):
        with patch.object(runner.subprocess, "run") as run:
            with self.assertRaises(ValueError):
                runner.run_test([])
            run.assert_not_called()

    def test_config_mode_only_prints_configuration(self):
        output = io.StringIO()
        with patch.object(runner.subprocess, "run") as run, redirect_stdout(output):
            self.assertEqual(runner.main(["--cargo-config"]), 0)
            run.assert_not_called()
        self.assertIn("cfg(all())", tomllib.loads(output.getvalue())["target"])


if __name__ == "__main__":
    unittest.main()
