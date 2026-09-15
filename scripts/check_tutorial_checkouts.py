#!/usr/bin/env python3
"""Replay tutorial checks under LF and CRLF Git checkout policies on any host.

Uses local, temporary clones sharing read-only objects; never builds Rust,
fetches references, initializes submodules, or modifies the source checkout.
The binary must already be built from the revision being checked. A local
attributes-only overlay is explicit and is not exact-commit acceptance.
"""

import argparse
import math
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
MODES = (("lf", "false", "lf"), ("crlf", "true", "crlf"))


def git(root, *args):
    return subprocess.run(
        ["git", "-C", str(root), "-c", f"core.hooksPath={Path(root).resolve() / '.disabled-hooks'}",
         "-c", "submodule.recurse=false", *args],
        check=True, capture_output=True,
    ).stdout


def prepare_checkout(source, destination, mode, attributes=None, revision=None):
    """Preserve Git history for review badges, without copying ignored caches."""
    name, autocrlf, eol = mode
    revision = revision or git(source, "rev-parse", "HEAD").decode().strip()
    git(source, "clone", "--shared", "--no-checkout", "--quiet",
        str(source.resolve()), str(destination))
    git(destination, "config", "core.autocrlf", autocrlf)
    git(destination, "config", "core.eol", eol)
    git(destination, "config", "core.longpaths", "true")
    git(destination, "update-ref", "--no-deref", "HEAD", revision)
    git(destination, "read-tree", revision)
    if attributes is not None:
        (destination / ".gitattributes").write_bytes(attributes)
        git(destination, "add", "--", ".gitattributes")
    # Export once, after the policy is in the index; otherwise Git may retain
    # already-written CRLF files whose stat cache still says they are unchanged.
    git(destination, "checkout-index", "--all", "--force")
    print(f"Tutorial checkout: {name}, revision={revision}, "
          f"attributes={'overlay' if attributes is not None else 'committed'}", flush=True)
    return revision


def check_checkout(binary, checkout, mode, timeout):
    env = os.environ.copy()
    # This gate checks committed offline teaching artifacts, not live services.
    env.pop("GENTLE_TEST_ONLINE", None)
    for argument in ("--check", "tutorial-check"):
        print(f"[{mode}] {binary.name} {argument}", flush=True)
        try:
            subprocess.run([str(binary), argument], cwd=checkout, env=env,
                           check=True, timeout=timeout)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
            raise RuntimeError(f"{mode} checkout failed: {argument}: {error}") from error


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--timeout-seconds", type=float, default=900)
    parser.add_argument("--attributes-from", type=Path,
                        help="Explicit local .gitattributes overlay; not exact-commit acceptance")
    args = parser.parse_args(argv)
    binary = args.binary.resolve()
    if not binary.is_file():
        parser.error(f"Build gentle_examples_docs first: {binary} does not exist")
    if not math.isfinite(args.timeout_seconds) or args.timeout_seconds <= 0:
        parser.error("--timeout-seconds must be positive")
    attributes = args.attributes_from.read_bytes() if args.attributes_from else None
    revision = git(args.repo_root, "rev-parse", "HEAD").decode().strip()
    for mode in MODES:
        # Release each checkout before creating the next one on small machines.
        with tempfile.TemporaryDirectory(prefix=f"gentle-tutorial-{mode[0]}-",
                                         dir=os.environ.get("RUNNER_TEMP")) as tmp:
            checkout = Path(tmp) / "repo"
            prepare_checkout(args.repo_root, checkout, mode, attributes, revision)
            check_checkout(binary, checkout, mode[0], args.timeout_seconds)
    print("Tutorial checkout portability: LF and CRLF passed", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as error:
        detail = error.stderr.decode(errors="replace") if error.stderr else ""
        raise SystemExit(f"{error}\n{detail}") from error
    except (OSError, RuntimeError) as error:
        raise SystemExit(str(error)) from error
