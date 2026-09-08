#!/usr/bin/env python3
"""Bind packaging checks to one candidate; never create tags or publish assets.

Uses only the Python standard library. The workflow owns publication authority;
this helper independently checks inputs, checkout identity and build receipts.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tomllib


def full_sha(value: str) -> str:
    if not re.fullmatch(r"[0-9a-f]{40}", value):
        raise ValueError("candidate_sha must be a full, lowercase 40-character commit SHA")
    return value


def release_tag(value: str) -> str:
    if not re.fullmatch(r"v[0-9][0-9A-Za-z.+-]*", value):
        raise ValueError("tag must be a version label such as v0.1.0-internal.10")
    return value


def publication_allowed(event: str, action: str, requested: str) -> bool:
    if requested not in ("true", "false"):
        raise ValueError("publish must be explicitly true or false")
    if event == "workflow_dispatch":
        return requested == "true"
    if event == "release" and action == "published":
        return True
    if event == "push" and requested == "false":
        return False
    raise ValueError("Unsupported release event or implicit publication request")


def git(root: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def tag_revision(root: Path, tag: str) -> str:
    # FETCH_HEAD observes the remote tag without moving/creating any local tag.
    subprocess.run(
        ["git", "-C", str(root), "fetch", "--no-tags", "origin", f"refs/tags/{release_tag(tag)}"],
        check=True,
    )
    return git(root, "rev-parse", "FETCH_HEAD^{commit}")


def validate_checkout(
    root: Path, tag: str, expected_revision: str, expected_lock: str = ""
) -> dict:
    tag = release_tag(tag)
    revision = full_sha(git(root, "rev-parse", "HEAD"))
    if revision != full_sha(expected_revision):
        raise ValueError("Checkout does not match the expected candidate SHA")
    manifest = tomllib.loads((root / "Cargo.toml").read_text())
    version = manifest["workspace"]["package"]["version"]
    if tag != f"v{version}":
        raise ValueError(f"Tag label {tag} does not match package version {version}")
    lock_digest = hashlib.sha256((root / "Cargo.lock").read_bytes()).hexdigest()
    if expected_lock and lock_digest != expected_lock:
        raise ValueError("Cargo.lock no longer matches the candidate receipt")
    if git(root, "status", "--porcelain", "--untracked-files=no", "--", "Cargo.toml", "Cargo.lock"):
        raise ValueError("Candidate Cargo metadata has uncommitted changes")
    return {"tag": tag, "revision": revision, "cargo_lock_sha256": lock_digest}


def prepare(root: Path, env: dict[str, str]) -> dict:
    event = env["CANDIDATE_EVENT"]
    publish = publication_allowed(event, env.get("CANDIDATE_ACTION", ""), env["PUBLISH_REQUESTED"])
    tag = release_tag(env["RELEASE_TAG"])
    expected = env.get("CANDIDATE_SHA", "")
    if event == "workflow_dispatch":
        full_sha(expected)  # No branch, tag or empty-input fallback on manual runs.
    else:
        expected = full_sha(env["CANDIDATE_EVENT_SHA"])
        if event == "push" and env.get("CANDIDATE_REF") != f"refs/tags/{tag}":
            raise ValueError("Only tag pushes are supported")
    record = validate_checkout(root, tag, expected)
    if publish and tag_revision(root, tag) != expected:
        raise ValueError("Release tag does not point to the approved candidate SHA")
    return {
        "schema": "gentle.release_candidate.v1",
        **record,
        "workflow_revision": full_sha(env["WORKFLOW_REVISION"]),
        "publish": publish,
        "mode": "publish" if publish else "validate_only",
    }


def collect_installers(root: Path, candidate: dict) -> dict:
    """Reject missing/mixed candidate receipts before retaining or publishing."""
    release_tag(candidate["tag"])
    full_sha(candidate["revision"])
    receipts = [json.loads(path.read_text()) for path in root.rglob("*.build.json")]
    platforms = {"linux": "tar.gz", "macos": "dmg", "windows": "zip"}
    if len(receipts) != 3 or {r.get("platform") for r in receipts} != set(platforms):
        raise ValueError("Expected exactly Linux, macOS and Windows build receipts")
    for receipt in receipts:
        if receipt.get("schema") != "gentle.release_build.v1":
            raise ValueError("Unsupported build receipt schema")
        for key in ("tag", "revision", "cargo_lock_sha256", "workflow_revision", "mode"):
            if receipt.get(key) != candidate[key]:
                raise ValueError(f"Build receipt {key} does not match the selected candidate")
        if receipt.get("profile") != "release" or receipt.get("features") != ["script-interfaces"]:
            raise ValueError("Build receipt does not describe the release-shaped script bundle")
        if not all(isinstance(receipt.get(key), str) and receipt[key] for key in ("rustc", "cargo")):
            raise ValueError("Build receipt is missing Rust/Cargo toolchain identity")
    artifacts = []
    for receipt in receipts:
        extension = platforms[receipt["platform"]]
        name = f"gentle-{candidate['tag']}-{receipt['platform']}-{receipt['arch']}.{extension}"
        paths = list(root.rglob(f"*.{extension}"))
        if len(paths) != 1 or paths[0].name != name or paths[0].stat().st_size == 0:
            raise ValueError(f"Expected one non-empty candidate artifact named {name}")
        path = paths[0]
        with path.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        if digest != receipt.get("artifact_sha256"):
            raise ValueError("Installer digest does not match its build receipt")
        artifacts.append({"name": name, "sha256": digest, "bytes": path.stat().st_size})
    return {
        "schema": "gentle.release_attributes.v1",
        **{key: candidate[key] for key in ("tag", "revision", "cargo_lock_sha256", "workflow_revision", "mode")},
        "linux_distribution": "tarball",
        "generated_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "artifacts": sorted(artifacts, key=lambda item: item["name"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("prepare", "verify", "collect"))
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--artifacts", type=Path, default=Path("dist"))
    args = parser.parse_args()
    if args.command == "prepare":
        record = prepare(args.root, os.environ)
    else:
        if not re.fullmatch(r"[0-9a-f]{64}", os.environ["EXPECTED_LOCK_SHA256"]):
            raise ValueError("A candidate lockfile SHA-256 is required")
        record = validate_checkout(
            args.root, os.environ["RELEASE_TAG"], os.environ["EXPECTED_REVISION"],
            os.environ["EXPECTED_LOCK_SHA256"],
        )
        record.update(workflow_revision=full_sha(os.environ["WORKFLOW_REVISION"]), mode=os.environ["CANDIDATE_MODE"])
        if record["mode"] not in ("validate_only", "publish"):
            raise ValueError("Unsupported candidate mode")
        if record["mode"] == "publish" and tag_revision(args.root, record["tag"]) != record["revision"]:
            raise ValueError("Release tag moved away from the approved candidate SHA")
        if args.command == "collect":
            record = collect_installers(args.artifacts, record)
    encoded = json.dumps(record, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    if args.command == "prepare" and os.environ.get("GITHUB_OUTPUT"):
        with open(os.environ["GITHUB_OUTPUT"], "a") as output:
            for key in ("tag", "revision", "cargo_lock_sha256", "workflow_revision", "mode", "publish"):
                value = record[key]
                output.write(f"{key}={str(value).lower() if isinstance(value, bool) else value}\n")
    print(encoded, end="")


if __name__ == "__main__":
    main()
