#!/usr/bin/env python3
"""Stage tracked desktop resources and smoke extracted native release packages."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import plistlib
import shutil
import subprocess


BINARIES = (
    "gentle", "gentle_cli", "gentle_mcp",
    "gentle_examples_docs", "gentle_publication_report",
)
# Match the Linux distribution inventory, never downloaded/untracked caches.
RESOURCE_PATHS = (
    "assets", "icons", "docs", "test_files", "integrations/python",
    "data/resources/affymetrix/platform_registry.json", "README.md",
    "CONTRIBUTING.md", "copyright", "Credits.rtf",
)


def layout(root: Path, platform: str) -> tuple[Path, Path]:
    if platform == "macos":
        return root / "Contents/MacOS", root / "Contents/Resources"
    return root / "bin", root


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def files(root: Path) -> list[Path]:
    return sorted(p for p in root.rglob("*") if p.is_file() and p != root / "SHA256SUMS")


def stage(repo: Path, binaries: Path, destination: Path, platform: str,
          tag: str, app_bundle: Path | None = None) -> None:
    repo, binaries, destination = repo.resolve(), binaries.resolve(), destination.resolve()
    if destination.exists():
        raise ValueError("Package staging destination must be new")
    suffix = ".exe" if platform == "windows" else ""
    for name in BINARIES:
        source = binaries / (name + suffix)
        if not source.is_file() or source.stat().st_size == 0:
            raise ValueError(f"Missing or empty release binary: {source}")
    tracked = subprocess.check_output(
        ["git", "ls-files", "-z", "--", *RESOURCE_PATHS], cwd=repo,
    ).decode("utf-8").split("\0")
    resources = [Path(p) for p in tracked if p]
    for relative in resources:
        source = repo / relative
        if not source.resolve().is_relative_to(repo) or not source.is_file():
            raise ValueError(f"Unavailable or escaping tracked resource: {relative}")

    app_icon = None
    if platform == "macos":
        if app_bundle is None or not (app_bundle / "Contents/Info.plist").is_file():
            raise ValueError("macOS staging requires a cargo-bundle app with Info.plist")
        app_bundle = app_bundle.resolve()
        with (app_bundle / "Contents/Info.plist").open("rb") as stream:
            icon_name = plistlib.load(stream).get("CFBundleIconFile")
        if icon_name:
            if not isinstance(icon_name, str) or Path(icon_name).name != icon_name:
                raise ValueError("Bundle icon must be a resource filename")
            app_icon = app_bundle / "Contents/Resources" / icon_name
            if not app_icon.suffix:
                app_icon = app_icon.with_suffix(".icns")
            if (not app_icon.is_file()
                    or not app_icon.resolve().is_relative_to(app_bundle / "Contents/Resources")):
                raise ValueError("Missing or escaping declared bundle icon")
        # Replace resources and executables with the reviewed inventory, not
        # cargo-bundle's globs or stale scripting/reproduction binaries.
        shutil.copytree(
            app_bundle, destination, symlinks=True,
            ignore=lambda directory, names: ["Resources", "MacOS"]
            if Path(directory) == app_bundle / "Contents" else [],
        )
    else:
        destination.mkdir(parents=True)
    binary_root, resource_root = layout(destination, platform)
    binary_root.mkdir(parents=True, exist_ok=True)
    if app_icon is not None:
        resource_root.mkdir(parents=True, exist_ok=True)
        shutil.copy2(app_icon, resource_root / app_icon.name)
    for name in BINARIES:
        shutil.copy2(binaries / (name + suffix), binary_root / (name + suffix))
    for relative in resources:
        target = resource_root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo / relative, target)

    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
    (destination / "REVISION").write_text(revision + "\n", encoding="utf-8")
    (destination / "VERSION").write_text(tag + "\n", encoding="utf-8")
    if platform == "linux":
        shutil.copy2(resource_root / "docs/linux_tarball.md", destination / "START-HERE.md")
        libraries = []
        for name in BINARIES:
            result = subprocess.run(["ldd", str(binary_root / name)],
                                    capture_output=True, text=True, check=True)
            if "not found" in result.stdout + result.stderr:
                raise ValueError(f"Unresolved runtime library: {name}")
            libraries.append(f"{name}:\n{result.stdout}{result.stderr}")
        (destination / "runtime-libraries.txt").write_text("\n".join(libraries), encoding="utf-8")
    (destination / "SHA256SUMS").write_text("".join(
        f"{digest(path)}  {path.relative_to(destination).as_posix()}\n"
        for path in files(destination)
    ), encoding="utf-8")


def smoke(root: Path, platform: str, revision: str, checkout: Path) -> None:
    root, checkout = root.resolve(), checkout.resolve()
    if root.is_relative_to(checkout):
        raise ValueError("Smoke the extracted package outside the checkout")
    if (root / "REVISION").read_text().strip() != revision:
        raise ValueError("Package revision does not match selected candidate")
    entries = (root / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
    expected = {}
    for line in entries:
        checksum, relative = line.split("  ", 1)
        path = root / relative
        if relative in expected or not path.resolve().is_relative_to(root):
            raise ValueError(f"Invalid checksum entry: {relative}")
        if not path.is_file() or digest(path) != checksum:
            raise ValueError(f"Package checksum mismatch: {relative}")
        expected[relative] = checksum
    if set(expected) != {p.relative_to(root).as_posix() for p in files(root)}:
        raise ValueError("Package file inventory differs from SHA256SUMS")

    binary_root, resource_root = layout(root, platform)
    suffix = ".exe" if platform == "windows" else ""
    for name in ("gentle_js", "gentle_lua"):
        if (binary_root / (name + suffix)).exists():
            raise ValueError(f"Unexpected packaged scripting binary: {name}")
    for relative in ("assets/genomes.json", "docs/tutorial/manifest.json",
                     "test_files", "integrations/python/gentle_py"):
        resource = resource_root / relative
        if not resource.exists() or (resource.is_file() and resource.stat().st_size == 0):
            raise ValueError(f"Missing packaged resource: {relative}")
    commands = (
        ("gentle", "--version"), ("gentle_cli", "capabilities"),
        ("gentle_mcp", "--help"), ("gentle_examples_docs", "--help"),
        ("gentle_publication_report", "--help"),
        ("gentle_examples_docs", "tutorial-manifest-check"),
    )
    for name, argument in commands:
        executable = binary_root / (name + suffix)
        if not executable.is_file() or executable.stat().st_size == 0:
            raise ValueError(f"Missing or empty packaged binary: {executable}")
        print(f"Package smoke: {name} {argument}", flush=True)
        result = subprocess.run([str(executable), argument], cwd=resource_root,
                                capture_output=True, text=True, encoding="utf-8",
                                errors="replace", timeout=120)
        if result.returncode:
            raise ValueError(f"{name} {argument} failed ({result.returncode}):\n"
                             f"{result.stdout[-4000:]}\n{result.stderr[-4000:]}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    staging = commands.add_parser("stage")
    staging.add_argument("--repo", type=Path, required=True)
    staging.add_argument("--binaries", type=Path, required=True)
    staging.add_argument("--destination", type=Path, required=True)
    staging.add_argument("--tag", required=True)
    staging.add_argument("--app-bundle", type=Path)
    checking = commands.add_parser("smoke")
    checking.add_argument("--root", type=Path, required=True)
    checking.add_argument("--revision", required=True)
    checking.add_argument("--checkout", type=Path, required=True)
    for command in (staging, checking):
        command.add_argument("--platform", choices=("macos", "windows", "linux"), required=True)
    args = vars(parser.parse_args())
    command = args.pop("command")
    (stage if command == "stage" else smoke)(**args)


if __name__ == "__main__":
    main()
