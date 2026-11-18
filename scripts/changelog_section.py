#!/usr/bin/env python3
"""Extract one version's section from CHANGELOG.md.

Used by the release workflow to put the changelog section into the GitHub
Release notes, and usable by hand when drafting a release.

    python scripts/changelog_section.py 0.5.1
    python scripts/changelog_section.py v0.5.1 --output release-notes.md
    python scripts/changelog_section.py --check 0.5.1

Exit status is 0 when the section exists and has content, 2 when it does not.
The heading itself is included by default; pass --no-heading to drop it.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


def extract(text: str, version: str) -> str | None:
    """Return the changelog section for `version`, or None if absent."""
    # Accept "## [0.5.1]", "## [0.5.1] - 2026-11-17", and "## 0.5.1".
    heading = re.compile(
        rf"^##\s+\[?{re.escape(version)}\]?(?:\s*-\s*\S.*)?\s*$", re.MULTILINE
    )
    match = heading.search(text)
    if match is None:
        return None

    start = match.end()
    following = re.search(r"^##\s+", text[start:], re.MULTILINE)
    end = start + following.start() if following else len(text)
    return text[match.start():end].strip() + "\n"


def has_content(section: str) -> bool:
    """True when the section has a non-heading, non-comment line."""
    for line in section.splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract a version section from CHANGELOG.md."
    )
    parser.add_argument("version", help="Version to extract, with or without a leading v.")
    parser.add_argument(
        "--file",
        default="CHANGELOG.md",
        type=Path,
        help="Changelog to read (default: CHANGELOG.md).",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        help="Write the section here instead of stdout.",
    )
    parser.add_argument(
        "--no-heading",
        action="store_true",
        help="Drop the '## [version]' line from the output.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Only report whether a non-empty section exists.",
    )
    args = parser.parse_args()

    version = args.version.lstrip("v")
    if not re.fullmatch(r"\d+\.\d+\.\d+", version):
        parser.error(f"not a semantic version: {args.version!r}")

    if not args.file.is_file():
        parser.error(f"changelog not found: {args.file}")

    section = extract(args.file.read_text(encoding="utf-8"), version)
    if section is None:
        print(f"ERROR: no '## [{version}]' section in {args.file}", file=sys.stderr)
        raise SystemExit(2)
    if not has_content(section):
        print(f"ERROR: the {version} changelog section is empty", file=sys.stderr)
        raise SystemExit(2)

    if args.check:
        print(f"OK: {version} has a changelog section", file=sys.stderr)
        return

    if args.no_heading:
        section = re.sub(r"^##\s+.*?\n+", "", section, count=1)

    if args.output:
        args.output.write_text(section, encoding="utf-8")
        print(f"wrote {args.output}", file=sys.stderr)
    else:
        sys.stdout.write(section)


if __name__ == "__main__":
    main()
