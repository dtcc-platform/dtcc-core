#!/usr/bin/env python3
"""
Turn a pytest JUnit XML report of live dataset tests into a Markdown table
with one row per test case and its status (passed, failed or skipped) plus
the first line of any failure or skip message.

Usage:
  python scripts/summarize_live_tests.py live-results.xml >> "$GITHUB_STEP_SUMMARY"
"""

from __future__ import annotations

import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

STATUS_ICONS = {"passed": "✅", "failed": "❌", "skipped": "⏭️"}


def _case_status(case: ET.Element) -> tuple[str, str]:
    for tag, status in (("failure", "failed"), ("error", "failed"), ("skipped", "skipped")):
        node = case.find(tag)
        if node is not None:
            message = node.get("message") or (node.text or "")
            return status, message.strip().splitlines()[0] if message.strip() else ""
    return "passed", ""


def _escape(text: str, limit: int = 200) -> str:
    text = text.replace("|", "\\|").replace("\n", " ")
    return text if len(text) <= limit else text[: limit - 1] + "…"


def summarize(report: Path) -> str:
    """Return a Markdown summary of the JUnit XML report."""
    rows = []
    counts = {"passed": 0, "failed": 0, "skipped": 0}
    for case in ET.parse(report).getroot().iter("testcase"):
        status, message = _case_status(case)
        counts[status] += 1
        module = case.get("classname", "").rsplit(".", 1)[-1]
        rows.append(
            f"| {STATUS_ICONS[status]} {status} | `{module}` | `{case.get('name', '')}` "
            f"| {_escape(message)} |"
        )

    lines = [
        "## Live dataset checks",
        "",
        f"{counts['passed']} passed, {counts['failed']} failed, {counts['skipped']} skipped",
        "",
        "| Status | File | Test | Message |",
        "|---|---|---|---|",
        *rows,
    ]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("report", type=Path, help="pytest --junitxml report")
    args = parser.parse_args(argv)
    if not args.report.is_file():
        print(f"## Live dataset checks\n\nNo report found at `{args.report}`.")
        return 0
    sys.stdout.write(summarize(args.report))
    return 0


if __name__ == "__main__":
    sys.exit(main())
