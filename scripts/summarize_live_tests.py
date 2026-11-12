#!/usr/bin/env python3
"""
Turn a pytest JUnit XML report of live dataset tests into a Markdown table
with one row per test case and its status (passed, failed or skipped) plus
the first line of any failure or skip message.

Usage:
  python scripts/summarize_live_tests.py live-results.xml >> "$GITHUB_STEP_SUMMARY"
  python scripts/summarize_live_tests.py --counts live-results.xml
"""

from __future__ import annotations

import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

STATUS_ICONS = {"passed": "✅", "failed": "❌", "skipped": "⏭️"}


def _case_status(case: ET.Element) -> tuple[str, str]:
    for tag, status in (
        ("failure", "failed"),
        ("error", "failed"),
        ("skipped", "skipped"),
    ):
        node = case.find(tag)
        if node is not None:
            message = node.get("message") or (node.text or "")
            return status, message.strip().splitlines()[0] if message.strip() else ""
    return "passed", ""


def _code_span(text: str, limit: int = 200) -> str:
    """Show upstream text verbatim in a table cell, without Markdown or HTML."""
    text = " ".join(text.split()).replace("`", "'").replace("|", "\\|")
    if len(text) > limit:
        text = text[: limit - 1] + "…"
    return f"`{text}`" if text else ""


def _read_cases(report: Path) -> list[tuple[str, str, str, str]]:
    """Return (status, module, test name, message) for each test case."""
    cases = []
    for case in ET.parse(report).getroot().iter("testcase"):
        status, message = _case_status(case)
        module = case.get("classname", "").rsplit(".", 1)[-1]
        cases.append((status, module, case.get("name", ""), message))
    return cases


def _count_line(cases: list[tuple[str, str, str, str]]) -> str:
    counts = {"passed": 0, "failed": 0, "skipped": 0}
    for status, *_ in cases:
        counts[status] += 1
    return (
        f"{counts['passed']} passed, {counts['failed']} failed, "
        f"{counts['skipped']} skipped"
    )


def summarize(report: Path) -> str:
    """Return a Markdown summary of the JUnit XML report."""
    cases = _read_cases(report)
    rows = [
        f"| {STATUS_ICONS[status]} {status} | `{module}` | `{name}` "
        f"| {_code_span(message)} |"
        for status, module, name, message in cases
    ]
    lines = [
        "## Live dataset checks",
        "",
        _count_line(cases),
        "",
        "| Status | File | Test | Message |",
        "|---|---|---|---|",
        *rows,
    ]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("report", type=Path, help="pytest --junitxml report")
    parser.add_argument(
        "--counts",
        action="store_true",
        help="print only the passed/failed/skipped counts as one line",
    )
    args = parser.parse_args(argv)

    problem = None
    if not args.report.is_file():
        problem = f"No report found at `{args.report}`."
    else:
        try:
            output = (
                _count_line(_read_cases(args.report)) + "\n"
                if args.counts
                else summarize(args.report)
            )
        except ET.ParseError as exc:
            problem = f"Could not read the report at `{args.report}`: {exc}."

    if problem is not None:
        if args.counts:
            print(problem.replace("`", ""))
        else:
            print(f"## Live dataset checks\n\n{problem}")
        return 0
    sys.stdout.write(output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
