"""The live-test summary renders upstream text inertly and survives bad reports."""

import importlib.util
from pathlib import Path

import pytest

REPORT = """<?xml version="1.0" encoding="utf-8"?>
<testsuites><testsuite name="pytest">
  <testcase classname="datasets.live.test_datasets_live" name="test_dataset_live[deso]"/>
  <testcase classname="datasets.live.test_datasets_live" name="test_dataset_live[trees]">
    <failure message="RuntimeError: &lt;b&gt;down&lt;/b&gt; [docs](https://x) ![i](https://y) `code` a|b"/>
  </testcase>
  <testcase classname="datasets.live.test_roads_live" name="test_roads">
    <skipped message="roads: transient upstream failure: timeout"/>
  </testcase>
</testsuite></testsuites>
"""


@pytest.fixture
def summarizer():
    script = Path(__file__).resolve().parents[1] / "scripts/summarize_live_tests.py"
    spec = importlib.util.spec_from_file_location("summarize_live_tests", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_summary_counts_and_code_spans_messages(summarizer, tmp_path):
    report = tmp_path / "live-results.xml"
    report.write_text(REPORT, encoding="utf-8")

    summary = summarizer.summarize(report)

    assert "1 passed, 1 failed, 1 skipped" in summary
    failed_row = next(line for line in summary.splitlines() if "trees" in line)
    message_cell = failed_row.split(" | ")[-1]
    assert message_cell == (
        "`RuntimeError: <b>down</b> [docs](https://x) ![i](https://y) 'code' a\\|b` |"
    )


def test_counts_mode_prints_one_line(summarizer, tmp_path, capsys):
    report = tmp_path / "live-results.xml"
    report.write_text(REPORT, encoding="utf-8")

    assert summarizer.main(["--counts", str(report)]) == 0
    assert capsys.readouterr().out == "1 passed, 1 failed, 1 skipped\n"


@pytest.mark.parametrize("content", ["", "<testsuites><testsuite>"])
@pytest.mark.parametrize("counts", [False, True])
def test_unreadable_report_prints_a_line_instead_of_crashing(
    summarizer, tmp_path, capsys, content, counts
):
    report = tmp_path / "live-results.xml"
    report.write_text(content, encoding="utf-8")

    argv = ["--counts", str(report)] if counts else [str(report)]
    assert summarizer.main(argv) == 0
    out = capsys.readouterr().out
    assert "Could not read the report" in out
    assert "Traceback" not in out


def test_missing_report_prints_a_line(summarizer, tmp_path, capsys):
    assert summarizer.main([str(tmp_path / "missing.xml")]) == 0
    assert "No report found" in capsys.readouterr().out
