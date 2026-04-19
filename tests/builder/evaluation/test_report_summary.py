import json
from pathlib import Path

from dtcc_core.builder.evaluation.report import summarize, write_summary_json
from tests.builder.evaluation.test_report_csv import _fake_result


def test_summarize_counts_stage_outcomes():
    results = [
        _fake_result(bid="a", outcome="success"),
        _fake_result(bid="b", outcome="success"),
        _fake_result(bid="c", outcome="insufficient_points"),
    ]
    summary = summarize(results)
    assert summary["count"] == 3
    assert summary["stage_outcome_counts"]["success"] == 2
    assert summary["stage_outcome_counts"]["insufficient_points"] == 1


def test_summarize_reports_roof_type_accuracy_on_success_subset():
    a = _fake_result(bid="a", outcome="success")
    a.roof_type_correct = True
    b = _fake_result(bid="b", outcome="success")
    b.roof_type_correct = False
    c = _fake_result(bid="c", outcome="insufficient_points")
    c.roof_type_correct = False
    summary = summarize([a, b, c])
    assert summary["roof_type_accuracy_on_success"] == 0.5


def test_summarize_reports_watertight_rate():
    a = _fake_result(bid="a", outcome="success")
    a.watertight = True
    b = _fake_result(bid="b", outcome="success")
    b.watertight = False
    summary = summarize([a, b])
    assert summary["watertight_rate"] == 0.5


def test_summarize_timing_percentiles():
    results = []
    for i, ms in enumerate([1.0, 2.0, 3.0, 4.0, 100.0]):
        r = _fake_result(bid=f"b{i}")
        r.total_ms = ms
        results.append(r)
    summary = summarize(results)
    assert summary["total_ms"]["median"] == 3.0
    assert summary["total_ms"]["max"] == 100.0


def test_summary_json_is_written(tmp_path: Path):
    out = tmp_path / "summary.json"
    write_summary_json([_fake_result()], out)
    with out.open() as f:
        data = json.load(f)
    assert data["count"] == 1
