"""Tests for static Dataset QA helpers."""

from __future__ import annotations

import json
from pathlib import Path

from dtcc_core.datasets import DatasetBaseArgs, DatasetDescriptor
from dtcc_core.datasets.qa import (
    VALID_FIELD_STATUSES,
    audit_dataset,
    audit_registered_datasets,
    audited_datasets,
    format_report,
    main as qa_main,
)


class MarkerArgs(DatasetBaseArgs):
    """Bounds-only args for QA marker tests."""


class MarkerDataset(DatasetDescriptor, register=False):
    """Descriptor with explicit non-present QA statuses."""

    name = "marker_dataset"
    description = "Dataset used by static QA tests."
    ArgsModel = MarkerArgs
    data_category = "raw"
    result_kind = "test_result"
    python_return_type = "object"
    supported_formats = ["pb"]
    default_crs = "EPSG:3006"
    provider = [{"name": "Example Provider", "role": "source_provider"}]
    source = ["Explicitly unknown upstream source."]
    license = "Requires review before redistribution."
    collection_period = "Not applicable: deterministic QA fixture."
    geographic_coverage = "Synthetic test bounds"
    update_frequency = "generated on demand"
    processing_steps = ["Create deterministic QA fixture"]
    presentation_summary = "Synthetic QA marker fixture."

    def build(self, args):
        raise AssertionError("static QA must not call build()")


class IncompleteArgs(DatasetBaseArgs):
    """Bounds-only args for missing-field tests."""


class IncompleteDataset(DatasetDescriptor, register=False):
    """Descriptor intentionally missing required QA contract fields."""

    name = "incomplete_dataset"
    ArgsModel = IncompleteArgs

    def build(self, args):
        raise AssertionError("static QA must not call build()")


def test_audit_registered_datasets_is_offline_and_has_no_contract_failures():
    report = audit_registered_datasets()

    assert report.dataset_count == len(audited_datasets())
    assert report.failures() == ()
    assert {finding.status for finding in report.findings} <= set(
        VALID_FIELD_STATUSES
    )
    assert any(finding.status == "requires_review" for finding in report.findings)
    assert any(finding.status == "not_applicable" for finding in report.findings)


def test_audit_dataset_detects_explicit_status_markers_without_building():
    report = audit_dataset(MarkerDataset())

    findings = {
        (finding.section, finding.field): finding for finding in report.findings
    }
    assert findings[("metadata", "source")].status == "explicitly_unknown"
    assert findings[("metadata", "license")].status == "requires_review"
    assert findings[("metadata", "collection_period")].status == "not_applicable"
    assert report.failures() == ()


def test_missing_required_contract_fields_fail_loudly():
    report = audit_dataset(IncompleteDataset())

    failing_fields = {
        (finding.section, finding.field)
        for finding in report.failures()
        if finding.status == "missing"
    }
    assert ("metadata", "provider") in failing_fields
    assert ("metadata", "source") in failing_fields
    assert ("metadata", "license") in failing_fields
    assert ("metadata", "crs") in failing_fields
    assert ("descriptor", "supported_formats") in failing_fields


def test_report_formatters_return_markdown_and_json():
    report = audit_registered_datasets(include="smoke")

    markdown = format_report(report, format="markdown")
    assert "# DTCC Dataset QA Report" in markdown
    assert "| Dataset | Section | Field | Status | Severity | Message |" in markdown
    assert "`smoke`" in markdown

    data = json.loads(format_report(report, format="json"))
    assert data["dataset_count"] == 1
    assert data["summary"]["failures"] == 0
    assert data["findings"]


def test_cli_returns_nonzero_only_for_strict_requires_review(capsys):
    assert qa_main(["--format", "json", "--include", "point_cloud"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data["dataset_count"] == 1
    assert data["summary"]["failures"] == 0
    assert data["summary"]["status"]["requires_review"] > 0

    assert (
        qa_main(["--format", "json", "--include", "point_cloud", "--strict"])
        == 1
    )


def test_qa_matrix_lists_each_builtin_dataset_once():
    matrix_path = Path(__file__).parents[2] / "docs" / "datasets" / "qa-matrix.md"
    rows = []
    for line in matrix_path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("| `"):
            continue
        rows.append(line.split("|")[1].strip().strip("`"))

    assert sorted(rows) == sorted(audited_datasets())
    assert len(rows) == len(set(rows))
