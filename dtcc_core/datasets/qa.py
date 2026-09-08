"""Offline QA helpers for registered Dataset v2 definitions."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import fnmatch
import json
from pathlib import Path
import sys
from typing import Any, Iterable, Literal, Sequence

from .dataset import DatasetDescriptor
from .registry import list_datasets


FieldStatus = Literal[
    "present",
    "missing",
    "explicitly_unknown",
    "not_applicable",
    "requires_review",
]
Severity = Literal["info", "warning", "error"]

VALID_FIELD_STATUSES: tuple[str, ...] = (
    "present",
    "missing",
    "explicitly_unknown",
    "not_applicable",
    "requires_review",
)
VALID_SEVERITIES: tuple[str, ...] = ("info", "warning", "error")

_QA_BOUNDS = (0.0, 0.0, 1.0, 1.0)
_ALLOWED_STATIC_REQUIRED_ARGS = {"bounds"}
_REPORT_FORMATS = ("text", "markdown", "json")


@dataclass(frozen=True)
class DatasetQAFinding:
    """One static QA finding for one dataset field."""

    dataset: str
    section: str
    field: str
    status: FieldStatus
    severity: Severity
    message: str


@dataclass(frozen=True)
class DatasetQAReport:
    """Static QA report for one or more registered datasets."""

    findings: tuple[DatasetQAFinding, ...]
    dataset_count: int = 0

    def failures(self, *, strict: bool = False) -> tuple[DatasetQAFinding, ...]:
        """Return findings that should make a CLI/check fail."""
        return tuple(
            finding
            for finding in self.findings
            if finding.severity == "error"
            or (strict and finding.status == "requires_review")
        )

    def has_failures(self, *, strict: bool = False) -> bool:
        """Whether the report has failing findings."""
        return bool(self.failures(strict=strict))

    def status_counts(self) -> dict[str, int]:
        """Return finding counts by QA field status."""
        counts = {status: 0 for status in VALID_FIELD_STATUSES}
        for finding in self.findings:
            counts[finding.status] += 1
        return counts

    def severity_counts(self) -> dict[str, int]:
        """Return finding counts by severity."""
        counts = {severity: 0 for severity in VALID_SEVERITIES}
        for finding in self.findings:
            counts[finding.severity] += 1
        return counts


def audited_datasets(
    *,
    include: str | Sequence[str] | None = None,
    exclude: str | Sequence[str] | None = None,
) -> dict[str, DatasetDescriptor]:
    """Return public built-in ``dtcc-core`` datasets selected for static QA.

    Test datasets, remote descriptors, and third-party registrations are not
    part of this audit. The first implementation pass is intentionally scoped
    to built-in dataset definitions shipped by ``dtcc-core``.
    """

    registered = list_datasets()
    builtins = {
        name: dataset
        for name, dataset in registered.items()
        if _is_builtin_public_dataset(dataset)
    }
    if not builtins:
        raise RuntimeError(
            "No built-in dtcc-core datasets were registered. Import "
            "dtcc_core.datasets before running static dataset QA."
        )

    include_patterns = _normalize_patterns(include)
    exclude_patterns = _normalize_patterns(exclude)
    selected = {
        name: dataset
        for name, dataset in builtins.items()
        if _matches_patterns(name, include_patterns, default=True)
        and not _matches_patterns(name, exclude_patterns, default=False)
    }
    if not selected:
        raise RuntimeError(
            "No built-in dtcc-core datasets matched the include/exclude filters."
        )
    return dict(sorted(selected.items()))


def audit_registered_datasets(
    *,
    include: str | Sequence[str] | None = None,
    exclude: str | Sequence[str] | None = None,
) -> DatasetQAReport:
    """Audit registered public built-in datasets without calling providers."""

    datasets = audited_datasets(include=include, exclude=exclude)
    findings: list[DatasetQAFinding] = []
    for dataset in datasets.values():
        findings.extend(audit_dataset(dataset).findings)
    return DatasetQAReport(tuple(findings), dataset_count=len(datasets))


def audit_dataset(dataset: DatasetDescriptor) -> DatasetQAReport:
    """Audit one dataset descriptor without executing ``dataset.build()``."""

    dataset_name = _dataset_label(dataset)
    findings: list[DatasetQAFinding] = []
    descriptor = _describe_dataset(dataset, dataset_name, findings)
    if descriptor is None:
        return DatasetQAReport(tuple(findings), dataset_count=1)

    _audit_descriptor(dataset_name, descriptor, findings)
    schema = descriptor.get("args_schema")
    required_args = _required_static_args(schema)
    extra_required_args = sorted(required_args - _ALLOWED_STATIC_REQUIRED_ARGS)
    if extra_required_args:
        findings.append(
            DatasetQAFinding(
                dataset=dataset_name,
                section="contract",
                field="static_context",
                status="requires_review",
                severity="warning",
                message=(
                    "Static QA cannot create a DatasetContext from bounds alone; "
                    "required argument(s): " + ", ".join(extra_required_args)
                ),
            )
        )
        return DatasetQAReport(tuple(findings), dataset_count=1)

    try:
        args = dataset.validate({"bounds": _QA_BOUNDS})
    except Exception as exc:  # pragma: no cover - exercised by regression tests
        findings.append(
            DatasetQAFinding(
                dataset=dataset_name,
                section="contract",
                field="static_context",
                status="missing",
                severity="error",
                message=(
                    "Static QA could not validate the minimal bounds-only "
                    f"request used for DatasetContext auditing: {exc}"
                ),
            )
        )
        return DatasetQAReport(tuple(findings), dataset_count=1)

    try:
        context = dataset.create_context(args)
    except Exception as exc:
        findings.append(
            DatasetQAFinding(
                dataset=dataset_name,
                section="contract",
                field="static_context",
                status="missing",
                severity="error",
                message=f"DatasetContext creation failed during static QA: {exc}",
            )
        )
        return DatasetQAReport(tuple(findings), dataset_count=1)

    _add_json_finding(findings, dataset_name, "contract", "context_json", context)
    _audit_context(dataset_name, context, findings)
    return DatasetQAReport(tuple(findings), dataset_count=1)


def format_report(report: DatasetQAReport, format: str = "text") -> str:
    """Format a QA report as text, Markdown, or JSON."""

    if format not in _REPORT_FORMATS:
        raise ValueError(
            f"Unsupported QA report format {format!r}. "
            f"Expected one of: {', '.join(_REPORT_FORMATS)}."
        )
    if format == "json":
        return _format_json_report(report)
    if format == "markdown":
        return _format_markdown_report(report)
    return _format_text_report(report)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the static dataset QA command-line interface."""

    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        report = audit_registered_datasets(
            include=args.include,
            exclude=args.exclude,
        )
        output = format_report(report, format=args.format)
    except RuntimeError as exc:
        parser.error(str(exc))
    except ValueError as exc:
        parser.error(str(exc))

    if args.output:
        Path(args.output).write_text(output + "\n", encoding="utf-8")
    else:
        print(output)

    return 1 if report.has_failures(strict=args.strict) else 0


def _build_parser() -> argparse.ArgumentParser:
    examples = """examples:
  python -m dtcc_core.datasets.qa --format markdown
  python -m dtcc_core.datasets.qa --include smoke --format json
  python -m dtcc_core.datasets.qa --strict --output dataset-qa.md
"""
    parser = argparse.ArgumentParser(
        description=(
            "Audit registered dtcc-core datasets offline and report Dataset v2 "
            "contract, metadata, provenance, and presentation QA findings."
        ),
        epilog=examples,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--format",
        choices=_REPORT_FORMATS,
        default="text",
        help="Report format to print or write. Default: text.",
    )
    parser.add_argument(
        "--output",
        metavar="PATH",
        help="Write the report to PATH instead of stdout.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Treat requires_review findings as failures.",
    )
    parser.add_argument(
        "--include",
        action="append",
        metavar="PATTERN",
        help="Only audit dataset names matching PATTERN. May be repeated.",
    )
    parser.add_argument(
        "--exclude",
        action="append",
        metavar="PATTERN",
        help="Exclude dataset names matching PATTERN. May be repeated.",
    )
    return parser


def _describe_dataset(
    dataset: DatasetDescriptor,
    dataset_name: str,
    findings: list[DatasetQAFinding],
) -> dict[str, Any] | None:
    try:
        descriptor = dataset.describe()
    except Exception as exc:
        findings.append(
            DatasetQAFinding(
                dataset=dataset_name,
                section="contract",
                field="describe",
                status="missing",
                severity="error",
                message=f"dataset.describe() failed during static QA: {exc}",
            )
        )
        return None

    if not isinstance(descriptor, dict):
        findings.append(
            DatasetQAFinding(
                dataset=dataset_name,
                section="contract",
                field="describe",
                status="missing",
                severity="error",
                message=(
                    "dataset.describe() must return a mapping, got "
                    f"{type(descriptor).__name__}."
                ),
            )
        )
        return None
    return descriptor


def _audit_descriptor(
    dataset_name: str,
    descriptor: dict[str, Any],
    findings: list[DatasetQAFinding],
) -> None:
    _add_json_finding(findings, dataset_name, "contract", "descriptor_json", descriptor)
    for field in (
        "name",
        "title",
        "description",
        "data_category",
        "result_kind",
        "python_return_type",
        "args_schema",
        "supported_formats",
    ):
        _add_value_finding(
            findings,
            dataset_name,
            "descriptor",
            field,
            descriptor.get(field),
        )

    schema = descriptor.get("args_schema")
    bounds_property = None
    if isinstance(schema, dict):
        properties = schema.get("properties")
        if isinstance(properties, dict):
            bounds_property = properties.get("bounds")
    _add_value_finding(
        findings,
        dataset_name,
        "args_schema",
        "bounds",
        bounds_property,
    )


def _audit_context(
    dataset_name: str,
    context: Any,
    findings: list[DatasetQAFinding],
) -> None:
    for field in ("name", "title"):
        _add_value_finding(
            findings,
            dataset_name,
            "identity",
            field,
            _get_path(context, "identity", field),
        )

    for field in (
        "description",
        "provider",
        "source",
        "license",
        "crs",
        "data_types",
        "formats",
        "geographic_coverage",
        "update_frequency",
    ):
        _add_value_finding(
            findings,
            dataset_name,
            "metadata",
            field,
            _get_path(context, "metadata", field),
            detect_markers=field in {"source", "license"},
        )

    _add_value_finding(
        findings,
        dataset_name,
        "metadata",
        "collection_period",
        _get_path(context, "metadata", "collection_period"),
        missing_status=_collection_period_missing_status(context),
        detect_markers=True,
    )

    for field in ("sources", "processing_steps", "generated_by"):
        _add_value_finding(
            findings,
            dataset_name,
            "provenance",
            field,
            _get_path(context, "provenance", field),
        )

    for field in ("headline", "summary"):
        _add_value_finding(
            findings,
            dataset_name,
            "presentation",
            field,
            _get_path(context, "presentation", field),
        )

    for field in ("dataset_name", "parameters", "bounds"):
        _add_value_finding(
            findings,
            dataset_name,
            "request",
            field,
            _get_path(context, "request", field),
        )


def _add_json_finding(
    findings: list[DatasetQAFinding],
    dataset_name: str,
    section: str,
    field: str,
    value: Any,
) -> None:
    try:
        json.dumps(_json_value(value), sort_keys=True)
    except TypeError as exc:
        findings.append(
            DatasetQAFinding(
                dataset=dataset_name,
                section=section,
                field=field,
                status="missing",
                severity="error",
                message=f"{section}.{field} is not JSON-serializable: {exc}",
            )
        )
        return
    findings.append(
        DatasetQAFinding(
            dataset=dataset_name,
            section=section,
            field=field,
            status="present",
            severity="info",
            message=f"{section}.{field} is JSON-serializable.",
        )
    )


def _add_value_finding(
    findings: list[DatasetQAFinding],
    dataset_name: str,
    section: str,
    field: str,
    value: Any,
    *,
    missing_status: FieldStatus = "missing",
    detect_markers: bool = False,
) -> None:
    if _has_value(value):
        status = _status_from_present_value(value) if detect_markers else "present"
        findings.append(
            DatasetQAFinding(
                dataset=dataset_name,
                section=section,
                field=field,
                status=status,
                severity=_severity_for_status(status),
                message=_present_message(section, field, status),
            )
        )
        return

    findings.append(
        DatasetQAFinding(
            dataset=dataset_name,
            section=section,
            field=field,
            status=missing_status,
            severity=_severity_for_status(missing_status),
            message=_missing_message(section, field, missing_status),
        )
    )


def _status_from_present_value(value: Any) -> FieldStatus:
    text_values = list(_flatten_text(value))
    normalized = " ".join(text_values).lower().replace("-", " ")
    if "not applicable" in normalized or " n/a " in f" {normalized} ":
        return "not_applicable"
    if "requires review" in normalized or "review " in normalized:
        return "requires_review"
    if "verify " in normalized or "unreviewed" in normalized:
        return "requires_review"
    if "explicitly unknown" in normalized or "unknown" in normalized:
        return "explicitly_unknown"
    return "present"


def _collection_period_missing_status(context: Any) -> FieldStatus:
    update_frequency = str(_get_path(context, "metadata", "update_frequency") or "")
    data_category = str(_get_path(context, "metadata", "data_category") or "")
    if "generated on demand" in update_frequency.lower():
        return "not_applicable"
    if data_category == "simulation":
        return "not_applicable"
    return "requires_review"


def _severity_for_status(status: FieldStatus) -> Severity:
    if status == "missing":
        return "error"
    if status in {"requires_review", "explicitly_unknown"}:
        return "warning"
    return "info"


def _present_message(section: str, field: str, status: FieldStatus) -> str:
    label = f"{section}.{field}"
    if status == "present":
        return f"{label} is populated."
    if status == "requires_review":
        return f"{label} is present but explicitly marked as requiring review."
    if status == "explicitly_unknown":
        return f"{label} is explicitly marked as unknown."
    return f"{label} is explicitly marked as not applicable."


def _missing_message(section: str, field: str, status: FieldStatus) -> str:
    label = f"{section}.{field}"
    if status == "missing":
        return f"Required field {label} is missing."
    if status == "requires_review":
        return f"{label} is not populated; QA review is required."
    if status == "explicitly_unknown":
        return f"{label} is explicitly unknown."
    return f"{label} is explicitly not applicable for this dataset."


def _format_json_report(report: DatasetQAReport) -> str:
    payload = {
        "dataset_count": report.dataset_count,
        "summary": {
            "status": report.status_counts(),
            "severity": report.severity_counts(),
            "failures": len(report.failures()),
        },
        "findings": [asdict(finding) for finding in report.findings],
    }
    return json.dumps(payload, indent=2, sort_keys=True)


def _format_markdown_report(report: DatasetQAReport) -> str:
    lines = [
        "# DTCC Dataset QA Report",
        "",
        f"Audited datasets: {report.dataset_count}",
        "",
        "## Summary",
        "",
        "| Status | Count |",
        "|---|---:|",
    ]
    lines.extend(
        f"| `{status}` | {count} |"
        for status, count in report.status_counts().items()
    )
    lines.extend(
        [
            "",
            "## Findings",
            "",
            "| Dataset | Section | Field | Status | Severity | Message |",
            "|---|---|---|---|---|---|",
        ]
    )
    for finding in report.findings:
        lines.append(
            "| "
            + " | ".join(
                (
                    f"`{_escape_markdown(finding.dataset)}`",
                    _escape_markdown(finding.section),
                    _escape_markdown(finding.field),
                    f"`{finding.status}`",
                    f"`{finding.severity}`",
                    _escape_markdown(finding.message),
                )
            )
            + " |"
        )
    return "\n".join(lines)


def _format_text_report(report: DatasetQAReport) -> str:
    lines = [
        "DTCC Dataset QA Report",
        f"Audited datasets: {report.dataset_count}",
        "",
        "Summary:",
    ]
    lines.extend(
        f"  {status}: {count}" for status, count in report.status_counts().items()
    )
    lines.append("")
    lines.append("Findings:")
    for finding in report.findings:
        lines.append(
            "  "
            f"[{finding.severity}] {finding.dataset} "
            f"{finding.section}.{finding.field}: "
            f"{finding.status} - {finding.message}"
        )
    return "\n".join(lines)


def _escape_markdown(value: str) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def _json_value(value: Any) -> Any:
    if hasattr(value, "model_dump"):
        return value.model_dump(mode="json")
    return value


def _required_static_args(schema: Any) -> set[str]:
    if not isinstance(schema, dict):
        return set()
    required = schema.get("required", ())
    if not isinstance(required, list):
        return set()
    return {str(item) for item in required}


def _has_value(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    if isinstance(value, (list, tuple, set, dict)):
        return bool(value)
    return True


def _flatten_text(value: Any) -> Iterable[str]:
    if value is None:
        return
    if isinstance(value, str):
        yield value
        return
    if isinstance(value, dict):
        for key, item in value.items():
            yield str(key)
            yield from _flatten_text(item)
        return
    if isinstance(value, (list, tuple, set)):
        for item in value:
            yield from _flatten_text(item)
        return
    yield str(value)


def _get_path(value: Any, *path: str) -> Any:
    current = value
    for item in path:
        if isinstance(current, dict):
            current = current.get(item)
        else:
            current = getattr(current, item, None)
        if current is None:
            return None
    return current


def _dataset_label(dataset: DatasetDescriptor) -> str:
    name = getattr(dataset, "name", None)
    if isinstance(name, str) and name.strip():
        return name
    return f"<unnamed:{type(dataset).__name__}>"


def _is_builtin_public_dataset(dataset: DatasetDescriptor) -> bool:
    module = dataset.__class__.__module__
    return (
        module.startswith("dtcc_core.datasets.")
        and module != "dtcc_core.datasets.remote"
    )


def _normalize_patterns(patterns: str | Sequence[str] | None) -> tuple[str, ...]:
    if patterns is None:
        return ()
    if isinstance(patterns, str):
        return (patterns,)
    return tuple(str(pattern) for pattern in patterns)


def _matches_patterns(
    value: str,
    patterns: Sequence[str],
    *,
    default: bool,
) -> bool:
    if not patterns:
        return default
    return any(fnmatch.fnmatchcase(value, pattern) for pattern in patterns)


__all__ = [
    "DatasetQAFinding",
    "DatasetQAReport",
    "FieldStatus",
    "Severity",
    "VALID_FIELD_STATUSES",
    "audit_dataset",
    "audit_registered_datasets",
    "audited_datasets",
    "format_report",
    "main",
]


if __name__ == "__main__":
    sys.exit(main())
