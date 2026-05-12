"""Tests for demos/smoke_table_cases.py.

Loaded by file path because demos/ is intentionally not a Python package.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
_SCRIPT_PATH = _REPO_ROOT / "demos" / "smoke_table_cases.py"
_spec = importlib.util.spec_from_file_location(
    "smoke_table_cases", _SCRIPT_PATH
)
script = importlib.util.module_from_spec(_spec)
sys.modules.setdefault("smoke_table_cases", script)
_spec.loader.exec_module(script)


_EXPECTED_NAMES = [
    "field_vtu",
    "field_pb",
    "slice_geojson",
    "streamlines_geojson",
    "slice_png",
    "streamlines_png",
    "streamlines_mp4",
]

_EXPECTED_FILENAMES = {
    "field_vtu": "smoke_field.vtu",
    "field_pb": "smoke_field.pb",
    "slice_geojson": "smoke_slice.geojson",
    "streamlines_geojson": "smoke_streamlines.geojson",
    "slice_png": "smoke_slice.png",
    "streamlines_png": "smoke_streamlines.png",
    "streamlines_mp4": "smoke_streamlines.mp4",
}


def test_bounds_is_stockholm_500m_square():
    assert script.BOUNDS == [319720, 6397660, 320220, 6398160]


def test_output_dir_is_under_output_smoke_table_cases():
    assert script.OUTPUT_DIR == Path("output/smoke/table_cases")


def test_cases_has_seven_entries_in_expected_order():
    names = [case["name"] for case in script.CASES]
    assert names == _EXPECTED_NAMES


def test_each_case_has_required_keys():
    required = {"name", "filename", "dataset_key", "params"}
    for case in script.CASES:
        assert required <= set(case.keys()), case


def test_filenames_match_spec():
    for case in script.CASES:
        assert case["filename"] == _EXPECTED_FILENAMES[case["name"]]


def test_dataset_keys_follow_table_smoke_pattern():
    for case in script.CASES:
        assert case["dataset_key"] == f"table-smoke-{case['name'].replace('_', '-')}"


def test_mp4_case_passes_format_explicitly():
    mp4 = next(c for c in script.CASES if c["name"] == "streamlines_mp4")
    assert mp4["params"].get("format") == "mp4"


def test_non_mp4_cases_omit_explicit_format():
    for case in script.CASES:
        if case["name"] == "streamlines_mp4":
            continue
        assert "format" not in case["params"], case["name"]


# ---------- resolve_publish_config ----------


def test_resolve_publish_config_both_missing_returns_disabled():
    url, token, enabled = script.resolve_publish_config({})
    assert url is None
    assert token is None
    assert enabled is False


def test_resolve_publish_config_blank_strings_disabled():
    cfg = script.resolve_publish_config(
        {"DTCC_UPLOAD_URL": "", "DTCC_UPLOAD_TOKEN": ""}
    )
    assert cfg == (None, None, False)


def test_resolve_publish_config_whitespace_only_disabled():
    cfg = script.resolve_publish_config(
        {"DTCC_UPLOAD_URL": "   ", "DTCC_UPLOAD_TOKEN": "\t\n"}
    )
    assert cfg == (None, None, False)


def test_resolve_publish_config_only_url_set_disabled():
    url, token, enabled = script.resolve_publish_config(
        {"DTCC_UPLOAD_URL": "http://x"}
    )
    assert url == "http://x"
    assert token is None
    assert enabled is False


def test_resolve_publish_config_only_token_set_disabled():
    url, token, enabled = script.resolve_publish_config(
        {"DTCC_UPLOAD_TOKEN": "tok"}
    )
    assert url is None
    assert token == "tok"
    assert enabled is False


def test_resolve_publish_config_both_set_enabled():
    cfg = script.resolve_publish_config(
        {"DTCC_UPLOAD_URL": "http://x", "DTCC_UPLOAD_TOKEN": "tok"}
    )
    assert cfg == ("http://x", "tok", True)


def test_resolve_publish_config_strips_outer_whitespace():
    cfg = script.resolve_publish_config(
        {"DTCC_UPLOAD_URL": "  http://x  ", "DTCC_UPLOAD_TOKEN": "  tok  "}
    )
    assert cfg == ("http://x", "tok", True)
