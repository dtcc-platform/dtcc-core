"""Tests for demos/smoke_table_cases.py.

Loaded by file path because demos/ is intentionally not a Python package.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

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


# ---------- run_case ----------


def _fake_case(
    name: str = "case_x",
    filename: str = "out.bin",
    params: dict | None = None,
):
    return {
        "name": name,
        "filename": filename,
        "dataset_key": f"table-smoke-{name.replace('_', '-')}",
        "params": params or {"product": "field"},
    }


def test_run_case_purges_stale_case_dir(tmp_path):
    case_dir = tmp_path / "case_x"
    case_dir.mkdir()
    (case_dir / "stale_artifact.bin").write_text("old")
    (case_dir / "stale.manifest.json").write_text("{}")

    dataset = MagicMock()
    dataset.export.return_value = MagicMock()

    exported, published, skipped = script.run_case(
        _fake_case(),
        bounds=[0, 0, 1, 1],
        output_dir=tmp_path,
        publish_config=(None, None, False),
        dataset=dataset,
    )

    assert (exported, published, skipped) == (True, False, False)
    # Case dir was recreated; stale files gone
    assert case_dir.is_dir()
    assert not (case_dir / "stale_artifact.bin").exists()
    assert not (case_dir / "stale.manifest.json").exists()


def test_run_case_calls_export_with_target_bounds_and_params(tmp_path):
    dataset = MagicMock()
    dataset.export.return_value = MagicMock()

    case = _fake_case(
        name="slice_demo",
        filename="smoke_slice.geojson",
        params={"product": "slice", "resolution": 31, "slice_axis": "z"},
    )
    script.run_case(
        case,
        bounds=[10, 20, 30, 40],
        output_dir=tmp_path,
        publish_config=(None, None, False),
        dataset=dataset,
    )

    assert dataset.export.call_count == 1
    call = dataset.export.call_args
    # First positional arg is the target path
    assert call.args[0] == tmp_path / "slice_demo" / "smoke_slice.geojson"
    assert call.kwargs["bounds"] == [10, 20, 30, 40]
    assert call.kwargs["product"] == "slice"
    assert call.kwargs["resolution"] == 31
    assert call.kwargs["slice_axis"] == "z"


def test_run_case_mp4_runtimeerror_is_skipped(tmp_path):
    dataset = MagicMock()
    dataset.export.side_effect = RuntimeError("ffmpeg not found on PATH")

    exported, published, skipped = script.run_case(
        _fake_case(name="streamlines_mp4", filename="smoke_streamlines.mp4"),
        bounds=[0, 0, 1, 1],
        output_dir=tmp_path,
        publish_config=(None, None, False),
        dataset=dataset,
    )

    assert (exported, published, skipped) == (False, False, True)


def test_run_case_non_mp4_runtimeerror_propagates(tmp_path):
    dataset = MagicMock()
    dataset.export.side_effect = RuntimeError("disk full")

    with pytest.raises(RuntimeError, match="disk full"):
        script.run_case(
            _fake_case(name="slice_png", filename="smoke_slice.png"),
            bounds=[0, 0, 1, 1],
            output_dir=tmp_path,
            publish_config=(None, None, False),
            dataset=dataset,
        )


def test_run_case_publish_called_when_enabled(tmp_path):
    package = MagicMock()
    package.publish.return_value = SimpleNamespace(
        dataset_key="table-smoke-case-x", version_number=7
    )
    dataset = MagicMock()
    dataset.export.return_value = package

    exported, published, skipped = script.run_case(
        _fake_case(name="case_x"),
        bounds=[0, 0, 1, 1],
        output_dir=tmp_path,
        publish_config=("http://upload.example", "tok-abc", True),
        dataset=dataset,
    )

    assert (exported, published, skipped) == (True, True, False)
    package.publish.assert_called_once_with(
        dataset_key="table-smoke-case-x",
        upload_url="http://upload.example",
        token="tok-abc",
    )


def test_run_case_publish_skipped_when_disabled(tmp_path):
    package = MagicMock()
    dataset = MagicMock()
    dataset.export.return_value = package

    exported, published, skipped = script.run_case(
        _fake_case(),
        bounds=[0, 0, 1, 1],
        output_dir=tmp_path,
        publish_config=(None, None, False),
        dataset=dataset,
    )

    assert (exported, published, skipped) == (True, False, False)
    package.publish.assert_not_called()
