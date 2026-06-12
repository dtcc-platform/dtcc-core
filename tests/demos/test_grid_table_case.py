"""Tests for demos/grid_table_case.py.

Loaded by file path because demos/ is intentionally not a Python package.
"""

from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

_REPO_ROOT = Path(__file__).resolve().parents[2]
_SCRIPT_PATH = _REPO_ROOT / "demos" / "grid_table_case.py"
_spec = importlib.util.spec_from_file_location("grid_table_case", _SCRIPT_PATH)
script = importlib.util.module_from_spec(_spec)
sys.modules.setdefault("grid_table_case", script)
_spec.loader.exec_module(script)

# Mirrors dtcc-upload's DATASET_KEY_RE server contract.
_SERVER_DATASET_KEY_RE = re.compile(r"^[a-z0-9][a-z0-9-]{1,62}[a-z0-9]$")


def test_bounds_matches_smoke_table_cases():
    """The grid must share the table bounds exactly: Atlas only reuses the
    physical calibration for tuple-equal bounds."""
    smoke_path = _REPO_ROOT / "demos" / "smoke_table_cases.py"
    smoke_spec = importlib.util.spec_from_file_location(
        "smoke_table_cases_for_grid_bounds", smoke_path
    )
    smoke = importlib.util.module_from_spec(smoke_spec)
    smoke_spec.loader.exec_module(smoke)

    assert script.BOUNDS == smoke.BOUNDS
    assert script.BOUNDS == [319720, 6397660, 320220, 6398160]


def test_output_dir_is_under_output_calibration_grid_table_case():
    assert script.OUTPUT_DIR == Path("output/calibration_grid/table_case")


def test_case_has_required_keys():
    required = {
        "name",
        "filename",
        "dataset_key",
        "title",
        "description",
        "publish_to_atlas",
        "params",
    }
    assert required <= set(script.CASE.keys())


def test_case_is_table_calibration_grid_geojson():
    assert script.CASE["name"] == "calibration_grid_geojson"
    assert script.CASE["filename"] == "calibration_grid.geojson"
    assert script.CASE["dataset_key"] == "table-calibration-grid-geojson"
    assert _SERVER_DATASET_KEY_RE.fullmatch(script.CASE["dataset_key"])
    assert script.CASE["publish_to_atlas"] is True


def test_case_requests_40_divisions():
    """40 cells map to exactly 1 cm per cell on the 40 cm printed model,
    drawing 41 lines per axis."""
    assert script.CASE["params"]["divisions"] == 40
    # Format comes from the .geojson suffix; no explicit format override.
    assert "format" not in script.CASE["params"]


# ---------- resolve_publish_config ----------


def test_resolve_publish_config_both_missing_returns_disabled():
    assert script.resolve_publish_config({}) == (None, None, False)


def test_resolve_publish_config_blank_strings_disabled():
    cfg = script.resolve_publish_config(
        {"DTCC_UPLOAD_URL": "  ", "DTCC_UPLOAD_TOKEN": ""}
    )
    assert cfg == (None, None, False)


def test_resolve_publish_config_both_set_enabled():
    cfg = script.resolve_publish_config(
        {"DTCC_UPLOAD_URL": "  http://x  ", "DTCC_UPLOAD_TOKEN": "  tok  "}
    )
    assert cfg == ("http://x", "tok", True)


# ---------- run_case ----------


def test_run_case_purges_stale_case_dir(tmp_path):
    case_dir = tmp_path / "calibration_grid_geojson"
    case_dir.mkdir()
    (case_dir / "stale.geojson").write_text("old")

    dataset = MagicMock()
    dataset.export.return_value = MagicMock()

    exported, published = script.run_case(
        script.CASE,
        bounds=[0, 0, 1, 1],
        output_dir=tmp_path,
        publish_config=(None, None, False),
        dataset=dataset,
    )

    assert (exported, published) == (True, False)
    assert case_dir.is_dir()
    assert not (case_dir / "stale.geojson").exists()


def test_run_case_calls_export_with_bounds_title_and_divisions(tmp_path):
    dataset = MagicMock()
    dataset.export.return_value = MagicMock()

    script.run_case(
        script.CASE,
        bounds=[10, 20, 30, 40],
        output_dir=tmp_path,
        publish_config=(None, None, False),
        dataset=dataset,
    )

    call = dataset.export.call_args
    assert call.args[0] == (
        tmp_path / "calibration_grid_geojson" / "calibration_grid.geojson"
    )
    assert call.kwargs["bounds"] == [10, 20, 30, 40]
    assert call.kwargs["title"] == script.CASE["title"]
    assert call.kwargs["description"] == script.CASE["description"]
    assert call.kwargs["divisions"] == 40


def test_run_case_publishes_when_enabled(tmp_path):
    package = MagicMock()
    package.publish.return_value = SimpleNamespace(
        dataset_key="table-calibration-grid-geojson", version_number=2
    )
    dataset = MagicMock()
    dataset.export.return_value = package

    exported, published = script.run_case(
        script.CASE,
        bounds=[0, 0, 1, 1],
        output_dir=tmp_path,
        publish_config=("http://upload.example", "tok-abc", True),
        dataset=dataset,
    )

    assert (exported, published) == (True, True)
    package.publish.assert_called_once_with(
        dataset_key="table-calibration-grid-geojson",
        upload_url="http://upload.example",
        token="tok-abc",
    )


def test_run_case_publish_skipped_when_disabled(tmp_path):
    package = MagicMock()
    dataset = MagicMock()
    dataset.export.return_value = package

    exported, published = script.run_case(
        script.CASE,
        bounds=[0, 0, 1, 1],
        output_dir=tmp_path,
        publish_config=(None, None, False),
        dataset=dataset,
    )

    assert (exported, published) == (True, False)
    package.publish.assert_not_called()


# ---------- main ----------


def test_main_exports_without_publish(tmp_path, capsys):
    dataset = MagicMock()
    dataset.export.return_value = MagicMock()

    exit_code = script.main(
        bounds=[0, 0, 1, 1],
        output_dir=tmp_path,
        env={},
        dataset=dataset,
    )

    assert exit_code == 0
    assert dataset.export.call_count == 1
    assert (tmp_path / "calibration_grid_geojson").is_dir()

    captured = capsys.readouterr()
    assert "DONE: 1 exported, 0 published" in captured.out
    assert "DTCC_UPLOAD_URL or DTCC_UPLOAD_TOKEN not set" in captured.out


def test_main_publishes_when_creds_set(tmp_path, capsys):
    package = MagicMock()
    package.publish.return_value = SimpleNamespace(
        dataset_key="table-calibration-grid-geojson", version_number=1
    )
    dataset = MagicMock()
    dataset.export.return_value = package

    exit_code = script.main(
        bounds=[0, 0, 1, 1],
        output_dir=tmp_path,
        env={
            "DTCC_UPLOAD_URL": "  http://upload.example  ",
            "DTCC_UPLOAD_TOKEN": "  tok-x  ",
        },
        dataset=dataset,
    )

    assert exit_code == 0
    package.publish.assert_called_once_with(
        dataset_key="table-calibration-grid-geojson",
        upload_url="http://upload.example",
        token="tok-x",
    )

    captured = capsys.readouterr()
    assert "DONE: 1 exported, 1 published" in captured.out
    assert "DTCC_UPLOAD_URL or DTCC_UPLOAD_TOKEN not set" not in captured.out
