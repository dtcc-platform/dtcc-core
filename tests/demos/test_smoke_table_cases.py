"""Tests for scripts/table_cases/smoke_table_cases.py.

Loaded by file path because scripts/ is intentionally not a Python package.
"""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_SCRIPT_PATH = _REPO_ROOT / "scripts" / "table_cases" / "smoke_table_cases.py"
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

_ATLAS_CASE_NAMES = [
    "slice_geojson",
    "streamlines_geojson",
    "slice_png",
    "streamlines_png",
    "streamlines_mp4",
]

# Mirrors dtcc-upload's DATASET_KEY_RE server contract.
_SERVER_DATASET_KEY_RE = re.compile(r"^[a-z0-9][a-z0-9-]{1,62}[a-z0-9]$")


def test_bounds_is_stockholm_500m_square():
    assert script.BOUNDS == [319720, 6397660, 320220, 6398160]


def test_output_dir_is_under_output_smoke_table_cases():
    assert script.OUTPUT_DIR == Path("output/smoke/table_cases")


def test_cases_has_seven_entries_in_expected_order():
    names = [case["name"] for case in script.CASES]
    assert names == _EXPECTED_NAMES


def test_each_case_has_required_keys():
    required = {
        "name",
        "filename",
        "dataset_key",
        "title",
        "description",
        "publish_to_atlas",
        "params",
    }
    for case in script.CASES:
        assert required <= set(case.keys()), case


def test_filenames_match_spec():
    for case in script.CASES:
        assert case["filename"] == _EXPECTED_FILENAMES[case["name"]]


def test_dataset_keys_follow_table_smoke_pattern():
    for case in script.CASES:
        assert case["dataset_key"] == f"table-smoke-{case['name'].replace('_', '-')}"
        assert _SERVER_DATASET_KEY_RE.fullmatch(case["dataset_key"]), case["name"]


def test_case_titles_are_unique_and_descriptions_are_present():
    titles = [case["title"] for case in script.CASES]
    assert len(set(titles)) == len(titles)
    for case in script.CASES:
        assert case["title"].strip()
        assert case["description"].strip()


def test_publish_to_atlas_subset_is_browser_renderable():
    atlas_cases = [case for case in script.CASES if case["publish_to_atlas"]]
    assert [case["name"] for case in atlas_cases] == _ATLAS_CASE_NAMES

    for case in atlas_cases:
        fmt = case["params"].get("format") or case["filename"].rsplit(".", 1)[1]
        assert fmt in {"geojson", "png", "mp4"}, case["name"]

    local_only_names = [
        case["name"] for case in script.CASES if not case["publish_to_atlas"]
    ]
    assert local_only_names == ["field_vtu", "field_pb"]


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
    publish_to_atlas: bool = True,
):
    return {
        "name": name,
        "filename": filename,
        "dataset_key": f"table-smoke-{name.replace('_', '-')}",
        "title": f"{name.replace('_', ' ').title()} Smoke",
        "description": f"{name} smoke dataset.",
        "publish_to_atlas": publish_to_atlas,
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
    assert call.kwargs["title"] == "Slice Demo Smoke"
    assert call.kwargs["description"] == "slice_demo smoke dataset."
    assert call.kwargs["product"] == "slice"
    assert call.kwargs["resolution"] == 31
    assert call.kwargs["slice_axis"] == "z"


def test_run_case_writes_manifest_title_and_description(tmp_path):
    case = _fake_case(
        name="slice_demo",
        filename="smoke_slice.geojson",
        params={"product": "slice", "resolution": 4},
    )

    class FakeDataset:
        def export(self, target, *, bounds, title, description, **params):
            target.write_text("{}")
            target.with_suffix(".manifest.json").write_text(
                json.dumps({"title": title, "description": description})
            )
            return MagicMock()

    script.run_case(
        case,
        bounds=[0, 0, 1, 1],
        output_dir=tmp_path,
        publish_config=(None, None, False),
        dataset=FakeDataset(),
    )

    manifest = json.loads(
        (tmp_path / "slice_demo" / "smoke_slice.manifest.json").read_text()
    )
    assert manifest["title"] == "Slice Demo Smoke"
    assert manifest["description"] == "slice_demo smoke dataset."


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


def test_run_case_publish_skipped_for_local_only_case_when_enabled(tmp_path):
    package = MagicMock()
    dataset = MagicMock()
    dataset.export.return_value = package

    exported, published, skipped = script.run_case(
        _fake_case(name="field_vtu", publish_to_atlas=False),
        bounds=[0, 0, 1, 1],
        output_dir=tmp_path,
        publish_config=("http://upload.example", "tok-abc", True),
        dataset=dataset,
    )

    assert (exported, published, skipped) == (True, False, False)
    package.publish.assert_not_called()


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


# ---------- main ----------


def test_main_runs_all_seven_cases_without_publish(tmp_path, capsys):
    dataset = MagicMock()
    dataset.export.return_value = MagicMock()

    exit_code = script.main(
        bounds=[0, 0, 1, 1],
        output_dir=tmp_path,
        env={},
        dataset=dataset,
    )

    assert exit_code == 0
    assert dataset.export.call_count == 7

    for case in script.CASES:
        assert (tmp_path / case["name"]).is_dir()

    captured = capsys.readouterr()
    assert "DONE: 7 exported, 0 published, 0 skipped" in captured.out
    assert "DTCC_UPLOAD_URL or DTCC_UPLOAD_TOKEN not set" in captured.out


def test_main_mp4_skip_does_not_stop_other_cases(tmp_path, capsys):
    # `streamlines_mp4` is the last entry in CASES, so a real cases list
    # cannot prove the loop continues after the skip. Pass a custom cases
    # list with cases on either side of the MP4 case so the assertion has
    # something to bite on.
    fake_cases = [
        {
            "name": "before_mp4",
            "filename": "before.bin",
            "dataset_key": "table-smoke-before",
            "title": "Before MP4",
            "description": "Before MP4.",
            "publish_to_atlas": True,
            "params": {"product": "field"},
        },
        {
            "name": "streamlines_mp4",
            "filename": "smoke_streamlines.mp4",
            "dataset_key": "table-smoke-streamlines-mp4",
            "title": "Streamlines MP4",
            "description": "Streamlines MP4.",
            "publish_to_atlas": True,
            "params": {"product": "streamlines", "format": "mp4"},
        },
        {
            "name": "after_mp4",
            "filename": "after.bin",
            "dataset_key": "table-smoke-after",
            "title": "After MP4",
            "description": "After MP4.",
            "publish_to_atlas": True,
            "params": {"product": "field"},
        },
    ]

    dataset = MagicMock()

    def export_side_effect(target, *, bounds, **params):
        if str(target).endswith(".mp4"):
            raise RuntimeError("ffmpeg unavailable in test env")
        return MagicMock()

    dataset.export.side_effect = export_side_effect

    exit_code = script.main(
        cases=fake_cases,
        bounds=[0, 0, 1, 1],
        output_dir=tmp_path,
        env={},
        dataset=dataset,
    )

    assert exit_code == 0
    # All three cases attempted
    assert dataset.export.call_count == 3

    # Subsequent case ran successfully — this is the "doesn't stop" assertion
    assert (tmp_path / "after_mp4").is_dir()

    mp4_dir = tmp_path / "streamlines_mp4"
    assert mp4_dir.is_dir()  # case dir was created and purged
    assert not (mp4_dir / "smoke_streamlines.mp4").exists()

    captured = capsys.readouterr()
    assert "DONE: 2 exported, 0 published, 1 skipped" in captured.out
    assert "SKIPPED" in captured.out


def test_main_publishes_only_atlas_cases_when_creds_set(tmp_path, capsys):
    packages_by_case: dict[str, MagicMock] = {}

    def export_side_effect(target, *, bounds, **params):
        pkg = MagicMock()
        pkg.publish.return_value = SimpleNamespace(
            dataset_key=f"key-{target.parent.name}", version_number=1
        )
        packages_by_case[target.parent.name] = pkg
        return pkg

    dataset = MagicMock()
    dataset.export.side_effect = export_side_effect

    exit_code = script.main(
        bounds=[0, 0, 1, 1],
        output_dir=tmp_path,
        env={"DTCC_UPLOAD_URL": "http://upload.example", "DTCC_UPLOAD_TOKEN": "tok"},
        dataset=dataset,
    )

    assert exit_code == 0
    assert len(packages_by_case) == 7

    for name in _ATLAS_CASE_NAMES:
        packages_by_case[name].publish.assert_called_once()
    for name in ["field_vtu", "field_pb"]:
        packages_by_case[name].publish.assert_not_called()

    captured = capsys.readouterr()
    assert "DONE: 7 exported, 5 published, 0 skipped" in captured.out
    assert "DTCC_UPLOAD_URL or DTCC_UPLOAD_TOKEN not set" not in captured.out


def test_main_passes_stripped_creds_to_published_cases(tmp_path):
    packages_by_case: dict[str, MagicMock] = {}

    def export_side_effect(target, *, bounds, **params):
        pkg = MagicMock()
        pkg.publish.return_value = SimpleNamespace(
            dataset_key="k", version_number=1
        )
        packages_by_case[target.parent.name] = pkg
        return pkg

    dataset = MagicMock()
    dataset.export.side_effect = export_side_effect

    script.main(
        bounds=[0, 0, 1, 1],
        output_dir=tmp_path,
        env={
            "DTCC_UPLOAD_URL": "  http://upload.example  ",
            "DTCC_UPLOAD_TOKEN": "  tok-x  ",
        },
        dataset=dataset,
    )

    expected_keys = {
        case["dataset_key"] for case in script.CASES if case["publish_to_atlas"]
    }
    seen_keys = set()
    for pkg in packages_by_case.values():
        if not pkg.publish.called:
            continue
        kwargs = pkg.publish.call_args.kwargs
        assert kwargs["upload_url"] == "http://upload.example"
        assert kwargs["token"] == "tok-x"
        seen_keys.add(kwargs["dataset_key"])
    assert seen_keys == expected_keys
