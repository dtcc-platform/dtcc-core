"""End-to-end smoke test: drive the CLI's main() directly against the fixture."""
import importlib.util
import sys
from pathlib import Path


FIXTURE_ROOT = Path(__file__).parent / "fixtures" / "minimal_dataset"
REPO_ROOT = Path(__file__).resolve().parents[3]
CLI_PATH = REPO_ROOT / "sandbox" / "evaluate_lod2.py"


def _load_cli_module():
    spec = importlib.util.spec_from_file_location("_eval_lod2_cli", CLI_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules["_eval_lod2_cli"] = module
    spec.loader.exec_module(module)
    return module


def test_cli_runs_end_to_end(tmp_path: Path):
    cli = _load_cli_module()
    cli.main([
        "--dataset", str(FIXTURE_ROOT),
        "--out", str(tmp_path),
    ])
    assert (tmp_path / "per_building.csv").exists()
    assert (tmp_path / "summary.json").exists()
