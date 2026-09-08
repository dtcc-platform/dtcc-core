"""Tests for the generated cross-repository dataset contract artifacts."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import zipfile
from pathlib import Path

from dtcc_core.datasets import DatasetManifest


ROOT = Path(__file__).resolve().parents[2]
GENERATOR = ROOT / "scripts" / "generate_dataset_contract.py"
ARTIFACT_NAMES = (
    "contract.json",
    "manifest-v2.schema.json",
    "golden.dtccpkg",
)


def test_contract_artifacts_are_complete_and_reproducible(tmp_path):
    first = tmp_path / "first"
    second = tmp_path / "second"

    _generate(first)
    _generate(second)

    for name in ARTIFACT_NAMES:
        assert (first / name).read_bytes() == (second / name).read_bytes()

    contract = json.loads((first / "contract.json").read_text(encoding="utf-8"))
    expected_names = _clean_registry_names(tmp_path)
    assert contract["schema_version"] == "dtcc-dataset-contract-v1"
    assert contract["dataset_names"] == expected_names
    assert sorted(contract["datasets"]) == expected_names
    assert all(
        contract["datasets"][name]["name"] == name for name in expected_names
    )

    manifest_schema = json.loads(
        (first / "manifest-v2.schema.json").read_text(encoding="utf-8")
    )
    assert manifest_schema == DatasetManifest.model_json_schema()

    with zipfile.ZipFile(first / "golden.dtccpkg") as archive:
        assert sorted(archive.namelist()) == [
            "artifacts/city.json",
            "manifest.json",
        ]
        manifest = DatasetManifest.model_validate_json(archive.read("manifest.json"))
        artifact_bytes = archive.read(manifest.artifacts[0].path)

    assert manifest.identity.name == "city"
    assert manifest.request.dataset_name == "city"
    assert manifest.artifacts[0].size == len(artifact_bytes)
    assert manifest.artifacts[0].sha256 == hashlib.sha256(artifact_bytes).hexdigest()


def _generate(output_directory: Path) -> None:
    subprocess.run(
        [sys.executable, str(GENERATOR), "--output", str(output_directory)],
        check=True,
        cwd=ROOT,
    )


def _clean_registry_names(working_directory: Path) -> list[str]:
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import json; import dtcc_core.datasets as datasets; "
                "print(json.dumps(sorted(datasets.list())))"
            ),
        ],
        check=True,
        cwd=working_directory,
        capture_output=True,
        text=True,
    )
    return json.loads(result.stdout)
