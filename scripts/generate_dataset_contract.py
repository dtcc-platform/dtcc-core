"""Generate the machine-readable dtcc-core dataset contract artifacts."""

from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path

import dtcc_core.datasets as datasets
from dtcc_core.datasets import DatasetManifest, attach_dataset_context
from dtcc_core.model import City


CONTRACT_SCHEMA_VERSION = "dtcc-dataset-contract-v1"
_FIXED_ZIP_TIMESTAMP = (1980, 1, 1, 0, 0, 0)


def generate_contract_artifacts(output_directory: Path) -> None:
    """Write the dataset registry contract, manifest schema, and golden package."""
    output_directory.mkdir(parents=True, exist_ok=True)

    registered = datasets.list()
    dataset_names = sorted(registered)
    descriptions = {}
    for name in dataset_names:
        description = registered[name].describe()
        if description.get("name") != name:
            raise ValueError(
                f"Dataset registry key {name!r} does not match "
                f"describe() name {description.get('name')!r}."
            )
        descriptions[name] = description

    _write_json(
        output_directory / "contract.json",
        {
            "schema_version": CONTRACT_SCHEMA_VERSION,
            "dataset_names": dataset_names,
            "datasets": descriptions,
        },
    )
    _write_json(
        output_directory / "manifest-v2.schema.json",
        DatasetManifest.model_json_schema(),
    )
    _write_golden_package(output_directory / "golden.dtccpkg")


def _write_json(path: Path, payload: dict) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _write_golden_package(path: Path) -> None:
    city = City(id="golden-city")
    args = datasets.city.validate({"bounds": (0.0, 0.0, 1.0, 1.0)})
    attach_dataset_context(city, datasets.city.create_context(args))
    city.export(path, format="json")
    _normalize_zip_metadata(path)


def _normalize_zip_metadata(path: Path) -> None:
    """Make the real export pipeline's ZIP envelope reproducible."""
    with zipfile.ZipFile(path) as archive:
        members = {
            info.filename: archive.read(info.filename)
            for info in archive.infolist()
            if not info.is_dir()
        }

    temporary_path = path.with_name(f"{path.name}.tmp")
    with zipfile.ZipFile(
        temporary_path,
        "w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=9,
    ) as archive:
        for name in sorted(members):
            info = zipfile.ZipInfo(name, date_time=_FIXED_ZIP_TIMESTAMP)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            archive.writestr(info, members[name], compresslevel=9)
    temporary_path.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate dtcc-core dataset contract artifacts."
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Directory in which to write the contract artifacts.",
    )
    args = parser.parse_args()
    generate_contract_artifacts(args.output)


if __name__ == "__main__":
    main()
