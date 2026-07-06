"""Publish the building-footprints alignment layer for the DTCC table.

Operational table-case script for the legacy building-footprints package. The
canonical table catalog generator will migrate this information into
declarative table model specs; until then this script preserves the existing
generation path.

When DTCC_UPLOAD_URL and DTCC_UPLOAD_TOKEN are set in the environment,
the exported package is also published to a dtcc-upload catalog under
the table-footprints-geojson dataset_key. Without upload credentials,
the script writes local files only and prints an INFO line.

Unlike the synthetic smoke cases, this demo downloads live footprint
data and therefore needs network access.
"""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path
from typing import Any, Mapping, Optional

BOUNDS: list[int] = [319720, 6397660, 320220, 6398160]
OUTPUT_DIR: Path = Path("output/footprints/table_case")

CASE: dict[str, Any] = {
    "name": "footprints_geojson",
    "filename": "building_footprints.geojson",
    "dataset_key": "table-footprints-geojson",
    "title": "Building Footprints GeoJSON",
    "description": (
        "Real building footprints over the table bounds for aligning "
        "projected content with the physical model."
    ),
    "publish_to_atlas": True,
    # The table validator rejects anything but EPSG:3006 GeoJSON.
    "params": {"crs": "EPSG:3006"},
}


def resolve_publish_config(
    env: Mapping[str, str],
) -> tuple[Optional[str], Optional[str], bool]:
    """Resolve publish credentials from an environment mapping.

    Blank or whitespace-only values are treated as unset, mirroring
    demos/smoke_table_cases.py and ``DatasetUploadClient.from_config``.
    """
    url = (env.get("DTCC_UPLOAD_URL") or "").strip() or None
    token = (env.get("DTCC_UPLOAD_TOKEN") or "").strip() or None
    return url, token, bool(url and token)


def run_case(
    case: Mapping[str, Any],
    *,
    bounds: list[int],
    output_dir: Path,
    publish_config: tuple[Optional[str], Optional[str], bool],
    dataset: Any,
) -> tuple[bool, bool]:
    """Execute the case end-to-end.

    Purges the case subdirectory first, exports the artifact via the
    injected dataset, and optionally publishes the resulting package.
    Returns (exported, published) flags.
    """
    case_dir = output_dir / case["name"]
    if case_dir.exists():
        shutil.rmtree(case_dir)
    case_dir.mkdir(parents=True)
    target = case_dir / case["filename"]
    print(f"==> {case['name']}: exporting {target}")

    package = dataset.export(
        target,
        bounds=bounds,
        title=case["title"],
        description=case["description"],
        **case["params"],
    )

    url, token, publish_enabled = publish_config
    if not publish_enabled or not case["publish_to_atlas"]:
        return True, False

    publication = package.publish(
        dataset_key=case["dataset_key"],
        upload_url=url,
        token=token,
    )
    print(
        f"    published {publication.dataset_key} "
        f"v{publication.version_number}"
    )
    return True, True


def main(
    *,
    bounds: list[int] | None = None,
    output_dir: Path | None = None,
    env: Mapping[str, str] | None = None,
    dataset: Any = None,
) -> int:
    """Run the case and print a one-line summary.

    All arguments are injectable for testing. By default we use the
    module-level constants and the live ``building_footprints``
    descriptor; the import is lazy so unit tests can skip the heavy
    dependency chain.
    """
    bounds = bounds if bounds is not None else BOUNDS
    output_dir = output_dir if output_dir is not None else OUTPUT_DIR
    env = env if env is not None else os.environ
    if dataset is None:
        import dtcc_core as dtcc

        dataset = dtcc.datasets.building_footprints

    output_dir.mkdir(parents=True, exist_ok=True)

    publish_config = resolve_publish_config(env)
    _, _, publish_enabled = publish_config
    if not publish_enabled:
        print(
            "INFO: DTCC_UPLOAD_URL or DTCC_UPLOAD_TOKEN not set; "
            "skipping publish step"
        )

    exported, published = run_case(
        CASE,
        bounds=bounds,
        output_dir=output_dir,
        publish_config=publish_config,
        dataset=dataset,
    )

    print(f"DONE: {int(exported)} exported, {int(published)} published")
    return 0


if __name__ == "__main__":
    sys.exit(main())
