"""Generate the seven-case smoke dataset table pack.

Companion to demos/smoke.py. Organizes the same family of products as
one (product, format) case per subdirectory under output/smoke/table_cases/
so that each case has a clean artifact + manifest pair.

When DTCC_UPLOAD_URL and DTCC_UPLOAD_TOKEN are set in the environment,
each case is also published to a dtcc-upload catalog under a
table-smoke-<name> dataset_key. Otherwise the script writes local files
only and prints an INFO line.

Requires Python 3.12+ (inherited from dtcc_core's tempfile usage).
"""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path
from typing import Any, Mapping, Optional

BOUNDS: list[int] = [319720, 6397660, 320220, 6398160]
OUTPUT_DIR: Path = Path("output/smoke/table_cases")

CASES: list[dict[str, Any]] = [
    {
        "name": "field_vtu",
        "filename": "smoke_field.vtu",
        "dataset_key": "table-smoke-field-vtu",
        "params": {"product": "field", "resolution": 17},
    },
    {
        "name": "field_pb",
        "filename": "smoke_field.pb",
        "dataset_key": "table-smoke-field-pb",
        "params": {"product": "field", "resolution": 17},
    },
    {
        "name": "slice_geojson",
        "filename": "smoke_slice.geojson",
        "dataset_key": "table-smoke-slice-geojson",
        "params": {
            "product": "slice",
            "resolution": 31,
            "slice_axis": "z",
            "slice_position": 0.5,
        },
    },
    {
        "name": "streamlines_geojson",
        "filename": "smoke_streamlines.geojson",
        "dataset_key": "table-smoke-streamlines-geojson",
        "params": {
            "product": "streamlines",
            "streamline_count": 25,
            "streamline_steps": 120,
        },
    },
    {
        "name": "slice_png",
        "filename": "smoke_slice.png",
        "dataset_key": "table-smoke-slice-png",
        "params": {
            "product": "slice",
            "resolution": 128,
            "slice_axis": "z",
            "slice_position": 0.5,
            "profile": "table",
            "width": 1920,
            "height": 1920,
        },
    },
    {
        "name": "streamlines_png",
        "filename": "smoke_streamlines.png",
        "dataset_key": "table-smoke-streamlines-png",
        "params": {
            "product": "streamlines",
            "streamline_count": 64,
            "streamline_steps": 180,
            "profile": "table",
            "width": 1920,
            "height": 1920,
        },
    },
    {
        "name": "streamlines_mp4",
        "filename": "smoke_streamlines.mp4",
        "dataset_key": "table-smoke-streamlines-mp4",
        "params": {
            "product": "streamlines",
            "streamline_count": 48,
            "streamline_steps": 150,
            "profile": "table",
            "format": "mp4",
            "width": 1920,
            "height": 1920,
            "fps": 30,
            "duration": 8.0,
            "period": 8.0,
        },
    },
]


def resolve_publish_config(
    env: Mapping[str, str],
) -> tuple[Optional[str], Optional[str], bool]:
    """Resolve publish credentials from an environment mapping.

    Blank or whitespace-only values are treated as unset, mirroring the
    validation in ``DatasetUploadClient.from_config`` so the local-only
    fallback fires at the gate instead of partway through publishing.
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
) -> tuple[bool, bool, bool]:
    """Execute one case end-to-end.

    Purges the case subdirectory first, exports the artifact via the
    injected dataset module, and optionally publishes the resulting
    package. Returns (exported, published, skipped) flags so the caller
    can keep summary counters.

    A RuntimeError from the MP4 case is caught and reported as skipped;
    any other RuntimeError propagates so the failure is loud.
    """
    case_dir = output_dir / case["name"]
    if case_dir.exists():
        shutil.rmtree(case_dir)
    case_dir.mkdir(parents=True)
    target = case_dir / case["filename"]
    print(f"==> {case['name']}: exporting {target}")

    try:
        package = dataset.export(target, bounds=bounds, **case["params"])
    except RuntimeError as exc:
        if case["name"] == "streamlines_mp4":
            print(f"    SKIPPED ({exc})")
            return False, False, True
        raise

    url, token, publish_enabled = publish_config
    if not publish_enabled:
        return True, False, False

    publication = package.publish(
        dataset_key=case["dataset_key"],
        upload_url=url,
        token=token,
    )
    print(
        f"    published {publication.dataset_key} "
        f"v{publication.version_number}"
    )
    return True, True, False


def main(
    *,
    cases: list[dict[str, Any]] | None = None,
    bounds: list[int] | None = None,
    output_dir: Path | None = None,
    env: Mapping[str, str] | None = None,
    dataset: Any = None,
) -> int:
    """Run every case and print a one-line summary.

    All arguments are injectable for testing. By default we use the
    module-level constants and the live ``dtcc_core.datasets.smoke``
    descriptor; the import is lazy so unit tests can skip the heavy
    dependency chain.
    """
    cases = list(cases) if cases is not None else CASES
    bounds = bounds if bounds is not None else BOUNDS
    output_dir = output_dir if output_dir is not None else OUTPUT_DIR
    env = env if env is not None else os.environ
    if dataset is None:
        import dtcc_core as dtcc

        dataset = dtcc.datasets.smoke

    output_dir.mkdir(parents=True, exist_ok=True)

    publish_config = resolve_publish_config(env)
    _, _, publish_enabled = publish_config
    if not publish_enabled:
        print(
            "INFO: DTCC_UPLOAD_URL or DTCC_UPLOAD_TOKEN not set; "
            "skipping publish step"
        )

    exported = 0
    published = 0
    skipped = 0

    for case in cases:
        e, p, s = run_case(
            case,
            bounds=bounds,
            output_dir=output_dir,
            publish_config=publish_config,
            dataset=dataset,
        )
        exported += int(e)
        published += int(p)
        skipped += int(s)

    print(
        f"DONE: {exported} exported, {published} published, "
        f"{skipped} skipped"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
