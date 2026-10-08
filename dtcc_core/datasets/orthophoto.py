"""Lantmäteriet orthophotos through the DTCC LM tile server."""

import math
import re
import sys
import tempfile
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, Optional, Sequence

from pydantic import Field, field_validator, model_validator

from dtcc_core.common.progress import ProgressTracker, report_progress
from dtcc_core.io.raster_tiles import (
    DEFAULT_MAX_MEMORY_BYTES,
    READ_OVERHEAD_BYTES,
    RasterTileLayoutError,
    WorkingMemoryError,
    build_mosaic,
    mosaic_minimum_bytes,
)
from dtcc_core.model import Bounds, Raster
from dtcc_core.model.values.raster_tiles import RasterTileCollection

from .context import attach_dataset_context
from .dataset import DatasetBaseArgs, DatasetDescriptor, DatasetUpstreamError
from .providers import provider_display_name, provider_entry

DEFAULT_SPEKTRALTYP = ("rgb", "rgbi")
BYTES_MESSAGE_STEP = 8 * 1024**2
CONTACT = "the DTCC development team at Chalmers"
MAX_MESSAGE_CHARS = 1024
_LM = provider_display_name("lantmateriet")
_IDENTIFIER = re.compile(r"[a-z0-9_-]+")
_TOO_LARGE = (
    'Use a coarser resolution, product="tiles", or a larger max_memory_bytes.'
)


class OrthophotoArgs(DatasetBaseArgs):
    """Arguments for the orthophoto dataset.

    Only the horizontal extent of ``bounds`` is used; six-value bounds may have a
    zero Z span. ``product="tiles"`` takes no ``resolution``, ``format`` or
    ``max_memory_bytes``, and the default ``product="raster"`` takes rgb and
    rgbi sources only.
    """

    bounds: Sequence[float] = Field(
        ...,
        description=(
            "Bounding box [minx, miny, maxx, maxy] or "
            "[minx, miny, minz, maxx, maxy, maxz] in EPSG:3006"
        ),
    )
    product: Literal["raster", "tiles"] = Field(
        "raster",
        description=(
            "'raster': an RGBA mosaic over the bounds; 'tiles': the original "
            "GeoTIFFs as a file-backed collection"
        ),
    )
    year: Optional[int] = Field(None, description="Acquisition year filter")
    collection: Optional[str] = Field(None, description="STAC collection filter")
    spektraltyp: Optional[list[Literal["rgb", "rgbi", "cir", "pan"]]] = Field(
        None, description="Spectral types; default rgb and rgbi"
    )
    resolution: Optional[float] = Field(
        None,
        description="Mosaic resolution in metres; default the finest source",
    )
    server_url: Optional[str] = Field(
        None,
        description=(
            "Tile server root; default the DTCC_ORTHOPHOTO_URL variable. The "
            "server is internal to Chalmers, reachable from the Chalmers network "
            f"or VPN only; ask {CONTACT} for the address."
        ),
    )
    max_memory_bytes: Optional[int] = Field(
        None, description="Mosaic working-memory budget in bytes; default 2 GiB"
    )
    connect_timeout: float = Field(10.0, description="Connect timeout in seconds")
    read_timeout: float = Field(
        150.0, description="Timeout in seconds for each wait for response bytes"
    )
    format: Optional[Literal["tif"]] = Field(
        None, description="Serialize the mosaic as a GeoTIFF with alpha"
    )

    @field_validator("bounds", mode="before")
    @classmethod
    def validate_bounds(cls, value):
        """Require 4 or 6 finite numbers in order; Z may have zero span."""
        values = list(value)
        if len(values) not in (4, 6):
            raise ValueError("Bounds must be 4 or 6 numbers")
        for number in values:
            if isinstance(number, bool) or not isinstance(number, (int, float)):
                raise ValueError(f"Bounds must be numbers, not {number!r}")
            if not math.isfinite(number):
                raise ValueError("Bounds must be finite")
        if len(values) == 4:
            xmin, ymin, xmax, ymax = values
        else:
            xmin, ymin, zmin, xmax, ymax, zmax = values
            if zmin > zmax:
                raise ValueError("Invalid bounds: zmin <= zmax")
        if not (xmin < xmax and ymin < ymax):
            raise ValueError("Invalid bounds: xmin < xmax, ymin < ymax")
        return [float(number) for number in values]

    @field_validator(
        "year",
        "resolution",
        "max_memory_bytes",
        "connect_timeout",
        "read_timeout",
        mode="before",
    )
    @classmethod
    def reject_bool(cls, value):
        """Booleans are not numbers here."""
        if isinstance(value, bool):
            raise ValueError("must be a number, not a boolean")
        return value

    @field_validator("resolution", "connect_timeout", "read_timeout")
    @classmethod
    def positive_finite(cls, value):
        """Require a finite positive number."""
        if value is not None and not (math.isfinite(value) and value > 0):
            raise ValueError("must be finite and positive")
        return value

    @field_validator("max_memory_bytes")
    @classmethod
    def positive_budget(cls, value):
        """Require a positive budget."""
        if value is not None and value <= 0:
            raise ValueError("must be positive")
        return value

    @field_validator("collection")
    @classmethod
    def collection_identifier(cls, value):
        """Require a server identifier: lower-case letters, digits, - and _."""
        if value is not None and not _IDENTIFIER.fullmatch(value):
            raise ValueError("must match [a-z0-9_-]+")
        return value

    @field_validator("spektraltyp")
    @classmethod
    def distinct_types(cls, value):
        """Require at least one type; drop repeats."""
        if value is None:
            return None
        if not value:
            raise ValueError("must name at least one spectral type")
        return list(dict.fromkeys(value))

    @model_validator(mode="after")
    def product_options(self):
        """Reject options the chosen product does not use."""
        if self.product == "tiles":
            for name in ("resolution", "format", "max_memory_bytes"):
                if getattr(self, name) is not None:
                    raise ValueError(f'{name} does not apply to product="tiles"')
        else:
            unsupported = {"cir", "pan"} & set(self.spektraltyp or ())
            if unsupported:
                raise ValueError(
                    f'product="raster" takes rgb and rgbi sources, not '
                    f"{sorted(unsupported)}; use product=\"tiles\""
                )
        return self


@dataclass
class _Report:
    """What one call listed, used and lost; never stored on the descriptor."""

    product: str
    requested_bounds: tuple[float, float, float, float]
    items_listed: int = 0
    errors: list = field(default_factory=list)
    failed_items: list = field(default_factory=list)
    sources: list = field(default_factory=list)
    grid: Optional[dict[str, Any]] = None


class OrthophotoDataset(DatasetDescriptor):
    """Lantmäteriet orthophotos as an RGBA mosaic or as original tiles."""

    name = "orthophoto"
    title = "Orthophoto"
    description = (
        "Lantmäteriet orthophotos through the DTCC LM tile server. The default "
        "product is an RGBA mosaic (dtcc_core.model.Raster) over the bounds at "
        "the finest source resolution, newest imagery first; product='tiles' "
        "returns the original GeoTIFFs as a RasterTileCollection. The tile "
        "server is internal to Chalmers, reachable from the Chalmers network or "
        "VPN only, and has no public default URL: set DTCC_ORTHOPHOTO_URL or pass "
        f"server_url, and ask {CONTACT} for the address and access. Use "
        "strict_live=True when testing the setup, so a failure raises instead of "
        "returning an empty result."
    )
    ArgsModel = OrthophotoArgs
    data_category = "raw"
    result_kind = "raster"
    python_return_type = "dtcc_core.model.Raster"
    default_crs = "EPSG:3006"
    provider = [
        provider_entry("lantmateriet", role="source_provider"),
        provider_entry("dtcc-platform", role="processor"),
    ]
    source = [
        {
            "name": f"{_LM} orthophotos",
            "role": "source",
            "selected_when": "always; served by the DTCC LM tile server",
            "source_terms_status": "requires_review",
        }
    ]
    license = (
        f"Requires review: verify {_LM} orthophoto terms before redistribution."
    )
    data_types = ["raster", "orthophoto"]
    geographic_coverage = "Sweden"
    update_frequency = f"varies by {_LM} acquisition year and area"
    processing_steps = [
        "List the original orthophotos covering the bounds (one manifest per call)",
        "Download each original GeoTIFF once and reuse it from the local cache",
        "Check each file's header before reading pixels",
        "Mosaic newest first at the finest resolution with nearest-neighbour "
        "sampling; older imagery fills gaps",
        "Write RGB from bands 1-3 and alpha from validity",
    ]
    presentation_headline = "Orthophoto"
    presentation_summary = (
        "Aerial imagery from Lantmäteriet as an RGBA mosaic or original tiles."
    )
    key_points = [
        "Transparent pixels have no imagery; valid black stays opaque.",
        "Mixed acquisition dates and resolutions are possible in one mosaic.",
    ]
    presentation_limitations = [
        "The tile server is internal to Chalmers, reachable from the Chalmers "
        "network or VPN only, with no public default URL; ask "
        f"{CONTACT} for the address and access.",
        "General availability needs a publicly reachable service or another "
        "supported source.",
        "A mosaic exported as a dataset package does not yet mark band 4 as "
        "alpha; use format='tif' for an RGBA GeoTIFF.",
        "Tiles partly imaged without declared nodata can cover older imagery "
        "with black.",
    ]

    def validate(self, kwargs):
        """Validate arguments; a Bounds object keeps its Z values for checking."""
        bounds = kwargs.get("bounds")
        if isinstance(bounds, Bounds):
            kwargs = dict(kwargs)
            kwargs["bounds"] = (
                bounds.xmin,
                bounds.ymin,
                bounds.zmin,
                bounds.xmax,
                bounds.ymax,
                bounds.zmax,
            )
        return super().validate(kwargs)

    def __call__(self, **kwargs):
        args = self.validate(kwargs)
        result, report = self._realize(args)
        result = self.prepare_result(result, args)
        return attach_dataset_context(result, self._context(args, report))

    def build(self, args):
        """Return the mosaic, the tile collection or GeoTIFF bytes."""
        return self._realize(args)[0]

    def _export_payload(self, args):
        """Build once and add this call's context to the export sidecar."""
        result, report = self._realize(args)
        context = self._context(args, report)
        return result, {"dataset_context": context.model_dump(mode="json")}

    def _realize(self, args):
        if args.product == "tiles":
            phases = {"discovery": 0.05, "transfer": 0.95}
        else:
            phases = {"discovery": 0.05, "transfer": 0.6, "headers": 0.05,
                      "mosaic": 0.25}
            if args.format is not None:
                phases["export"] = 0.05
        with ProgressTracker(phases=phases) as tracker:
            return self._realize_tracked(args, tracker)

    def _realize_tracked(self, args, tracker):
        from dtcc_core.io.data import cache as data_cache
        from dtcc_core.io.data import orthophoto as client

        bounds = _horizontal(args.bounds)
        report = _Report(product=args.product, requested_bounds=bounds)
        server_url = client.resolve_server_url(args.server_url)
        timeout = (args.connect_timeout, args.read_timeout)
        budget = args.max_memory_bytes or DEFAULT_MAX_MEMORY_BYTES
        raster = args.product == "raster"
        if raster and args.resolution is not None:
            _refuse_minimum(bounds, args.resolution, budget)
        # Phases are entered outside the client-error handlers, so an error
        # raised by a progress subscriber never reaches them.
        with _phase(tracker, "discovery", "Listing orthophotos"):
            try:
                items = client.fetch_manifest(
                    bounds,
                    server_url=server_url,
                    year=args.year,
                    collection=args.collection,
                    spektraltyp=tuple(args.spektraltyp or DEFAULT_SPEKTRALTYP),
                    timeout=timeout,
                )
            except client.OrthophotoClientError as error:
                upstream = _upstream(error)
                if args.strict_live:
                    raise upstream from error
                report.errors.append(upstream)
                return _empty(args, upstream), report
        report.items_listed = len(items)
        resolutions = [item.resolution for item in items]
        if raster and args.resolution is None and items and None not in resolutions:
            _refuse_minimum(bounds, min(resolutions), budget)
        failures = None if args.strict_live else []
        raised = []
        with _phase(tracker, "transfer", "Downloading original orthophotos"):
            try:
                tiles = client.download_items(
                    items,
                    server_url=server_url,
                    timeout=timeout,
                    cache_root=data_cache.cache_dir,
                    failures=failures,
                    progress=_transfer_progress(items, tracker, failures, raised),
                )
            except client.OrthophotoClientError as error:
                # A subscriber's own error is not an upstream failure.
                if any(error is other for other in raised):
                    raise
                raise _upstream(error) from error
        for item, error in failures or ():
            report.errors.append(_upstream(error))
            report.failed_items.append((item.id, item.datetime, Bounds(*item.bbox)))
        if not raster:
            report.sources = [(tile, None) for tile in tiles]
            return tiles, report
        with ExitStack() as phase:
            phase.enter_context(_phase(tracker, "headers", "Checking headers"))
            stage = ["headers"]
            reporters = {
                name: _Reporter(tracker, name) for name in ("headers", "mosaic")
            }

            def enter_mosaic():
                phase.close()
                phase.enter_context(_phase(tracker, "mosaic", "Compositing"))
                stage[0] = "mosaic"

            def mosaic_progress(event, done, total):
                if event == "mosaic" and stage[0] == "headers":
                    enter_mosaic()
                percent = 100.0 * done / total if total else 100.0
                # The last event of a stage is always reported.
                reporters[event](percent, f"{event}: {done} of {total}", done == total)

            mosaic = build_mosaic(
                tiles.tiles,
                bounds,
                resolution=args.resolution,
                max_memory_bytes=budget,
                progress=mosaic_progress,
            )
            if stage[0] == "headers":
                enter_mosaic()
        source_errors = [_source_error(failure) for failure in mosaic.failures]
        if args.strict_live and source_errors:
            raise source_errors[0]
        report.errors.extend(source_errors)
        report.failed_items.extend(
            (failure.tile.id, failure.tile.datetime, failure.tile.extent)
            for failure in mosaic.failures
        )
        report.sources = list(mosaic.sources)
        report.grid = _grid(mosaic)
        if args.format is None:
            return mosaic.raster, report
        if not mosaic.valid_pixels:
            refusal = (
                "The orthophoto request has no valid pixel, so no GeoTIFF is written"
            )
            if not report.errors:
                raise ValueError(f"{refusal}.")
            first = report.errors[0]
            raise ValueError(f"{refusal}; the first upstream error: {first}") from first
        with _phase(tracker, "export", "Writing the GeoTIFF"):
            with tempfile.TemporaryDirectory() as directory:
                path = _write_geotiff(mosaic, budget, Path(directory))
                # Only the written file is needed from here: release the pixels.
                mosaic = None
                size = _written_size(path)
                if size > budget:
                    raise WorkingMemoryError(
                        size,
                        budget,
                        f"Reading the {size}-byte GeoTIFF",
                        _TOO_LARGE,
                    )
                return _read_written(path), report

    def _context(self, args, report: _Report):
        context = self.create_context(args)
        content = any(pixels != 0 for _, pixels in report.sources)
        if report.errors:
            status = "partial" if content else "failed"
        else:
            status = "complete" if content else "empty"
        health = {
            "status": status,
            "partial_result": bool(report.errors),
            "upstream_error_count": len(report.errors),
            "upstream_errors": [
                self.serialize_upstream_error(error) for error in report.errors
            ],
            "items_listed": report.items_listed,
            "items_used": len(report.sources),
        }
        if report.product == "raster":
            grid = report.grid or {}
            total = grid.get("total_pixels")
            valid = grid.get("valid_pixels", 0)
            health.update(
                requested_bounds=list(report.requested_bounds),
                bounds=grid.get("bounds"),
                resolution=grid.get("resolution"),
                valid_pixels=valid,
                total_pixels=total,
                coverage_complete=bool(total) and valid == total,
            )
        warnings = _warnings(report, health)
        sources = [_source_entry(tile, pixels) for tile, pixels in report.sources]
        dates = sorted(tile.datetime for tile, _ in report.sources)
        period = (
            f"{dates[0].date().isoformat()}/{dates[-1].date().isoformat()}"
            if dates
            else None
        )
        steps = list(context.provenance.processing_steps)
        if report.grid and report.grid.get("resolution") is not None:
            steps.append(
                f"Mosaic grid {report.grid['bounds']} at "
                f"{report.grid['resolution']} m"
            )
        return context.model_copy(
            update={
                "metadata": context.metadata.model_copy(
                    update={"collection_period": period}
                ),
                "provenance": context.provenance.model_copy(
                    update={
                        "sources": list(context.provenance.sources) + sources,
                        "processing_steps": steps,
                    }
                ),
                "presentation": context.presentation.model_copy(
                    update={
                        "warnings": list(context.presentation.warnings) + warnings
                    }
                ),
                "health": health,
                "warnings": warnings,
            }
        )


@contextmanager
def _phase(tracker, name: str, message: str):
    """Run the body in a tracker phase, keeping the body's exception.

    The tracker reports again while closing a phase, so a subscriber that
    raised in the body may raise again there; that second error is dropped.
    """
    phase = tracker.phase(name, message)
    phase.__enter__()
    try:
        yield
    except BaseException:
        try:
            phase.__exit__(*sys.exc_info())
        except BaseException:
            pass
        raise
    phase.__exit__(None, None, None)


class _Reporter:
    """Reports progress in one tracker phase when the call's whole percent changes.

    A report with a different ``detail`` goes out even within the same percent.
    """

    def __init__(self, tracker, name: str):
        phases = tracker.state.phases
        names = list(phases)
        self._start = sum(phases[other].weight for other in names[: names.index(name)])
        self._weight = phases[name].weight
        self._last = None

    def __call__(self, percent: float, message: str, detail=None) -> None:
        key = (int(100 * self._start + self._weight * percent), detail)
        if key == self._last:
            return
        self._last = key
        report_progress(percent=percent, message=message)


def _transfer_progress(items, tracker, failures, raised: list):
    """A download_items callback reporting per item, and per byte with a length.

    The phase advances by item, and within an item by the share of its
    Content-Length received; without one it holds and the message gives the
    bytes received. Reports go out when the call's whole percent changes and,
    without a Content-Length, every BYTES_MESSAGE_STEP bytes; the last one
    counts the originals ready, leaving out the ``failures`` recorded. Errors
    raised while reporting are collected in ``raised``.
    """
    report = _Reporter(tracker, "transfer")

    def progress(index, count, written, expected):
        fraction = min(written, expected) / expected if expected else 0.0
        percent = 100.0 * (index + fraction) / count if count else 100.0
        if index < count:
            item = items[index]
            message = f"{item.collection}/{item.id}: {written / 1024**2:.1f} MiB"
            detail = written // BYTES_MESSAGE_STEP if expected is None else 0
        else:
            ready = count - len(failures or ())
            message = f"{ready} of {count} original orthophotos ready"
            detail = "ready"
        try:
            report(percent, message, detail)
        except BaseException as error:
            raised.append(error)
            raise

    return progress


def _horizontal(bounds) -> tuple[float, float, float, float]:
    if len(bounds) == 6:
        return (bounds[0], bounds[1], bounds[3], bounds[4])
    return tuple(bounds)


def _refuse_minimum(bounds, resolution, budget) -> None:
    minimum = mosaic_minimum_bytes(bounds, resolution)
    if minimum > budget:
        raise WorkingMemoryError(
            minimum,
            budget,
            f"An orthophoto mosaic at {resolution} m over {list(bounds)}",
            _TOO_LARGE,
        )


def _empty(args, error: DatasetUpstreamError):
    """The result of a call whose item list failed with ``error``."""
    if args.product == "tiles":
        return RasterTileCollection()
    if args.format is not None:
        raise ValueError(
            "The orthophoto request has no valid pixel, so no GeoTIFF is "
            f"written; the item list could not be fetched: {error}"
        ) from error
    return Raster(crs="EPSG:3006")


def _upstream(error) -> DatasetUpstreamError:
    return DatasetUpstreamError(
        dataset="orthophoto",
        operation=error.operation,
        target=error.target,
        failure_class=error.failure_class,
        status_code=error.status_code,
        message=str(error)[:MAX_MESSAGE_CHARS],
    )


def _source_error(failure) -> DatasetUpstreamError:
    """A mosaic source failure, described without the exception's text.

    GDAL messages name the cached file, so only the item, the layout reason
    (which names no path) and the pixels kept are reported.
    """
    tile = failure.tile
    target = f"{tile.collection}/{tile.id}"
    if isinstance(failure.error, RasterTileLayoutError):
        reason = failure.error.reason
    else:
        reason = "the file could not be decoded"
    message = (
        f"{target}: {reason}; {failure.pixels} pixels were used before the "
        "failure"
    )
    return DatasetUpstreamError(
        dataset="orthophoto",
        operation="read",
        target=target,
        failure_class="invalid_payload",
        status_code=None,
        message=message[:MAX_MESSAGE_CHARS],
    )


def _grid(mosaic) -> dict[str, Any]:
    data = mosaic.raster.data
    total = data.shape[0] * data.shape[1] if data.ndim == 3 else None
    return {
        "bounds": None if mosaic.bounds is None else list(mosaic.bounds),
        "resolution": mosaic.resolution,
        "valid_pixels": mosaic.valid_pixels,
        "total_pixels": total,
    }


def _source_entry(tile, pixels) -> dict[str, Any]:
    entry = {
        "name": f"{tile.collection}/{tile.id}",
        "role": "source_item",
        "collection": tile.collection,
        "id": tile.id,
        "datetime": tile.datetime.isoformat(),
        "spektraltyp": tile.spektraltyp,
        "resolution": tile.resolution,
    }
    if pixels is not None:
        entry["pixels"] = pixels
    return entry


def _overlaps(a, b) -> bool:
    return (
        a.xmin < b.xmax and b.xmin < a.xmax and a.ymin < b.ymax and b.ymin < a.ymax
    )


def _warnings(report: _Report, health) -> list[str]:
    warnings = [
        f"{error.target}: {error.failure_class}: {error.message}"
        for error in report.errors
    ]
    used = [tile for tile, pixels in report.sources if pixels != 0]
    for failed_id, failed_datetime, failed_extent in report.failed_items:
        older = [
            tile
            for tile in used
            if tile.datetime < failed_datetime and _overlaps(tile.extent, failed_extent)
        ]
        if older:
            warnings.append(
                f"Older imagery ({', '.join(t.id for t in older)}) overlaps "
                f"{failed_id}, which failed."
            )
    if report.product == "raster":
        total = health.get("total_pixels")
        if total and not health["coverage_complete"]:
            share = 1 - health["valid_pixels"] / total
            warnings.append(f"{share:.1%} of the mosaic has no imagery (transparent).")
    if not report.sources and not report.errors:
        warnings.append("No orthophoto has imagery in the requested bounds.")
    dates = {tile.datetime.date() for tile in used}
    resolutions = {tile.resolution for tile in used}
    if len(dates) > 1:
        warnings.append("The imagery has mixed acquisition dates.")
    if len(resolutions) > 1:
        warnings.append("The imagery has mixed source resolutions.")
    return warnings


def _write_geotiff(mosaic, budget, directory: Path) -> Path:
    """Write the mosaic as an RGBA GeoTIFF after checking the write budget."""
    import rasterio

    from dtcc_core.io import raster as raster_io

    height, width = mosaic.raster.data.shape[:2]
    needed = (
        height * width * 4
        + min(raster_io.WRITE_STRIP_BYTES, height * width * 4)
        + READ_OVERHEAD_BYTES
    )
    if needed > budget:
        raise WorkingMemoryError(
            needed, budget, f"Writing a {width} x {height} GeoTIFF", _TOO_LARGE
        )
    path = directory / "orthophoto.tif"
    with rasterio.Env(GDAL_CACHEMAX=READ_OVERHEAD_BYTES):
        raster_io._save_geotif(mosaic.raster, path, alpha=True)
    return path


def _written_size(path: Path) -> int:
    return path.stat().st_size


def _read_written(path: Path) -> bytes:
    return path.read_bytes()
