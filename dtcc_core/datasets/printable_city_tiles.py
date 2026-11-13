"""Printable city tiles and their georeferenced table layout."""

import json
import math
import shutil
from collections.abc import Sequence
from contextlib import suppress
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal

import numpy as np
from pydantic import Field, field_validator
from pyproj import Transformer

from ..builder.geometry_builders.printable_tiles import (
    BASE_HEIGHT_MM,
    BASE_SIZE_MM,
    MAGNET_DEPTH_MM,
    MAGNET_INSET_MM,
    MAGNET_RADIUS_MM,
    MAXIMUM_TILE_SIZE_MM,
    MINIMUM_TILE_SIZE_MM,
    _build_printable_tile,
    _resize_base,
)
from ..common.progress import ProgressTracker, report_progress
from ..model import Mesh, Object, Transform
from .city_surface_mesh import CitySurfaceMeshDataset
from .dataset import DatasetBaseArgs, DatasetDescriptor
from .geospatial import validate_2d_bounds
from .providers import provider_entry

NAME = "printable_city_tiles"
WORKING_CRS = "EPSG:3006"
MAX_TILES = 64
MAX_SPAN_M = 5000
SOURCE_MARGIN_M = 6
BASE_FILE = Path(__file__).with_name("data") / "printable_base_200.stl"


class PrintableCityTilesArgs(DatasetBaseArgs):
    """Area, scale and geometry settings for a printable table model."""

    bounds: Sequence[float] = Field(
        description="Area as [west, south, east, north] in the declared CRS"
    )
    crs: Literal["EPSG:3006", "EPSG:4326"] = Field(
        WORKING_CRS, description="CRS of the bounds, with longitude first in WGS84"
    )
    scale: float = Field(1250, ge=1, description="Denominator of the model scale")
    tile_size_mm: float = Field(
        BASE_SIZE_MM,
        ge=MINIMUM_TILE_SIZE_MM,
        le=MAXIMUM_TILE_SIZE_MM,
        description="Square tile width and grid pitch in millimetres",
    )
    lod: Literal["lod1", "lod2"] = Field("lod1", description="Building level of detail")
    grid_origin: tuple[float, float] | None = Field(
        None, description="Optional fixed grid origin in EPSG:3006 metres"
    )
    anchor_z: float = Field(
        0, description="Reference elevation for base steps in metres"
    )

    @field_validator("bounds")
    @classmethod
    def validate_area(cls, bounds):
        return validate_2d_bounds(bounds)

    @field_validator("scale", "anchor_z")
    @classmethod
    def validate_finite(cls, value):
        if not math.isfinite(value):
            raise ValueError("value must be finite")
        return value

    @field_validator("grid_origin")
    @classmethod
    def validate_origin(cls, value):
        if value is not None and not all(map(math.isfinite, value)):
            raise ValueError("grid_origin must be finite")
        return value


class PrintableCityTilesDataset(DatasetDescriptor):
    """Terrain and buildings in millimetre STLs with matching magnetic bases."""

    name = NAME
    title = "Printable City Tiles"
    description = "City tiles with magnet pockets, a matching base and a map overlay."
    ArgsModel = PrintableCityTilesArgs
    data_category = "derived"
    result_kind = "mesh"
    python_return_type = "dtcc_core.model.Object"
    supported_formats = ["stl", "geojson"]
    multi_file_formats = ("stl",)
    default_crs = WORKING_CRS
    geographic_coverage = "Sweden, subject to terrain and building source coverage"
    update_frequency = "derived on demand from upstream source data"
    provider = [
        provider_entry("lantmateriet", role="source_provider"),
        provider_entry("dtcc-platform", role="processor"),
    ]
    source = [
        {
            "name": "DTCC/Lantmäteriet terrain and building footprints",
            "role": "upstream_dataset",
            "source_terms_status": "requires_review",
        },
        {
            "name": "DTCC Table ribbed base",
            "role": "base_tile_geometry",
            "url": "https://github.com/dtcc-platform/dtcc-table/tree/main/base_tiles",
        },
    ]
    derived_from = [
        {"name": "city_surface_mesh", "relationship": "terrain and buildings"}
    ]
    license = "Review upstream point cloud and footprint terms before redistribution."
    export_capabilities = [
        {
            "product": None,
            "canonical": True,
            "web_formats": ["geojson"],
            "constraints": {},
            "options": {},
            "exclusions": [],
        }
    ]

    def layout(self, **kwargs):
        """Preview the resolved tile footprints as WGS84 GeoJSON without fetching data."""
        args = self.validate(kwargs)
        return _preview(_layout(args), ())

    def describe(self):
        """Describe the Dataset and the formats of its assembly exports."""
        descriptor = super().describe()
        descriptor["serialization"]["bytes_when_format_is_set"] = False
        return descriptor

    def export(self, path, **kwargs):
        """Direct callers to the multi-file export or canonical assembly package."""
        raise ValueError(
            "Use export_stls(folder, **selection) for STLs or "
            "export(path, canonical=True) on the built assembly for a package"
        )

    def export_stls(
        self, directory: str | Path, *, tiles: Object | None = None, **kwargs
    ) -> Path:
        """Write the tile and base STLs with layout.json to a new directory.

        Pass a realized tile assembly as `tiles`, or supply the Dataset
        arguments to build it before export. The output directory must be new.
        """
        if tiles is not None and kwargs:
            raise ValueError("Pass either built tiles or Dataset arguments")
        if tiles is None:
            _print_destination(directory)
            tiles = self(**kwargs)
        return tiles.export_stl(folder=directory)

    def create_context(self, args):
        """Record the assembly CRS alongside the original area parameters."""
        context = super().create_context(args)
        context.metadata.crs = [WORKING_CRS]
        return context

    def build(self, args: PrintableCityTilesArgs):
        layout = _layout(args)
        bounds = layout["bounds"]
        xmin, ymin, xmax, ymax = bounds
        source_bounds = (
            xmin - SOURCE_MARGIN_M,
            ymin - SOURCE_MARGIN_M,
            xmax + SOURCE_MARGIN_M,
            ymax + SOURCE_MARGIN_M,
        )
        with ProgressTracker(phases={"surface": 0.65, "tiles": 0.35}) as progress:
            with progress.phase("surface", "Building the city surface"):
                surface = CitySurfaceMeshDataset()(
                    bounds=source_bounds,
                    lod=args.lod,
                    strict_live=args.strict_live,
                    report_mesh_quality=False,
                    stage_audit_enabled=args.lod == "lod2",
                )
            with progress.phase("tiles", "Building printable tiles"):
                return _realization(layout, args, surface)


def _snap_index(value, origin, pitch, rounding):
    index = (value - origin) / pitch
    nearest = round(index)
    return nearest if abs(index - nearest) <= 1e-9 else rounding(index)


def _layout(args):
    bounds = args.bounds
    if args.crs == "EPSG:4326":
        west, south, east, north = bounds
        if not (-180 <= west < east <= 180 and -90 < south < north < 90):
            raise ValueError("WGS84 bounds must be an ordered longitude/latitude area")
        transformer = Transformer.from_crs(args.crs, WORKING_CRS, always_xy=True)
        bounds = transformer.transform_bounds(*bounds, densify_pts=21, errcheck=True)
    xmin, ymin, xmax, ymax = validate_2d_bounds(bounds)
    pitch = args.tile_size_mm * (args.scale / 1000)
    if pitch > MAX_SPAN_M or max(xmax - xmin, ymax - ymin) > MAX_SPAN_M:
        raise ValueError(
            f"printable layout may span at most {MAX_SPAN_M} metres per axis"
        )
    ox, oy = args.grid_origin or (xmin, ymin)
    col0 = _snap_index(xmin, ox, pitch, math.floor)
    row0 = _snap_index(ymin, oy, pitch, math.floor)
    col1 = max(col0 + 1, _snap_index(xmax, ox, pitch, math.ceil))
    row1 = max(row0 + 1, _snap_index(ymax, oy, pitch, math.ceil))
    columns, rows = col1 - col0, row1 - row0
    if columns * rows > MAX_TILES or max(columns, rows) * pitch > MAX_SPAN_M:
        raise ValueError(
            f"layout is limited to {MAX_TILES} tiles and {MAX_SPAN_M} metres per axis"
        )
    x0, y0 = ox + col0 * pitch, oy + row0 * pitch
    resolved = validate_2d_bounds((x0, y0, x0 + columns * pitch, y0 + rows * pitch))
    tolerance = max(1e-8, pitch * 1e-9)
    if (
        x0 > xmin + tolerance
        or y0 > ymin + tolerance
        or resolved[2] < xmax - tolerance
        or resolved[3] < ymax - tolerance
    ):
        raise ValueError(
            "grid_origin is too far from the selection for precise placement"
        )
    tiles = []
    for row in range(rows):
        for col in range(columns):
            tile_id = f"tile_{col0 + col}_{row0 + row}"
            tiles.append(
                {
                    "id": tile_id,
                    "column": col0 + col,
                    "row": row0 + row,
                    "stl": f"{tile_id}.stl",
                    "bounds": [
                        x0 + col * pitch,
                        y0 + row * pitch,
                        x0 + (col + 1) * pitch,
                        y0 + (row + 1) * pitch,
                    ],
                }
            )
    return {
        "crs": WORKING_CRS,
        "requested_crs": args.crs,
        "requested_bounds": list(args.bounds),
        "bounds": list(resolved),
        "grid_origin": [ox, oy],
        "axes": ["east", "grid north", "up"],
        "scale": args.scale,
        "lod": args.lod,
        "tile_size_mm": args.tile_size_mm,
        "width_mm": columns * args.tile_size_mm,
        "height_mm": rows * args.tile_size_mm,
        "tile_pitch_m": pitch,
        "units": "mm",
        "base_tile": {"stl": "base_tile.stl", "height_mm": BASE_HEIGHT_MM},
        "magnet_pocket": {
            "radius_mm": MAGNET_RADIUS_MM,
            "depth_mm": MAGNET_DEPTH_MM,
            "inset_mm": MAGNET_INSET_MM,
        },
        "tiles": tiles,
    }


def _realization(layout, args, surface):
    from ..io import load_mesh

    x0, y0 = layout["bounds"][:2]
    pitch = layout["tile_pitch_m"]
    size = args.tile_size_mm
    tiles = layout["tiles"]
    xyz = surface.transform(np.asarray(surface.vertices, dtype=float))
    xyz[:, :2] -= [x0, y0]
    triangles = xyz[surface.faces]
    assembly = Object(id=NAME)
    base = _resize_base(load_mesh(BASE_FILE), size)
    base.vertices = base.vertices.astype(float)
    base.offset((0, 0, -BASE_HEIGHT_MM))
    assembly.add_geometry(base, id="base_tile", role="print_template")
    bottoms = []
    for i, tile in enumerate(tiles):
        col = tile["column"] - tiles[0]["column"]
        row = tile["row"] - tiles[0]["row"]
        bounds = (col * pitch, row * pitch, (col + 1) * pitch, (row + 1) * pitch)
        try:
            solid, bottom, buildings = _build_printable_tile(
                triangles,
                np.asarray(surface.markers),
                bounds,
                size,
                args.scale,
                args.anchor_z,
            )
        except ValueError as exc:
            raise ValueError(f"{tile['id']}: {exc}") from exc
        child = Object(id=tile["id"], attributes={"buildings_mm": buildings})
        child.add_geometry(solid, id="print_mesh", lod=args.lod)
        assembly.add_child(child)
        bottoms.append(bottom)
        report_progress(current=i + 1, total=len(tiles), message=f"Built {tile['id']}")
    _place_tiles(assembly, layout, args, bottoms)
    roof_audit = getattr(surface, "stage_audit", {}).get("roof_catalogue", {})
    if "sources" in roof_audit:
        layout["lod1_fallback_count"] = sum(
            source.get("final_outcome") == "lod1_fallback"
            for source in roof_audit["sources"]
        )
    return assembly


def _place_tiles(assembly, layout, args, bottoms):
    x0, y0 = layout["bounds"][:2]
    tiles = layout["tiles"]
    datum = min(bottoms)
    metres_per_mm = args.scale / 1000
    z0 = args.anchor_z + datum * BASE_HEIGHT_MM * metres_per_mm
    affine = np.diag([metres_per_mm, metres_per_mm, metres_per_mm, 1.0])
    affine[:3, 3] = [x0, y0, z0]
    assembly.transform = Transform(srs=WORKING_CRS, affine=affine)
    layout["mounting_elevation_m"] = z0 - BASE_HEIGHT_MM * metres_per_mm
    layout["base_tile"]["count"] = sum(bottom - datum + 1 for bottom in bottoms)
    for child, tile, bottom in zip(
        assembly.get_children(Object), tiles, bottoms, strict=True
    ):
        offset = [
            (tile["column"] - tiles[0]["column"]) * args.tile_size_mm,
            (tile["row"] - tiles[0]["row"]) * args.tile_size_mm,
            (bottom - datum) * BASE_HEIGHT_MM,
        ]
        tile["print_origin_mm"] = offset
        tile["base_elevation_m"] = z0 + offset[2] * metres_per_mm
        tile["base_stack_count"] = bottom - datum + 1
        child.get_geometry(id="print_mesh").offset(offset)
    assembly.attributes["layout"] = layout


def _preview(layout, children):
    from ..io.geojson import feature_collection, polygon_rings, positions

    def feature(rings, **properties):
        polygons = [positions(ring, WORKING_CRS, "EPSG:4326") for ring in rings]
        return {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": polygon_rings(polygons[0], polygons[1:]),
            },
            "properties": properties,
        }

    features = []
    for tile in layout["tiles"]:
        x, y, east, north = tile["bounds"]
        outline = [(x, y), (east, y), (east, north), (x, north)]
        features.append(feature([outline], kind="tile", **tile))
    for child, tile in zip(children, layout["tiles"], strict=False):
        origin = tile["bounds"][:2]
        for polygon in child.attributes["buildings_mm"]:
            rings = [
                np.asarray(ring) * (layout["scale"] / 1000) + origin for ring in polygon
            ]
            features.append(feature(rings, kind="building", tile_id=tile["id"]))
    return feature_collection(features, source_crs=WORKING_CRS, crs="EPSG:4326")


def _print_destination(folder):
    path = Path(folder)
    if path.exists() or path.is_symlink():
        raise FileExistsError(f"Choose a new print output folder: {path}")
    return path


def _export_stl(assembly, folder):
    path = _print_destination(folder)
    _print_layout(assembly)
    path.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=path.parent, prefix=".dtcc-print-") as temporary:
        stage = Path(temporary) / "tiles"
        stage.mkdir()
        _write_print_files(assembly, stage)
        files = sorted(stage.iterdir(), key=lambda file: file.suffix != ".stl")
        published = []
        path.mkdir(exist_ok=False)
        try:
            for file in files:
                destination = path / file.name
                with file.open("rb") as source, destination.open("xb") as target:
                    published.append(destination)
                    shutil.copyfileobj(source, target)
        except OSError:
            for file in published:
                with suppress(OSError):
                    file.unlink()
            with suppress(OSError):
                path.rmdir()
            raise
    return path


def _write_print_files(assembly, directory):
    from ..io import save_mesh

    layout = _print_layout(assembly)
    children = _print_children(assembly, layout)
    (directory / "layout.json").write_text(
        json.dumps(layout, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    base = assembly.get_geometry(id="base_tile")
    save_mesh(
        Mesh(vertices=base.vertices + [0, 0, BASE_HEIGHT_MM], faces=base.faces),
        directory / "base_tile.stl",
    )
    for child, tile in zip(children, layout["tiles"], strict=True):
        mesh = child.get_geometry(id="print_mesh")
        save_mesh(
            Mesh(vertices=mesh.vertices - tile["print_origin_mm"], faces=mesh.faces),
            directory / tile["stl"],
        )


def _print_preview(assembly):
    layout = _print_layout(assembly)
    return _preview(layout, _print_children(assembly, layout))


def _print_children(assembly, layout):
    children = assembly.get_children(Object)
    by_id = {child.id: child for child in children}
    ids = []
    for tile in layout["tiles"]:
        if not isinstance(tile, dict):
            raise ValueError("Invalid tile entry in print layout")
        col, row = tile.get("column"), tile.get("row")
        tile_id = f"tile_{col}_{row}"
        if (
            type(col) is not int
            or type(row) is not int
            or tile.get("id") != tile_id
            or tile.get("stl") != f"{tile_id}.stl"
        ):
            raise ValueError("Invalid tile coordinates, ID or STL filename")
        ids.append(tile_id)
    if (
        len(set(ids)) != len(ids)
        or len(by_id) != len(children)
        or set(ids) != by_id.keys()
    ):
        raise ValueError("Print layout requires one unique mesh per tile ID")
    return [by_id[tile_id] for tile_id in ids]


def _print_layout(assembly):
    layout = getattr(assembly, "attributes", {}).get("layout")
    if (
        not isinstance(layout, dict)
        or not isinstance(layout.get("tiles"), list)
        or not layout.get("tiles")
        or not isinstance(layout.get("base_tile"), dict)
        or layout["base_tile"].get("stl") != "base_tile.stl"
        or not isinstance(assembly.get_geometry(id="base_tile"), Mesh)
    ):
        raise ValueError("Printable city tiles require a realized print layout")
    return layout
