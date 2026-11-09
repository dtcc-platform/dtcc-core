"""Build the offline canonical packages that index.html displays.

Writes four directory packages and ``expected_controls.json`` into ``out/``.
"""

import json
import shutil
from pathlib import Path

import matplotlib
import numpy as np
from affine import Affine
from pyproj import Transformer

import dtcc_core.datasets as datasets
from dtcc_core.datasets import attach_dataset_context
from dtcc_core.model import (
    Field,
    FootprintCollection,
    Object,
    Point,
    Raster,
    SensorCollection,
    Surface,
)

OUT = Path(__file__).resolve().parent / "out"
BOUNDS = (319700.0, 6397600.0, 320300.0, 6398200.0)
# Source grid: 40 x 40 pixels of 10 m, rotated 30 degrees about its origin.
GEOREF = (
    Affine.translation(319732.0, 6397973.0)
    * Affine.rotation(30.0)
    * Affine.scale(10.0, -10.0)
)
NODATA = -1
# Label -> (row, column) of the top-left pixel of its 2 x 2 block.
BLOCKS = {2: (6, 8), 3: (6, 28), 4: (18, 18), 5: (28, 6), 6: (30, 30), 7: (12, 14)}


def with_context(model, name):
    # Borrows a registered Dataset's context so the artifact file takes its name;
    # the manifest then carries that Dataset's metadata while the data is synthetic.
    dataset = datasets.get_dataset(name)
    context = dataset.create_context(dataset.validate({"bounds": BOUNDS}))
    attach_dataset_context(model, context)
    return model


def sensors():
    collection = SensorCollection(attributes={"crs": "EPSG:3006"})
    for i, (dx, dy, value) in enumerate(
        [(-150, -100, 4.5), (0, 50, 7.0), (120, 100, 9.5)]
    ):
        point = Point(x=320000.0 + dx, y=6397900.0 + dy, z=10.0)
        point.fields = [
            Field(
                name="air_temperature",
                unit="degC",
                association="sample",
                values=np.array([value]),
            )
        ]
        station = Object(attributes={"station_id": str(i + 1)})
        station.add_geometry(point, "location")
        collection.add_station(station)
    return collection


def footprints():
    x, y = 319900.0, 6397800.0
    outer = np.array(
        [[x, y, 0], [x + 120, y, 0], [x + 120, y + 120, 0], [x, y + 120, 0]]
    )
    hole = np.array(
        [
            [x + 40, y + 40, 0],
            [x + 80, y + 40, 0],
            [x + 80, y + 80, 0],
            [x + 40, y + 80, 0],
        ]
    )
    small = outer + [-180.0, 200.0, 0.0]
    return FootprintCollection(
        footprints=[Surface(vertices=outer, holes=[hole]), Surface(vertices=small)],
        source_ids=["building-1", "building-2"],
        source_indices=[0, 1],
    )


def labels():
    data = np.ones((40, 40), dtype=np.int32)
    data[0:6, 30:40] = NODATA
    for label, (row, col) in BLOCKS.items():
        data[row : row + 2, col : col + 2] = label
    return Raster(data=data, georef=GEOREF, crs="EPSG:3006", nodata=NODATA)


def expected_controls(raster, placement, rendering):
    """Pick, per block, the destination pixel centre that maps nearest its middle.

    The pixel must land at least 0.3 source pixels inside the block's edges.

    The inverse mapping runs from the destination grid in the manifest to the
    source labels; the produced PNG is never read.
    """
    to_source = Transformer.from_crs("EPSG:3857", "EPSG:3006", always_xy=True)
    to_wgs84 = Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)
    grid = Affine(*placement.transform)
    rows, cols = np.mgrid[0 : placement.height, 0 : placement.width]
    x, y = grid * (cols + 0.5, rows + 0.5)
    u, v = ~raster.georef * to_source.transform(x, y)
    colormap = matplotlib.colormaps[rendering.colormap]
    low, high = rendering.value_range
    controls = []
    for label, (row, col) in BLOCKS.items():
        distance = np.maximum(abs(u - (col + 1.0)), abs(v - (row + 1.0)))
        r, c = np.unravel_index(np.argmin(distance), distance.shape)
        assert distance[r, c] < 0.7, f"block {label} has no well-inside pixel"
        assert raster.data[int(np.floor(v[r, c])), int(np.floor(u[r, c]))] == label
        lon, lat = to_wgs84.transform(x[r, c], y[r, c])
        rgba = colormap((label - low) / (high - low), bytes=True)
        controls.append(
            {
                "label": label,
                "lon": lon,
                "lat": lat,
                "rgba": [int(value) for value in rgba],
            }
        )
    return controls


def main():
    OUT.mkdir(exist_ok=True)
    models = {
        "calibration_grid": datasets.calibration_grid(bounds=BOUNDS, divisions=3),
        "sensors": with_context(sensors(), "air_quality"),
        "footprints": with_context(footprints(), "building_footprints"),
        "labels": with_context(labels(), "city"),
    }
    names = []
    for name, model in models.items():
        target = OUT / name
        if target.exists():
            shutil.rmtree(target)
        options = {"format": "png" if name == "labels" else "geojson"}
        package = model.export(target, canonical=True, **options)
        names.append(name)
        if name == "labels":
            [artifact] = [a for a in package.artifacts if a.format == "png"]
            controls = expected_controls(model, artifact.placement, artifact.rendering)
    (OUT / "index.json").write_text(json.dumps(names) + "\n")
    (OUT / "expected_controls.json").write_text(json.dumps(controls, indent=2) + "\n")
    print(f"Wrote {len(names)} packages to {OUT}")


if __name__ == "__main__":
    main()
