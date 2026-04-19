from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator, Optional

import fiona
import numpy as np
from shapely.geometry import shape as shapely_shape

from dtcc_core.model.geometry.pointcloud import PointCloud
from dtcc_core.model.geometry.surface import MultiSurface, Surface


@dataclass
class GroundTruth:
    """Ground-truth attributes for a single building in the eval dataset."""

    roof_type: str
    ground_height: float = 0.0
    ridge_height: Optional[float] = None
    eave_height: Optional[float] = None
    geometry: Optional[MultiSurface] = None
    expected_plane_count: Optional[int] = None
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass
class EvalSample:
    """A single building (input + ground truth) for the eval harness to process."""

    building_id: str
    footprint: Surface
    point_cloud: PointCloud
    ground_truth: GroundTruth


class EvalDataset:
    """Iterates over buildings on disk using a folder-per-building convention.

    Directory layout:
        <root>/
          <building_id>/
            footprint.geojson   OR   footprint.shp   OR   footprint.gpkg
            points.npy          OR   points.las      OR   points.laz
            ground_truth.json
    """

    def __init__(self, root: Path | str):
        self._root = Path(root)
        if not self._root.is_dir():
            raise FileNotFoundError(f"Dataset root not found: {self._root}")
        self._building_dirs = sorted(
            d for d in self._root.iterdir() if d.is_dir()
        )

    def __len__(self) -> int:
        return len(self._building_dirs)

    def __iter__(self) -> Iterator[EvalSample]:
        for bdir in self._building_dirs:
            yield self._load(bdir)

    def _load(self, bdir: Path) -> EvalSample:
        gt = self._load_ground_truth(bdir / "ground_truth.json")
        footprint = self._load_footprint(bdir, gt)
        points = self._load_points(bdir)
        pc = PointCloud(points=points)
        return EvalSample(
            building_id=bdir.name,
            footprint=footprint,
            point_cloud=pc,
            ground_truth=gt,
        )

    @staticmethod
    def _load_ground_truth(path: Path) -> GroundTruth:
        with open(path) as f:
            raw = json.load(f)
        known = {
            "roof_type", "ground_height", "ridge_height",
            "eave_height", "expected_plane_count",
        }
        return GroundTruth(
            roof_type=raw["roof_type"],
            ground_height=float(raw.get("ground_height", 0.0)),
            ridge_height=raw.get("ridge_height"),
            eave_height=raw.get("eave_height"),
            expected_plane_count=raw.get("expected_plane_count"),
            extra={k: v for k, v in raw.items() if k not in known},
        )

    @staticmethod
    def _load_footprint(bdir: Path, gt: GroundTruth) -> Surface:
        for candidate in ("footprint.geojson", "footprint.shp", "footprint.gpkg"):
            fp_path = bdir / candidate
            if fp_path.exists():
                break
        else:
            raise FileNotFoundError(f"No footprint found in {bdir}")
        z = gt.ground_height
        with fiona.open(fp_path) as src:
            feat = next(iter(src))
            geom = shapely_shape(feat["geometry"])
            coords = list(geom.exterior.coords)
            if coords[0] == coords[-1]:
                coords = coords[:-1]
            verts = np.array([[x, y, z] for x, y in coords], dtype=float)
        return Surface(vertices=verts)

    @staticmethod
    def _load_points(bdir: Path) -> np.ndarray:
        npy = bdir / "points.npy"
        if npy.exists():
            return np.load(npy)
        for ext in ("points.las", "points.laz"):
            las_path = bdir / ext
            if las_path.exists():
                import laspy
                with laspy.open(las_path) as src:
                    las = src.read()
                    return np.column_stack([las.x, las.y, las.z])
        raise FileNotFoundError(f"No point cloud found in {bdir}")
