"""Cleaning handoff and phase measurements for the city benchmark.

The cleaning phase saves the raw city, the cleaned city (``clean_footprints``)
and the cleaned footprints. The meshing phase meshes the cleaned city, whose
cleaning record lets the mesh builders skip cleaning. Saved coordinates are in
EPSG:3006 (metres), never GeoJSON's implied WGS84.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from affine import Affine
from shapely import make_valid
from shapely.geometry import mapping, shape
from shapely.ops import unary_union

from dtcc_core.model import Bounds, City, GeometryType, Raster
from dtcc_core.io.model import load_model, save_model
from dtcc_core.builder.cleaning.clean import MERGED_FROM_KEY, RECORD_KEY, cleaning_record
from dtcc_core.builder.geometry_builders import meshes
from dtcc_core.datasets._city_mesh_common import (
    prepare_city_from_bounds,
    prepare_footprint_city_from_bounds,
)
from benchmarks.benchmark_catalog import CLEANING_PARAMETERS

# Format version of saved cleaning artifacts; loading rejects any other.
ARTIFACT_VERSION = 4
PREPARATION_PARAMETERS = (
    "raster_cell_size",
    "raster_radius",
    "remove_outliers",
    "outlier_threshold",
)


def write_json(path: Path, value) -> None:
    def json_ready(item):
        if isinstance(item, dict):
            return {str(key): json_ready(child) for key, child in item.items()}
        if isinstance(item, (list, tuple)):
            return [json_ready(child) for child in item]
        if isinstance(item, np.ndarray):
            return json_ready(item.tolist())
        if isinstance(item, np.generic):
            return item.item()
        return item

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(json_ready(value), indent=2, allow_nan=False) + "\n"
    )
    temporary.replace(path)


def footprint_metrics(raw, cleaned, source_map, scale, *, epsilon=None):
    """Measure output and drift independently of repair operators/trace keys."""
    from dtcc_core.builder.cleaning.contract import check_cleaning_contract

    raw_union = unary_union([make_valid(p) for p in raw])
    clean_union = unary_union(cleaned)
    represented = {i for group in source_map for i in group}
    return {
        "geometric_contract": (
            check_cleaning_contract(raw, cleaned, delta=scale, epsilon=epsilon)
            if scale > 0
            else {"status": "not_checked", "reason": "No positive resolution declared"}
        ),
        "input_count": len(raw),
        "output_count": len(cleaned),
        "invalid_input_count": sum(not p.is_valid for p in raw),
        "invalid_output_count": sum(not p.is_valid for p in cleaned),
        "unrepresented_source_count": len(raw) - len(represented),
        "input_hole_count": sum(len(p.interiors) for p in raw),
        "output_hole_count": sum(len(p.interiors) for p in cleaned),
        "input_area": raw_union.area,
        "output_area": clean_union.area,
        "removed_area": raw_union.difference(clean_union).area,
        "added_area": clean_union.difference(raw_union).area,
        "boundary_displacement_max": (
            raw_union.boundary.hausdorff_distance(clean_union.boundary)
            if not raw_union.is_empty and not clean_union.is_empty
            else None
        ),
        "overlap_area": max(0.0, sum(p.area for p in cleaned) - clean_union.area),
        "vertex_count": sum(
            len(ring.coords) - 1
            for polygon in cleaned
            for ring in [polygon.exterior, *polygon.interiors]
        ),
    }


def _footprint_polygon(building):
    footprint = building.get_geometry(GeometryType.LOD0)
    return None if footprint is None else footprint.to_polygon(simplify=0.0)


def _raw_source_map(raw_city, cleaned_city, kept):
    """Raw building indices behind each kept cleaned building."""
    index = {building.id: i for i, building in enumerate(raw_city.buildings)}
    source_map = []
    for position in kept:
        building = cleaned_city.buildings[position]
        ids = building.attributes.get(MERGED_FROM_KEY) or [building.id]
        source_map.append(sorted({index[i] for i in ids if i in index}))
    return source_map


def save_cleaning(directory, raw_city, cleaned_city, footprints, task):
    """Save the raw and cleaned cities and the cleaning stage metadata."""
    directory.mkdir(parents=True, exist_ok=True)
    save_model(raw_city, directory / "raw.dtcc")
    save_model(cleaned_city, directory / "cleaned.dtcc")
    payload = {
        "type": "FeatureCollection",
        "metadata": {
            "version": ARTIFACT_VERSION,
            "crs": "EPSG:3006",
            "bounds": task["case"]["bounds"],
            "parameters": {k: task["parameters"][k] for k in CLEANING_PARAMETERS},
            "declared_scale": footprints.declared_scale,
            "cleaning": cleaned_city.attributes[RECORD_KEY],
            "before_selection_contract": footprints.diagnostics.get(
                "before_selection_contract"
            ),
            "selection": footprints.diagnostics.get("selection"),
        },
        "features": [
            {
                "type": "Feature",
                "geometry": mapping(surface.to_polygon(simplify=0.0)),
                "properties": {
                    "building_index": indices[0],
                    "building_id": cleaned_city.buildings[indices[0]].id,
                },
            }
            for surface, indices in zip(footprints.surfaces, footprints.source_map)
        ],
    }
    path = directory / "footprints.geojson"
    write_json(path, payload)
    return path


def load_cleaning(path, task):
    """Admit a saved cleaning; never silently fix its geometry or identity.

    Returns the raw city and the cleaned city. The cleaned city must carry a
    cleaning record that is valid for its buildings, so meshing reuses it.
    """
    path = Path(path)
    try:
        payload = json.loads(path.read_text())
        metadata = payload["metadata"]
        if (
            payload["type"] != "FeatureCollection"
            or metadata["version"] != ARTIFACT_VERSION
            or metadata["crs"] != "EPSG:3006"
        ):
            raise ValueError("unsupported cleaning artifact format or CRS")
        if metadata["bounds"] != list(task["case"]["bounds"]):
            raise ValueError("cleaning artifact bounds do not match the selected case")
        parameters = {k: task["parameters"][k] for k in CLEANING_PARAMETERS}
        if metadata["parameters"] != parameters:
            raise ValueError(
                "cleaning artifact parameters do not match the selected case"
            )
        raw_city = load_model(path.parent / "raw.dtcc", expected_type=City)
        cleaned_city = load_model(path.parent / "cleaned.dtcc", expected_type=City)
        record = cleaning_record(
            cleaned_city, min_building_detail=parameters["min_building_detail"]
        )
        if record is None:
            raise ValueError(
                "cleaned city has no valid cleaning record for these parameters"
            )
        if metadata.get("cleaning") != record:
            raise ValueError("cleaning metadata does not match the cleaned city")
        buildings = cleaned_city.buildings
        for feature in payload["features"]:
            properties = feature["properties"]
            index = properties["building_index"]
            if type(index) is not int or not 0 <= index < len(buildings):
                raise ValueError("invalid building index")
            if properties["building_id"] != buildings[index].id or not shape(
                feature["geometry"]
            ).equals_exact(_footprint_polygon(buildings[index]), 0.0):
                raise ValueError("saved footprints do not match the cleaned city")
        return raw_city, cleaned_city
    except (KeyError, TypeError, OSError, ValueError) as exc:
        raise ValueError(f"Invalid cleaning artifact {path}: {exc}") from exc


def _cleaning_stage(cleaned_city, parameters):
    """The footprint stage the mesh builders will use, from the record."""
    return meshes._meshing_footprints(
        cleaned_city,
        lod=None,
        min_building_detail=parameters["min_building_detail"],
        min_building_area=parameters["min_building_area"],
        max_mesh_size=parameters["max_mesh_size"],
    )


def _flat_city(city, bounds):
    # Flat meshing only needs a rectangular domain, not LiDAR or roof heights.
    raster = Raster()
    raster.data = np.zeros((2, 2))
    raster.georef = Affine(
        (bounds.xmax - bounds.xmin) / 2,
        0,
        bounds.xmin,
        0,
        -(bounds.ymax - bounds.ymin) / 2,
        bounds.ymax,
    )
    city.add_terrain(raster)
    return city


def execute_phases(task, metrics, artifacts):
    """Mutate phase results as they complete so later failures retain evidence."""
    from benchmarks.benchmark_datasets import (
        _mesh_metrics,
        _metric_contract_warnings,
        _split_contract_warnings,
    )

    directory = Path(task["task_dir"])
    parameters = task["parameters"]
    bounds = Bounds(*task["case"]["bounds"])
    phase = task.get("phase", "both")
    dataset = task["dataset"]
    cleaning_dir = directory / "cleaning"
    started = time.perf_counter()
    metrics["active_phase"] = "input"
    if phase == "meshing":
        source_path = Path(task["cleaning_input"])
        raw_city, cleaned_city = load_cleaning(source_path, task)
        metrics["cleaning"] = {"status": "reused", "input": str(source_path)}
        # Make the new result self-contained for another replay.
        path = save_cleaning(
            cleaning_dir,
            raw_city,
            cleaned_city,
            _cleaning_stage(cleaned_city, parameters),
            task,
        )
    else:
        raw_city = prepare_footprint_city_from_bounds(bounds)
        metrics["input"] = {
            "status": "success",
            "seconds": time.perf_counter() - started,
            "source": "provider/cache",
        }
        raw = [
            polygon
            for b in raw_city.buildings
            if (geometry := b.flatten_geometry(GeometryType.LOD0)) is not None
            if (polygon := geometry.to_polygon(simplify=0.0)) is not None
            and not polygon.is_empty
        ]
        metrics["active_phase"] = "cleaning"
        started = time.perf_counter()
        cleaned_city = raw_city.clean_footprints(
            min_building_detail=parameters["min_building_detail"],
            lod=GeometryType.LOD0,
        )
        stage = _cleaning_stage(cleaned_city, parameters)
        seconds = time.perf_counter() - started
        cleaned = [s.to_polygon(simplify=0.0) for s in stage.surfaces]
        cleaning_warnings, _ = _split_contract_warnings(
            _metric_contract_warnings({"contract": stage.contract})
        )
        metrics["cleaning"] = {
            "status": "warning" if cleaning_warnings else "success",
            "seconds": seconds,
            **footprint_metrics(
                raw,
                cleaned,
                _raw_source_map(
                    raw_city, cleaned_city, [i for (i,) in stage.source_map]
                ),
                stage.declared_scale,
                epsilon=stage.diagnostics.get("fidelity_budget_used"),
            ),
            "contract": stage.contract,
            "before_selection_contract": stage.diagnostics.get(
                "before_selection_contract"
            ),
            "selection": stage.diagnostics.get("selection"),
        }
        path = save_cleaning(cleaning_dir, raw_city, cleaned_city, stage, task)
    artifacts["cleaning_input"] = {"format": "geojson", "path": str(path)}
    write_json(cleaning_dir / "metrics.json", metrics["cleaning"])
    if phase == "cleaning" or dataset == "city_footprints":
        metrics.pop("active_phase", None)
        return None

    metrics["active_phase"] = "mesh_input"
    started = time.perf_counter()
    mesh_dir = directory / "meshing"
    mesh_dir.mkdir(parents=True, exist_ok=True)
    preparation = {k: parameters[k] for k in PREPARATION_PARAMETERS}
    saved_city = task.get("prepared_city")
    if dataset == "city_flat_mesh":
        city = _flat_city(cleaned_city, bounds)
        source = "bounds"
    elif saved_city:
        prepared = load_model(saved_city, expected_type=City)
        if len(prepared.buildings) != len(cleaned_city.buildings) or any(
            original.id != enriched.id
            or not _footprint_polygon(original).equals_exact(
                _footprint_polygon(enriched), 0.0
            )
            for original, enriched in zip(cleaned_city.buildings, prepared.buildings)
        ):
            raise ValueError(
                "Saved mesh input does not match the cleaned buildings"
            )
        city = prepared
        source = "saved"
    else:
        city = prepare_city_from_bounds(
            bounds, buildings=cleaned_city.buildings, **preparation
        )
        # Heights change only z, so the cleaning record stays valid.
        city.attributes[RECORD_KEY] = cleaned_city.attributes[RECORD_KEY]
        source = "provider/cache"
    city_path = mesh_dir / "city.dtcc"
    save_model(city, city_path)
    artifacts["prepared_city"] = {"format": "dtcc", "path": str(city_path)}
    metrics["mesh_input"] = {
        "status": "success",
        "seconds": time.perf_counter() - started,
        "source": source,
        "parameters": preparation,
    }
    metrics["active_phase"] = "meshing"
    started = time.perf_counter()
    build = getattr(meshes, "build_" + dataset)
    kwargs = {
        k: parameters[k]
        for k in ("max_mesh_size", "min_mesh_angle", *CLEANING_PARAMETERS)
    }
    if dataset == "city_volume_mesh":
        kwargs["domain_height"] = parameters["domain_height"]
    mesh = build(
        city,
        report_mesh_quality=False,
        stage_audit={},
        **kwargs,
    )
    seconds = time.perf_counter() - started
    metrics["active_phase"] = "meshing"
    mesh_metrics = _mesh_metrics(mesh)
    metrics["meshing"] = {"status": "success", "seconds": seconds, **mesh_metrics}
    # Preserve the existing report's mesh columns and contract-warning handling.
    metrics.update(mesh_metrics)
    write_json(mesh_dir / "metrics.json", metrics["meshing"])
    metrics.pop("active_phase", None)
    return mesh
