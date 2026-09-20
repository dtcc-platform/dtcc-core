"""Cleaning handoff and phase measurements for the city benchmark.

Saved coordinates are in EPSG:3006 (metres), never GeoJSON's implied WGS84.
The mesh builders remain the authority for consuming conditioned footprints.
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

from dtcc_core.model import Bounds, City, GeometryType, Raster, Surface
from dtcc_core.io.model import load_model, save_model
from dtcc_core.builder.geometry_builders import meshes
from dtcc_core.datasets._city_mesh_common import (
    prepare_city_from_bounds,
    prepare_footprint_city_from_bounds,
)
from benchmarks.benchmark_catalog import CLEANING_PARAMETERS

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
    lengths = [
        np.linalg.norm(np.diff(np.asarray(ring.coords)[:, :2], axis=0), axis=1)
        for polygon in cleaned
        for ring in [polygon.exterior, *polygon.interiors]
    ]
    lengths = np.concatenate(lengths) if lengths else np.array([])
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
        "min_clearance": min((p.minimum_clearance for p in cleaned), default=None),
        "min_edge_length": float(lengths.min()) if lengths.size else None,
        "short_edge_count": int(np.count_nonzero(lengths < scale - 1e-9)),
        "vertex_count": int(lengths.size),
    }


def save_cleaning(directory, city, conditioned, task):
    directory.mkdir(parents=True, exist_ok=True)
    save_model(city, directory / "raw.dtcc")
    payload = {
        "type": "FeatureCollection",
        "metadata": {
            "version": 1,
            "crs": "EPSG:3006",
            "bounds": task["case"]["bounds"],
            "parameters": {k: task["parameters"][k] for k in CLEANING_PARAMETERS},
            "declared_scale": conditioned.declared_scale,
            "before_selection_contract": conditioned.diagnostics.get(
                "before_selection_contract"
            ),
            "selection": conditioned.diagnostics.get("selection"),
            "final_handoff_contract": conditioned.diagnostics.get(
                "final_handoff_contract"
            ),
            "policy_exclusions": conditioned.diagnostics.get(
                "policy_exclusions", []
            ),
            "diagnostics": {
                k: conditioned.diagnostics[k]
                for k in (
                    "precision_grid",
                    "output_grid",
                    "mesher_ready_coverage_revalidation_output_grid",
                    "mesher_ready_coverage_revalidation_applied",
                    "geos_exception_count",
                )
                if k in conditioned.diagnostics
            },
        },
        "features": [
            {
                "type": "Feature",
                "geometry": mapping(surface.to_polygon(simplify=0.0)),
                "properties": {"source_indices": indices},
            }
            for surface, indices in zip(conditioned.surfaces, conditioned.source_map)
        ],
    }
    path = directory / "footprints.geojson"
    write_json(path, payload)
    return path


def load_cleaning(path, task):
    """Admit a persisted handoff; never silently fix its geometry or identity."""
    path = Path(path)
    try:
        payload = json.loads(path.read_text())
        metadata = payload["metadata"]
        if (
            payload["type"] != "FeatureCollection"
            or metadata["version"] != 1
            or metadata["crs"] != "EPSG:3006"
        ):
            raise ValueError("unsupported cleaning artifact format or CRS")
        if metadata["bounds"] != list(task["case"]["bounds"]):
            raise ValueError("cleaning artifact bounds do not match the selected case")
        if metadata["parameters"] != {
            k: task["parameters"][k] for k in CLEANING_PARAMETERS
        }:
            raise ValueError(
                "cleaning artifact parameters do not match the selected case"
            )
        scale = float(metadata["declared_scale"])
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError("invalid declared scale")
        diagnostics = metadata["diagnostics"]
        for key, value in diagnostics.items():
            if key == "mesher_ready_coverage_revalidation_applied":
                if type(value) is not bool:
                    raise ValueError("invalid revalidation flag")
            elif type(value) not in (float, int) or not np.isfinite(value) or value < 0:
                raise ValueError("invalid cleaning diagnostic")
        before_selection_contract = metadata.get("before_selection_contract")
        selection = metadata.get("selection")
        final_handoff_contract = metadata.get("final_handoff_contract")
        policy_exclusions = metadata.get("policy_exclusions", [])
        if not isinstance(before_selection_contract, dict):
            raise ValueError("missing before-selection contract")
        if not isinstance(selection, dict):
            raise ValueError("missing selection metadata")
        if final_handoff_contract is not None and not isinstance(
            final_handoff_contract, dict
        ):
            raise ValueError("invalid final handoff contract")
        if not isinstance(policy_exclusions, list):
            raise ValueError("invalid policy exclusions")
        city = load_model(path.parent / "raw.dtcc", expected_type=City)
        surfaces, source_map = [], []
        for feature in payload["features"]:
            polygon = shape(feature["geometry"])
            indices = feature["properties"]["source_indices"]
            if (
                polygon.geom_type != "Polygon"
                or polygon.is_empty
                or not polygon.is_valid
            ):
                raise ValueError("invalid cleaned polygon")
            if not all(
                np.isfinite(np.asarray(r.coords)).all()
                for r in [polygon.exterior, *polygon.interiors]
            ):
                raise ValueError("nonfinite cleaned coordinates")
            if (
                not isinstance(indices, list)
                or not indices
                or any(
                    type(i) is not int or not 0 <= i < len(city.buildings)
                    for i in indices
                )
                or indices != sorted(set(indices))
            ):
                raise ValueError("invalid source indices")
            surface = Surface()
            surface.from_polygon(polygon)
            surfaces.append(surface)
            source_map.append(indices)
        required_selection_counts = (
            "input_count",
            "output_count",
            "removed_count",
        )
        if any(
            type(selection.get(key)) is not int or selection[key] < 0
            for key in required_selection_counts
        ):
            raise ValueError("invalid selection counts")
        if (
            selection["output_count"] != len(surfaces)
            or selection["input_count"]
            != selection["output_count"] + selection["removed_count"]
        ):
            raise ValueError("selection counts do not match saved geometry")
        for field in ("excluded_regions",):
            if not isinstance(selection.get(field), list):
                raise ValueError(f"invalid selection {field}")
        for exclusion in [*selection["excluded_regions"], *policy_exclusions]:
            if not isinstance(exclusion, dict):
                raise ValueError("invalid selection or policy exclusion")
            indices = exclusion.get("source_indices")
            area = exclusion.get("area")
            reason = exclusion.get("reason")
            if (
                not isinstance(indices, list)
                or indices != sorted(set(indices))
                or any(
                    type(index) is not int or not 0 <= index < len(city.buildings)
                    for index in indices
                )
                or type(area) not in (float, int)
                or not np.isfinite(area)
                or area < 0
                or not isinstance(reason, str)
                or not reason
            ):
                raise ValueError("invalid selection or policy exclusion")
        for field in (
            "represented_source_indices",
            "unrepresented_source_indices",
        ):
            indices = selection.get(field)
            if (
                not isinstance(indices, list)
                or indices != sorted(set(indices))
                or any(
                    type(index) is not int or not 0 <= index < len(city.buildings)
                    for index in indices
                )
            ):
                raise ValueError(f"invalid selection {field}")
        represented = sorted({index for indices in source_map for index in indices})
        if selection["represented_source_indices"] != represented:
            raise ValueError("selection represented sources do not match geometry")
        diagnostics = {
            **diagnostics,
            "before_selection_contract": before_selection_contract,
            "selection": selection,
            "policy_exclusions": policy_exclusions,
        }
        if final_handoff_contract is not None:
            diagnostics["final_handoff_contract"] = final_handoff_contract
        conditioned = meshes.ConditionedFootprints(
            surfaces,
            source_map,
            [],
            diagnostics,
            scale,
            {},
        )
        # Graph and scale validation happens at the mesh-builder input boundary.
        return city, conditioned
    except (KeyError, TypeError, OSError, ValueError) as exc:
        raise ValueError(f"Invalid cleaning artifact {path}: {exc}") from exc


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
        city, conditioned = load_cleaning(source_path, task)
        metrics["cleaning"] = {"status": "reused", "input": str(source_path)}
        # Make the new result self-contained for another replay.
        path = save_cleaning(cleaning_dir, city, conditioned, task)
    else:
        city = prepare_footprint_city_from_bounds(bounds)
        metrics["input"] = {
            "status": "success",
            "seconds": time.perf_counter() - started,
            "source": "provider/cache",
        }
        raw = [
            polygon
            for b in city.buildings
            if (geometry := b.flatten_geometry(GeometryType.LOD0)) is not None
            if (polygon := geometry.to_polygon(simplify=0.0)) is not None
            and not polygon.is_empty
        ]
        metrics["active_phase"] = "cleaning"
        started = time.perf_counter()
        conditioned = meshes.build_conditioned_footprints(
            city.buildings,
            lod=GeometryType.LOD0,
            **{k: parameters[k] for k in CLEANING_PARAMETERS},
            max_mesh_size=parameters["max_mesh_size"],
            cleaning_diagnostics=False,
        )
        seconds = time.perf_counter() - started
        cleaned = [s.to_polygon(simplify=0.0) for s in conditioned.surfaces]
        cleaning_warnings, _ = _split_contract_warnings(
            _metric_contract_warnings({"contract": conditioned.contract})
        )
        metrics["cleaning"] = {
            "status": "warning" if cleaning_warnings else "success",
            "seconds": seconds,
            **footprint_metrics(
                raw, cleaned, conditioned.source_map, conditioned.declared_scale
            ),
            "contract": conditioned.contract,
            "before_selection_contract": conditioned.diagnostics.get(
                "before_selection_contract"
            ),
            "selection": conditioned.diagnostics.get("selection"),
        }
        path = save_cleaning(cleaning_dir, city, conditioned, task)
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
        city = _flat_city(city, bounds)
        source = "bounds"
    elif saved_city:
        prepared = load_model(saved_city, expected_type=City)
        if len(prepared.buildings) != len(city.buildings) or any(
            original.id != enriched.id
            or not original.footprint()
            .to_polygon(simplify=0.0)
            .equals_exact(enriched.footprint().to_polygon(simplify=0.0), 0.0)
            for original, enriched in zip(city.buildings, prepared.buildings)
        ):
            raise ValueError(
                "Saved mesh input does not match the cleaning source buildings"
            )
        city = prepared
        source = "saved"
    else:
        city = prepare_city_from_bounds(bounds, buildings=city.buildings, **preparation)
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
        lod=GeometryType.LOD0,
        conditioned_footprints=conditioned,
        report_mesh_quality=False,
        cleaning_diagnostics=False,
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
