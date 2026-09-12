"""Build native DTCC data and project only its semantic facts for the experiment."""

import argparse
import json
from pathlib import Path

import numpy as np

from dtcc_core import io
from dtcc_core.datasets import load_model_package
from dtcc_core.datasets.schema import (
    DatasetContext, DatasetIdentity, DatasetMetadata, DatasetPresentation,
    DatasetProvenance, DatasetRequest,
)
from dtcc_core.model import Building, BuildingPart, City, Field, Mesh, Object, Point


NAMESPACE = "https://example.org/dtcc/"


def project(root):
    """Project native semantic facts and derive parent IDs from actual containment.

    Numerical arrays stay in DTCC. This projection is only a profile-validator input.
    """
    records = []
    seen = set()
    seen_ids = set()
    stack = [(root, None)]
    while stack:
        obj, parent = stack.pop()
        if id(obj) in seen:
            raise ValueError(f"Repeated object in containment at {obj.id!r}")
        seen.add(id(obj))
        if not obj.id or obj.id in seen_ids:
            raise ValueError(f"Missing or duplicate model ID: {obj.id!r}")
        seen_ids.add(obj.id)
        if obj.semantic_type is None:
            raise ValueError(f"Missing explicit semantic type for {obj.id!r}")
        references = {name: list(targets) for name, targets in obj.relations.items()}
        if "parent" in references:
            raise ValueError("Containment must come from Object.children")
        if parent is not None:
            references["parent"] = parent
        records.append({
            "id": obj.id,
            "semantic_type": obj.semantic_type,
            "attributes": dict(obj.attributes),
            "relations": references,
        })
        children = [child for group in obj.children.values() for child in group]
        stack.extend((child, obj.id) for child in reversed(children))
    return records


def build_example():
    city = City(id="city-1", semantic_type=NAMESPACE + "City",
                profile_id=NAMESPACE + "profiles/city", profile_version="0.1.0", attributes={"name": "Profile experiment"})
    building = Building(id="building-1", semantic_type=NAMESPACE + "Building", attributes={"height": 12.5, "usage": "residential"})
    part = BuildingPart(id="part-1", semantic_type=NAMESPACE + "BuildingPart")
    sensor = Object(id="sensor-1", semantic_type=NAMESPACE + "Sensor",
                    relations={"observes": ["building-1"]}, attributes={"phenomenon": "temperature"})
    city.add_child(building)
    building.add_child(part)
    city.add_child(sensor)
    mesh = Mesh(
        vertices=np.array([
            [325000.123, 6400000.235, 12.125],
            [325002.123, 6400000.235, 12.125],
            [325000.123, 6400002.235, 12.125],
        ], dtype=np.float64),
        faces=np.array([[0, 1, 2]], dtype=np.int64),
        markers=np.array([-2], dtype=np.int32),
        normals=np.array([[0., 0., 1.]], dtype=np.float64),
    )
    mesh.transform.srs = "EPSG:3006"
    field = Field(name="temperature", unit="K", association="vertex", values=np.array([293.15, 294.25, 295.35]))
    mesh.add_field(field)
    part.add_mesh(mesh)
    sensor.add_geometry(Point(x=325001.123, y=6400001.235, z=14), "location")
    return city, mesh


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    city, mesh = build_example()
    objects = project(city)
    file_path = args.output.with_suffix(".dtcc")
    city.save(file_path)
    restored = io.load_city(file_path)
    assert project(restored) == objects
    restored_mesh = restored.buildings[0].building_parts[0].mesh
    for name in ("vertices", "faces", "markers", "normals"):
        original, result = getattr(mesh, name), getattr(restored_mesh, name)
        np.testing.assert_array_equal(original, result)
        assert original.dtype == result.dtype and original.shape == result.shape
    np.testing.assert_array_equal(mesh.fields[0].values, restored_mesh.fields[0].values)
    assert restored_mesh.fields[0].association == "vertex"
    assert restored_mesh.fields[0].unit == "K"
    assert mesh.transform.srs == restored_mesh.transform.srs
    np.testing.assert_array_equal(mesh.transform.affine, restored_mesh.transform.affine)
    city.dataset_context = DatasetContext(
        identity=DatasetIdentity(name="profile-example", title="Profile example"),
        metadata=DatasetMetadata(crs=["EPSG:3006"]),
        provenance=DatasetProvenance(sources=["synthetic"]),
        presentation=DatasetPresentation(),
        request=DatasetRequest(dataset_name="profile-example"),
        health={"status": "ok"}, warnings=["Synthetic demonstration"],
    )
    package = city.export(args.output.with_suffix(".dtccpkg"), canonical=True)
    packaged = load_model_package(package.path)
    assert project(packaged) == objects
    assert packaged.dataset_context == city.dataset_context
    extension = restored.copy()
    extension.profile_id = NAMESPACE + "profiles/low-rise"
    extension.buildings[0].attributes["height"] = 8.0
    extension.children[Object][0].attributes["accuracy"] = 0.1
    extension.children[Object][0].semantic_type = NAMESPACE + "WeatherStation"
    assert type(extension.children[Object][0]) is Object
    payload = {
        "objects": project(restored),
        "extension_objects": project(extension),
        "numerical_example": {
            "association": restored_mesh.fields[0].association,
            "vertices_dtype": str(restored_mesh.vertices.dtype),
            "vertices_shape": list(restored_mesh.vertices.shape),
            "field_shape": list(restored_mesh.fields[0].values.shape),
            "crs": restored_mesh.transform.srs,
        },
        "canonical_exchange_observation": {
            "max_coordinate_error_m": float(np.max(np.abs(mesh.vertices - restored_mesh.vertices))),
            "markers_preserved": True,
            "normals_preserved": True,
            "field_shape_preserved": mesh.fields[0].values.shape == restored_mesh.fields[0].values.shape,
            "max_field_value_error": float(np.max(np.abs(
                mesh.fields[0].values - restored_mesh.fields[0].values
            ))),
            "semantic_facts_preserved": project(restored) == objects,
            "package_context_preserved": packaged.dataset_context == city.dataset_context,
        },
    }
    args.output.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(f"Wrote canonical {file_path} and {package.path}; projected {len(objects)} restored objects.")
    print(json.dumps(payload["canonical_exchange_observation"], indent=2))


if __name__ == "__main__":
    main()
