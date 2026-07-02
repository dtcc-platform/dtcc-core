import xml.etree.ElementTree as ET

import h5py
import numpy as np
import pytest

from dtcc_core import io
from dtcc_core.model import Field, VolumeMesh


def _single_tet_mesh():
    return VolumeMesh(
        vertices=np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            dtype=float,
        ),
        cells=np.array([[0, 1, 2, 3]], dtype=np.int64),
    )


def _xdmf_attributes(path):
    tree = ET.parse(path)
    return {
        attribute.attrib["Name"]: attribute.attrib
        for attribute in tree.getroot().iter("Attribute")
    }


def _field_dataset_by_name(h5_file, field_name):
    fields_group = h5_file["Mesh/mesh/fields"]
    for dataset in fields_group.values():
        if dataset.attrs.get("name") == field_name:
            return dataset
    raise AssertionError(f"HDF5 field dataset not found for {field_name!r}")


def test_save_xdmf_volume_mesh_writes_node_cell_and_vector_fields(tmp_path):
    mesh = _single_tet_mesh()
    temperature = np.array([18.0, 19.0, 20.0, 21.0])
    cell_quality = np.array([0.75])
    velocity = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ]
    )
    pollutant = np.array([30.0, 31.0, 32.0, 33.0])

    mesh.add_field(Field(name="temperature", unit="degC", values=temperature, dim=1))
    mesh.add_field(Field(name="cell_quality", values=cell_quality, dim=1))
    mesh.add_field(Field(name="velocity", unit="m/s", values=velocity, dim=3))
    mesh.add_field(Field(name="NO2/urban mean", values=pollutant, dim=1))

    path = tmp_path / "field_mesh.xdmf"
    mesh.save(path)

    h5_path = tmp_path / "field_mesh.h5"
    assert path.exists()
    assert h5_path.exists()

    attributes = _xdmf_attributes(path)
    assert attributes["temperature"]["Center"] == "Node"
    assert attributes["temperature"]["AttributeType"] == "Scalar"
    assert attributes["cell_quality"]["Center"] == "Cell"
    assert attributes["velocity"]["Center"] == "Node"
    assert attributes["velocity"]["AttributeType"] == "Vector"
    assert attributes["NO2/urban mean"]["Center"] == "Node"

    with h5py.File(h5_path, "r") as h5_file:
        np.testing.assert_allclose(
            _field_dataset_by_name(h5_file, "temperature")[()], temperature
        )
        np.testing.assert_allclose(
            _field_dataset_by_name(h5_file, "cell_quality")[()], cell_quality
        )
        np.testing.assert_allclose(
            _field_dataset_by_name(h5_file, "velocity")[()], velocity
        )
        pollutant_dataset = _field_dataset_by_name(h5_file, "NO2/urban mean")
        np.testing.assert_allclose(pollutant_dataset[()], pollutant)
        assert pollutant_dataset.attrs["name"] == "NO2/urban mean"


def test_save_xdmf_volume_mesh_rejects_invalid_field_shape(tmp_path):
    mesh = _single_tet_mesh()
    mesh.add_field(Field(name="temperature", values=np.array([1.0, 2.0]), dim=1))

    with pytest.raises(ValueError) as exc_info:
        mesh.save(tmp_path / "invalid_field.xdmf")

    message = str(exc_info.value)
    assert "temperature" in message
    assert "actual shape=(2,)" in message
    assert "expected vertex count=4" in message
    assert "expected cell count=1" in message


def test_load_xdmf_volume_mesh_preserves_native_fields(tmp_path):
    mesh = _single_tet_mesh()
    mesh.add_field(
        Field(name="temperature", values=np.array([18.0, 19.0, 20.0, 21.0]), dim=1)
    )
    mesh.add_field(
        Field(
            name="velocity",
            values=np.array(
                [
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [1.0, 1.0, 1.0],
                ]
            ),
            dim=3,
        )
    )

    path = tmp_path / "roundtrip.xdmf"
    mesh.save(path)

    loaded = io.load_volume_mesh(path)

    assert [field.name for field in loaded.fields] == ["temperature", "velocity"]
    np.testing.assert_allclose(
        loaded.fields[0].values.reshape(-1), [18.0, 19.0, 20.0, 21.0]
    )
    np.testing.assert_allclose(
        loaded.fields[1].values,
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ],
    )
