import pytest
from pathlib import Path
from dtcc_core import builder, io


@pytest.fixture
def data_dir():
    return (Path(__file__).parent / "../data").resolve()


@pytest.fixture
def minimal_case_footprints(data_dir):
    return data_dir / "MinimalCase" / "PropertyMap.shp"


def test_convert_surface(minimal_case_footprints):
    footprints = io.load_footprints(minimal_case_footprints)
    footprint = footprints[0].get_footprint()
    builder_surface = builder.model_conversion.create_builder_surface(footprint)
    assert len(builder_surface.vertices) == 4
    assert footprint.vertices[0][0] == builder_surface.vertices[0].x


if __name__ == "__main__":
    pytest.main()
