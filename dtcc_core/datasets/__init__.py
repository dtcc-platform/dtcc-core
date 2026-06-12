from __future__ import annotations

from .dataset import DatasetDescriptor, DatasetBaseArgs, DatasetExportResult
from .publish import (
    DatasetPackageError,
    DatasetPublication,
    DatasetPublishConfigurationError,
    DatasetPublishError,
    DatasetUploadClient,
    DatasetUploadConflictError,
    DatasetUploadError,
    DatasetUploadInProgressError,
    DatasetUploadRateLimitError,
    PublishedFile,
)
from dtcc_core.common import info as log_info, log_table

# Import dataset classes to trigger auto-registration
from .pointcloud import PointCloudDataset
from .buildings import BuildingDataset
from .city import CityDataset
from .terrain_surface_mesh import TerrainSurfaceMeshDataset
from .city_surface_mesh import CitySurfaceMeshDataset
from .city_flat_mesh import CityFlatMeshDataset
from .city_footprints import CityFootprintsDataset
from .city_volume_mesh import CityVolumeMeshDataset
from .air_quality import AirQualityDataset
from .trees import TreesDataset
from .weather import WeatherDataset
from .hydrology import HydrologyDataset
from .ocean import OceanDataset
from .footprints import FootprintsDataset
from .roads import RoadsDataset
from .space_syntax import SpaceSyntaxDataset
from .transit_vehicles import (
    TransitVehiclesDataset,
    BusesDataset,
    TramsDataset,
    TrainsDataset,
    MetrosDataset,
    FerriesDataset,
)
from .deso import DeSODataset
from .smoke import SmokeDataset
from .calibration_grid import CalibrationGridDataset

# Remote dataset support (optional, requires httpx)
try:
    from .remote import (
        RemoteDatasetDescriptor,
        RemoteValidationError,
        register_remote_service,
        register_remote_descriptors_from_cache,
        get_cached_discoveries,
    )
except ImportError:
    pass

# Import registry infrastructure from separate module to avoid circular imports
from .registry import (
    _register_dataset_class,
    register,
    register_class,
    unregister,
    list_datasets as list,
    get_dataset,
)

# Override module attributes with registered dataset instances
# This is necessary because importing from submodules makes the submodules
# themselves available as attributes, which would shadow the registered instances
point_cloud = get_dataset("point_cloud")
buildings = get_dataset("buildings")
building_footprints = get_dataset("building_footprints")
city = get_dataset("city")
city_footprints = get_dataset("city_footprints")

terrain_surface_mesh = get_dataset("terrain_surface_mesh")
city_surface_mesh = get_dataset("city_surface_mesh")
city_flat_mesh = get_dataset("city_flat_mesh")
city_volume_mesh = get_dataset("city_volume_mesh")
air_quality = get_dataset("air_quality")
trees = get_dataset("trees")
weather = get_dataset("weather")
hydrology = get_dataset("hydrology")
ocean = get_dataset("ocean")
roads = get_dataset("roads")
space_syntax = get_dataset("space_syntax")
transit_vehicles = get_dataset("transit_vehicles")
buses = get_dataset("buses")
trams = get_dataset("trams")
trains = get_dataset("trains")
metros = get_dataset("metros")
ferries = get_dataset("ferries")
deso = get_dataset("deso")
smoke = get_dataset("smoke")
calibration_grid = get_dataset("calibration_grid")


_CATEGORY_ORDER = ("raw", "derived", "simulation", "remote", "unknown")


def info():
    """
    Print information about all available datasets.

    This function iterates over all registered datasets and prints
    a nicely formatted summary of each one.

    Example:
        >>> import dtcc_core.datasets as datasets
        >>> datasets.info()
    """
    datasets_dict = list()

    if not datasets_dict:
        log_info("DTCC Datasets: no datasets are currently registered.")
        return

    grouped = _group_dataset_summary_rows(datasets_dict)
    log_info(f"DTCC Datasets ({len(datasets_dict)} available)")
    log_info("Use datasets.<name>() to access a dataset.")
    log_info("Use print(datasets.<name>) to see dataset parameters.")

    columns = [
        ("Dataset", "left"),
        ("Product", "left"),
        ("Formats", "left"),
        ("Source", "left"),
    ]
    for category, rows in grouped.items():
        log_table(
            log_info,
            f"{_format_category_title(category)} ({len(rows)})",
            columns,
            rows,
        )


def __getattr__(name):
    """
    Dynamic attribute lookup for registered datasets.

    This enables access like `datasets.pointcloud`.

    Args:
        name (str): The attribute name being accessed.

    Returns:
        The registered dataset instance, if found.

    Raises:
        AttributeError: If the attribute is not a registered dataset.
    """
    try:
        return get_dataset(name)
    except KeyError:
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


def _group_dataset_summary_rows(
    datasets_dict: dict[str, DatasetDescriptor],
) -> dict[str, list[tuple[str, str, str, str]]]:
    grouped: dict[str, list[tuple[str, str, str, str]]] = {}
    for name, dataset in sorted(datasets_dict.items()):
        meta = dataset.describe()
        category = str(meta.get("data_category") or "unknown")
        grouped.setdefault(category, []).append(
            (
                name,
                _format_identifier(str(meta.get("result_kind") or "unknown")),
                _format_supported_formats(meta.get("supported_formats")),
                _dataset_source(dataset),
            )
        )
    return {
        category: grouped[category]
        for category in _ordered_categories(grouped)
    }


def _ordered_categories(grouped: dict[str, list]) -> list[str]:
    known = [category for category in _CATEGORY_ORDER if category in grouped]
    extra = sorted(category for category in grouped if category not in _CATEGORY_ORDER)
    return known + extra


def _format_supported_formats(formats) -> str:
    if not formats:
        return "python object"
    return ", ".join(str(fmt) for fmt in formats)


def _dataset_source(dataset: DatasetDescriptor) -> str:
    return str(getattr(dataset, "source_service", None) or "dtcc-core")


def _format_category_title(category: str) -> str:
    return f"{_format_identifier(category).title()} Datasets"


def _format_identifier(value: str) -> str:
    return value.replace("_", " ")
