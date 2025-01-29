from pathlib import Path
import numpy as np

from dtcc_core.model import PointCloud, Bounds
from .container import PointCloudDirectory
from . import generic
from .logging import info, warning, error
from typing import List, Union

from .pointcloud import las_file_bounds



def load(las_directory: Union[str, Path]) -> PointCloudDirectory:
    """
    create a PointCloudContainer object from a directory of las files
    """
    las_directory = Path(las_directory)
    if not las_directory.exists():
        raise ValueError(f"Path {las_directory} does not exist")
    if not las_directory.is_dir():
        raise ValueError(f"Path {las_directory} is not a directory")

    las_files = list(las_directory.glob("*.la[sz]"))
    if len(las_files) == 0:
        warning(f"No valid LAS files found in {las_directory}, returning empty PointCloudContainer")
        return PointCloudDirectory()
    info(f"Found {len(las_files)} LAS files in {las_directory}")
    files = []
    bounds = []
    for lf in las_files:
        try:
            bounds.append(las_file_bounds(lf))
            files.append(lf)
        except Exception as e:
            error(f"Failed to load {lf}: {e}")
    return PointCloudDirectory(files, bounds)


