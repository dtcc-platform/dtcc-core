import os
from functools import partial
from pathlib import Path
from typing import List, Union

import numpy as np
import rasterio
import rasterio.merge
from PIL import Image
from rasterio.transform import from_origin

from ..model import Raster
from . import generic
from .logging import error, info, warning
from .model import load_model, save_model


def _load_rasterio(path: Path | list, **kwargs):
    raster = Raster()
    if not isinstance(path, list):
        path = Path(path)
        if not path.is_file():
            error(f"File {path} does not exist or is not a file.")
            return raster
        with rasterio.open(path) as src:
            data = src.read()
            data = data.squeeze()
            if data.ndim == 3:
                # rasterio returns (channels, height, width)
                # we want (width, heigh, channels)
                data = np.moveaxis(data, 0, -1)

            raster.data = np.squeeze(data)
            raster.georef = src.transform
            raster.crs = str(src.crs)
    else:
        if "merge_method" in kwargs:
            merge_method = kwargs["merge_method"]
            if merge_method not in ["first", "last", "min", "max"]:
                warning(f"Invalid merge method: {merge_method}. Using 'first' instead.")
                merge_method = "first"
        else:
            merge_method = "first"
        data, transform = rasterio.merge.merge(path, method=merge_method)
        raster.data = data.squeeze()
        raster.georef = transform
        with rasterio.open(path[0]) as src:
            raster.crs = str(src.crs)

    return raster


def _load_csv(path, delimiter=",", **kwargs):
    raster = Raster()
    data = np.loadtxt(path, delimiter=delimiter)
    raster.data = data
    return raster


def load(path, delimiter=",", *, validate_schema=True) -> Raster:
    """
    Load a raster file as a `Raster` object.

    Args:
        path (str): The path to the raster file.
        delimiter (str): The delimiter used in case of a CSV file (default ",").

    Returns:
        Raster: A `Raster` object representing the raster file loaded.
    """

    if isinstance(path, (str, Path)) and Path(path).suffix.lower() == ".dtcc":
        return load_model(path, expected_type=Raster, validate_schema=validate_schema)
    if validate_schema is not True:
        raise ValueError("validate_schema applies to .dtcc input")
    return generic.load(path, "raster", Raster, _load_formats, delimiter=delimiter)


def _save_json_raster(raster, path):
    path.write_text(raster.to_json())


def _save_geotif(raster, path):
    data = raster.data
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=raster.height,
        width=raster.width,
        count=raster.channels,
        dtype=data.dtype,
        transform=raster.georef,
        compress="DEFLATE",
    ) as dst:
        if raster.channels == 1:
            dst.write(data, 1)
        else:
            dst.write(data, list(range(1, raster.channels + 1)))
    return True


def _save_image(raster, path):
    suffix = path.suffix.lower()
    wld_suffix = f"{suffix[0]}{suffix[-1]}w"
    wld_path = path.with_suffix(wld_suffix)
    data = raster.data
    if raster.channels == 1:
        data = np.repeat(data[:, :, np.newaxis], 3, axis=2)
    with open(wld_path, "w") as f:
        f.write(
            f"{raster.georef.a}\n{raster.georef.b}\n{raster.georef.d}\n{raster.georef.e}\n{raster.georef.xoff}\n{raster.georef.yoff}"
        )
    im = Image.fromarray(data)
    im.save(path)
    return True


def _save_csv(raster, path):
    np.savetxt(path, raster.data, delimiter=",")
    return True


def save(raster: Raster, path, **kwargs):
    """
    Save a `Raster` object to a file.

    Parameters
    ----------
    raster : Raster
        The `Raster` object to save.
    path : str
        The path to the output file.
    """

    path = Path(path)
    return generic.save(raster, path, "raster", _save_formats, **kwargs)


_load_formats = {
    Raster: {
        ".dtcc": partial(load_model, expected_type=Raster),
        ".tif": _load_rasterio,
        ".geotif": _load_rasterio,
        ".png": _load_rasterio,
        ".asc": _load_rasterio,
        ".csv": _load_csv,
        ".txt": _load_csv,
    }
}

_save_formats = {
    Raster: {
        ".dtcc": save_model,
        ".json": _save_json_raster,
        ".tif": _save_geotif,
        ".png": _save_image,
        ".jpg": _save_image,
        ".csv": _save_csv,
    }
}
