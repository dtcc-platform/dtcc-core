from functools import partial
from .model import load_model, save_model
import rasterio
import rasterio.merge
from rasterio.transform import from_origin
import os

from typing import Union, List
import numpy as np
from pathlib import Path
from PIL import Image

from . import generic

from ..model import Raster
from .logging import info, error, warning




def _load_rasterio(path: Union[Path, List], **kwargs):
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

    Parameters
    ----------
    path : str
        The path to the raster file.
    delimiter : str
        The delimiter used in case of a CSV file (default ",").

    Returns
    -------
    Raster
        A `Raster` object representing the raster file loaded.
    """




    if isinstance(path, (str, Path)) and Path(path).suffix.lower() == '.dtcc':
        return load_model(path, expected_type=Raster, validate_schema=validate_schema)
    if validate_schema is not True:
        raise ValueError('validate_schema applies to .dtcc input')
    return generic.load(path, "raster", Raster, _load_formats, delimiter=delimiter)




def _save_json_raster(raster, path):
    path.write_text(raster.to_json())


# Largest array copied at once when writing (height, width, channels) data.
WRITE_STRIP_BYTES = 8 * 1024**2


def _save_geotif(raster, path, alpha=False):
    """Write a GeoTIFF with the raster's CRS and nodata, if set.

    (height, width, channels) data is written band by band in windows of at
    most WRITE_STRIP_BYTES, so no band-first copy of the whole array is made.
    With ``alpha=True`` the last band is labelled alpha, which GeoTIFF supports
    for two bands (gray, alpha) or four (red, green, blue, alpha); otherwise no
    band is, whatever the band count. A CRS of "" or "None" (how the loader
    records a file without one) is not written.
    """
    if alpha and raster.channels not in (2, 4):
        raise ValueError(
            f"alpha needs a raster with 2 or 4 bands, not {raster.channels}"
        )
    data = raster.data
    profile = dict(
        driver="GTiff",
        height=raster.height,
        width=raster.width,
        count=raster.channels,
        dtype=data.dtype,
        transform=raster.georef,
        compress="DEFLATE",
    )
    if raster.crs and raster.crs != "None":
        profile["crs"] = raster.crs
    if not np.isnan(raster.nodata):
        profile["nodata"] = raster.nodata
    if raster.channels > 1:
        profile["alpha"] = "YES" if alpha else "UNSPECIFIED"
    with rasterio.open(path, "w", **profile) as dst:
        if raster.channels == 1:
            dst.write(data, 1)
        else:
            _write_bands_in_windows(dst, data)
    return True


def _write_bands_in_windows(dst, data):
    height, width, channels = data.shape
    itemsize = data.dtype.itemsize
    cols = max(1, min(width, WRITE_STRIP_BYTES // itemsize))
    rows = max(1, min(height, WRITE_STRIP_BYTES // (cols * itemsize)))
    # One buffer for every window; each window uses a contiguous prefix of it.
    # Rasterio writes a C-contiguous (1, rows, cols) array with a list of band
    # indexes without copying it; a 2-D array would be copied.
    buffer = np.empty(rows * cols, dtype=data.dtype)
    for band in range(channels):
        for row in range(0, height, rows):
            for col in range(0, width, cols):
                source = data[row : row + rows, col : col + cols, band]
                block = buffer[: source.size].reshape((1, *source.shape))
                np.copyto(block[0], source)
                dst.write(
                    block,
                    [band + 1],
                    window=rasterio.windows.Window(
                        col, row, source.shape[1], source.shape[0]
                    ),
                )


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
