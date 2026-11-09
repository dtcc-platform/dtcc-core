import warnings
from functools import cache, partial
from pathlib import Path

import numpy as np
import rasterio
import rasterio.merge
from PIL import Image

from ..model import Raster
from . import generic
from .logging import error, warning
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

    if isinstance(path, (str, Path)) and Path(path).suffix.lower() == ".dtcc":
        return load_model(path, expected_type=Raster, validate_schema=validate_schema)
    if validate_schema is not True:
        raise ValueError("validate_schema applies to .dtcc input")
    return generic.load(path, "raster", Raster, _load_formats, delimiter=delimiter)


def _save_json_raster(raster, path):
    path.write_text(raster.to_json())


# Largest array copied at once when writing (height, width, channels) data.
WRITE_STRIP_BYTES = 8 * 1024**2


def _save_geotif(raster, path, alpha=False):
    """Write a GeoTIFF with the raster's CRS and nodata, if set.

    (height, width, channels) data is written in windows of at most
    WRITE_STRIP_BYTES, all bands of a window at once, so no band-first copy of
    the whole array is made and each compressed strip is written once.
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
    pixel_bytes = channels * data.dtype.itemsize
    cols = max(1, min(width, WRITE_STRIP_BYTES // pixel_bytes))
    rows = max(1, min(height, WRITE_STRIP_BYTES // (cols * pixel_bytes)))
    # One buffer for every window; each window uses a contiguous prefix of it.
    # Rasterio writes a C-contiguous (bands, rows, cols) array with a list of
    # band indexes without copying it. Writing all bands of a window together
    # keeps GDAL from rewriting a pixel-interleaved strip once per band when its
    # block cache cannot hold the image.
    buffer = np.empty(channels * rows * cols, dtype=data.dtype)
    indexes = list(range(1, channels + 1))
    for row in range(0, height, rows):
        for col in range(0, width, cols):
            source = data[row : row + rows, col : col + cols]
            block = buffer[: source.size].reshape((channels, *source.shape[:2]))
            np.copyto(block, np.moveaxis(source, -1, 0))
            dst.write(
                block,
                indexes,
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


def valid_pixels(data, nodata):
    """Mask of (height, width) pixels that hold data.

    A pixel is invalid when it is NaN or equals a finite ``nodata``; a pixel
    of a multi-band array is invalid when every band is. Integer data is
    compared exactly, so only an integer-valued ``nodata`` within the dtype's
    range matches any pixel.
    """
    invalid = np.isnan(data) if data.dtype.kind == "f" else np.zeros(data.shape, bool)
    if data.dtype.kind in "iu":
        if np.isfinite(nodata) and float(nodata).is_integer():
            value = int(nodata)
            info = np.iinfo(data.dtype)
            if info.min <= value <= info.max:
                invalid |= data == value
    elif not np.isnan(nodata):
        invalid |= data == nodata
    return ~invalid.all(axis=-1) if data.ndim == 3 else ~invalid


def _placement(crs, transform, width, height, nodata):
    from rasterio.warp import transform_bounds

    corners = [
        transform * c for c in ((0, 0), (width, 0), (0, height), (width, height))
    ]
    xs, ys = zip(*corners, strict=True)
    bounds = [min(xs), min(ys), max(xs), max(ys)]
    return {
        "crs": crs,
        "transform": [float(value) for value in list(transform)[:6]],
        "width": int(width),
        "height": int(height),
        "bounds": [float(value) for value in bounds],
        "bounds_wgs84": list(transform_bounds(crs, "EPSG:4326", *bounds)),
        "nodata": "NaN" if np.isnan(nodata) else float(nodata),
    }


@cache
def _warn_approximate_placement():
    """Warn, once per process, that PNG placement is approximate."""
    warnings.warn(
        f"rasterio {rasterio.__version__} cannot place PNG pixels exactly; "
        "placement is approximate to within 1/8 source pixel. Python 3.12+ with "
        "rasterio 1.5+ places exactly.",
        UserWarning,
        stacklevel=3,
    )


def write_placed_png(raster, path, *, band=None, value_range=None, colormap=None):
    """Write a north-up EPSG:3857 RGBA PNG of ``raster`` and return its placement.

    The grid comes from ``rasterio.warp.calculate_default_transform``. With
    ``band=None`` the uint8 RGB(A) bands are written as they are; otherwise
    that band is mapped through the named matplotlib ``colormap`` over
    ``value_range``, with values at or below the low end in the low colour.
    Integer data is resampled nearest and floating data bilinear; the
    validity mask is resampled nearest and becomes alpha (existing alpha is
    kept), and invalid or uncovered pixels are fully transparent black.
    """
    from inspect import signature

    import matplotlib
    from rasterio.control import GroundControlPoint
    from rasterio.warp import Resampling, calculate_default_transform, reproject

    # Exact transformation: the approximate default moves pixel picks. Older
    # rasterio lacks the parameter and would pass it on to GDAL unchecked.
    exact = {"tolerance": 0} if "tolerance" in signature(reproject).parameters else {}
    if not exact:
        _warn_approximate_placement()
    width, height = raster.width, raster.height
    corners = ((0, 0), (width, 0), (0, height), (width, height))
    gcps = [
        GroundControlPoint(row, col, *(raster.georef * (col, row)), z=0.0)
        for col, row in corners
    ]
    transform, dst_width, dst_height = calculate_default_transform(
        raster.crs, "EPSG:3857", width, height, gcps=gcps
    )

    def warp(source, resampling, nodata=None):
        fill = np.nan if source.dtype.kind == "f" else 0
        destination = np.full((dst_height, dst_width), fill, dtype=source.dtype)
        reproject(
            source,
            destination,
            src_transform=raster.georef,
            src_crs=raster.crs,
            src_nodata=nodata,
            dst_transform=transform,
            dst_crs="EPSG:3857",
            dst_nodata=nodata,
            resampling=resampling,
            **exact,
        )
        return destination

    rgb = band is None
    data = raster.data if rgb or raster.data.ndim == 2 else raster.data[..., band]
    valid = valid_pixels(data, raster.nodata)
    mask = warp(valid.astype(np.uint8), Resampling.nearest) == 1
    if rgb:
        bands = [warp(data[..., i], Resampling.nearest) for i in range(data.shape[2])]
        if len(bands) == 3:
            bands.append(np.full_like(bands[0], 255))
        rgba = np.stack(bands, axis=-1)
    else:
        source = np.where(valid, data, np.nan).astype(np.float64)
        resampling = (
            Resampling.nearest if data.dtype.kind in "iu" else Resampling.bilinear
        )
        values = warp(source, resampling, nodata=np.nan)
        low, high = value_range
        if high > low:
            values = np.clip((values - low) / (high - low), 0.0, 1.0)
        else:
            values = np.zeros_like(values)
        rgba = matplotlib.colormaps[colormap](values, bytes=True)
    rgba[~mask] = 0
    Image.fromarray(rgba).save(path, format="PNG")
    return _placement("EPSG:3857", transform, dst_width, dst_height, raster.nodata)


def write_georeferenced_tif(raster, path, *, unit=None, description=None):
    """Write ``raster`` on its own grid as a GeoTIFF and return its placement.

    The file carries the CRS, geotransform, every band and the nodata value
    (NaN only for floating data); ``unit`` and ``description`` are set on every
    band when given.
    """
    nodata = raster.nodata
    profile = dict(
        driver="GTiff",
        height=raster.height,
        width=raster.width,
        count=raster.channels,
        dtype=raster.data.dtype,
        crs=raster.crs,
        transform=raster.georef,
        compress="DEFLATE",
    )
    if not np.isnan(nodata) or raster.data.dtype.kind == "f":
        profile["nodata"] = nodata
    if raster.channels > 1:
        profile["alpha"] = "UNSPECIFIED"
    with rasterio.open(path, "w", **profile) as dst:
        if raster.channels == 1:
            dst.write(raster.data.reshape(raster.height, raster.width), 1)
        else:
            _write_bands_in_windows(dst, raster.data)
        for index in range(1, raster.channels + 1):
            if unit is not None:
                dst.set_band_unit(index, unit)
            if description is not None:
                dst.set_band_description(index, description)
    return _placement(raster.crs, raster.georef, raster.width, raster.height, nodata)


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
