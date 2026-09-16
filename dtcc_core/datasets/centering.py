"""Centring exported geometry on the world origin.

National reference systems put geometry a long way from the origin: a Swedish
northing in SWEREF99 TM is about 6.4 million metres. Binary STL stores 32-bit
floats, which resolve only about half a metre at that magnitude, so a building
is already distorted inside the file. CAD tools such as Rhino, AutoCAD, 3ds Max
and Blender also lose precision the further geometry sits from the origin.

Centring uses the centre of the *requested bounds* rather than the extent of
the result. The requested bounds are known to the caller and recorded in the
dataset manifest, so the translation can be undone without extra metadata, and
datasets requested for the same bounds stay aligned with each other.
"""

from __future__ import annotations

from typing import Sequence

from ..common import warning

CENTER_ON_ORIGIN_DESCRIPTION = (
    "Translate exported geometry so the centre of the requested bounds sits at "
    "x=0, y=0. Heights are left alone. Keeps coordinates small for CAD tools and "
    "for the 32-bit floats used by STL. Ignored for raster output."
)

CENTER_ON_ORIGIN_PROCESSING_STEP = (
    "For center_on_origin=True, translate exported geometry so the centre of the "
    "requested bounds is at x=0, y=0, leaving heights unchanged"
)

# Beyond this magnitude a 32-bit float can no longer resolve a centimetre.
LARGE_COORDINATE_LIMIT = 100_000.0


def bounds_center_offset(bounds: Sequence[float]) -> tuple[float, float, float]:
    """Return the offset that moves the centre of ``bounds`` to the origin."""
    values = [float(value) for value in bounds]
    if len(values) == 4:
        xmin, ymin, xmax, ymax = values
    elif len(values) == 6:
        xmin, ymin, _zmin, xmax, ymax, _zmax = values
    else:
        raise ValueError("Bounds must be 4 or 6 floats")
    return (-(xmin + xmax) / 2.0, -(ymin + ymax) / 2.0, 0.0)


def center_geometry_on_origin(geometry, bounds: Sequence[float]):
    """Translate ``geometry`` in place onto the centre of ``bounds``."""
    return geometry.offset(bounds_center_offset(bounds))


def center_result_if_requested(result, args):
    """Centre a mesh result when the caller asked for it.

    Results without vertices to move, such as a raster, are returned as they
    are: centring only means something for geometry.
    """
    if not getattr(args, "center_on_origin", False):
        return result
    if not (hasattr(result, "offset") and hasattr(result, "vertices")):
        return result
    return center_geometry_on_origin(result, args.bounds)


def warn_if_far_from_origin(result, format: str | None) -> None:
    """Warn when STL output carries coordinates too large for 32-bit floats."""
    if format != "stl":
        return
    try:
        bounds = result.bounds
        largest = max(
            abs(bounds.xmin), abs(bounds.xmax), abs(bounds.ymin), abs(bounds.ymax)
        )
    except (AttributeError, TypeError, ValueError):
        return
    if largest > LARGE_COORDINATE_LIMIT:
        warning(
            f"Exporting STL with coordinates up to {largest:.0f} m from the origin. "
            "STL stores 32-bit floats, so vertices may shift by decimetres and CAD "
            "tools lose precision. Pass center_on_origin=True to avoid this."
        )
