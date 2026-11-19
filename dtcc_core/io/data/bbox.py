"""Convert bounding boxes passed to data downloads into :class:`Bounds`."""

from __future__ import annotations

import math
from collections.abc import Sequence

from ...model import Bounds


def as_bounds(bounds: Bounds | Sequence[float]) -> Bounds:
    """Return ``bounds`` as a checked :class:`Bounds`.

    Accepts a ``Bounds`` or an ordered sequence of numbers, such as a list,
    tuple or NumPy array, in the order the datasets use: ``(minx, miny, maxx,
    maxy)`` or ``(minx, miny, minz, maxx, maxy, maxz)``. Six values are not in
    the positional order of ``Bounds(xmin, ymin, xmax, ymax, zmin, zmax)``.
    A box may have zero width or height, as for a point lookup.

    Raises
    ------
    TypeError
        If ``bounds`` is not a ``Bounds`` or an ordered sequence of 4 or 6
        numbers. Unordered collections such as sets are rejected.
    ValueError
        If a value is not finite or a minimum is above its maximum, which is
        also what a six-value box in ``Bounds`` order usually produces.
    """
    if isinstance(bounds, Bounds):
        result = bounds
    else:
        values = _numbers(bounds)
        if values is None or len(values) not in (4, 6):
            raise TypeError(
                "bounds must be a Bounds or 4 values (minx, miny, maxx, maxy) or "
                "6 values (minx, miny, minz, maxx, maxy, maxz)."
            )
        if len(values) == 4:
            xmin, ymin, xmax, ymax = values
            zmin = zmax = 0.0
        else:
            xmin, ymin, zmin, xmax, ymax, zmax = values
        result = Bounds(
            xmin=xmin, ymin=ymin, zmin=zmin, xmax=xmax, ymax=ymax, zmax=zmax
        )
    _check(result)
    return result


def _numbers(bounds) -> list[float] | None:
    ordered = (
        isinstance(bounds, Sequence) and not isinstance(bounds, (str, bytes))
    ) or (hasattr(bounds, "__array__") and getattr(bounds, "ndim", 1) == 1)
    if not ordered:
        return None
    try:
        return [float(value) for value in bounds]
    except (TypeError, ValueError):
        return None


def _check(bounds: Bounds) -> None:
    values = (
        bounds.xmin,
        bounds.ymin,
        bounds.zmin,
        bounds.xmax,
        bounds.ymax,
        bounds.zmax,
    )
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"bounds {values} must be finite numbers.")
    if (
        bounds.xmin > bounds.xmax
        or bounds.ymin > bounds.ymax
        or bounds.zmin > bounds.zmax
    ):
        raise ValueError(
            f"bounds (xmin={bounds.xmin}, ymin={bounds.ymin}, zmin={bounds.zmin}, "
            f"xmax={bounds.xmax}, ymax={bounds.ymax}, zmax={bounds.zmax}) must have "
            "each minimum at or below its maximum; six values go in the order "
            "(minx, miny, minz, maxx, maxy, maxz)."
        )
