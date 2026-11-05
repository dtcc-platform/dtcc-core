"""Convert bounding boxes passed to data downloads into :class:`Bounds`."""

from __future__ import annotations

from collections.abc import Iterable, Mapping

from ...model import Bounds


def as_bounds(bounds: Bounds | Iterable[float]) -> Bounds:
    """Return ``bounds`` as a :class:`Bounds`.

    Accepts a ``Bounds`` or any sequence of numbers, such as a list, tuple or
    NumPy array, in the order the datasets use: ``(minx, miny, maxx, maxy)``
    or ``(minx, miny, minz, maxx, maxy, maxz)``. Six values are not in the
    positional order of ``Bounds(xmin, ymin, xmax, ymax, zmin, zmax)``.

    Raises
    ------
    TypeError
        If ``bounds`` is not a ``Bounds`` or a sequence of 4 or 6 numbers.
    ValueError
        If a minimum is not below its maximum, which is also what a six-value
        box in ``Bounds`` order usually produces.
    """
    if isinstance(bounds, Bounds):
        return bounds
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
    if not (xmin < xmax and ymin < ymax and zmin <= zmax):
        raise ValueError(
            f"bounds {tuple(values)} must have each minimum below its maximum, "
            "in the order (minx, miny[, minz], maxx, maxy[, maxz])."
        )
    return Bounds(xmin=xmin, ymin=ymin, zmin=zmin, xmax=xmax, ymax=ymax, zmax=zmax)


def _numbers(bounds) -> list[float] | None:
    if isinstance(bounds, (str, bytes, Mapping)) or not isinstance(bounds, Iterable):
        return None
    try:
        return [float(value) for value in bounds]
    except (TypeError, ValueError):
        return None
