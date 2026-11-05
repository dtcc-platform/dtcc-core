"""Convert bounding boxes passed to data downloads into :class:`Bounds`."""

from __future__ import annotations

from collections.abc import Sequence

from ...model import Bounds


def as_bounds(bounds: Bounds | Sequence[float]) -> Bounds:
    """Return ``bounds`` as a :class:`Bounds`.

    Accepts a ``Bounds``, ``[minx, miny, maxx, maxy]``, or
    ``[minx, miny, minz, maxx, maxy, maxz]``, the two orders the datasets use.
    """
    if isinstance(bounds, Bounds):
        return bounds
    if isinstance(bounds, Sequence) and not isinstance(bounds, str):
        values = [float(value) for value in bounds]
        if len(values) == 4:
            xmin, ymin, xmax, ymax = values
            return Bounds(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax)
        if len(values) == 6:
            xmin, ymin, zmin, xmax, ymax, zmax = values
            return Bounds(
                xmin=xmin, ymin=ymin, zmin=zmin, xmax=xmax, ymax=ymax, zmax=zmax
            )
    raise TypeError(
        "bounds must be a Bounds or a sequence of 4 values (minx, miny, maxx, maxy) "
        "or 6 values (minx, miny, minz, maxx, maxy, maxz)."
    )
