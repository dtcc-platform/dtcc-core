"""Checks for attached data that an export format has no mapping for."""


def reject_time_series(geometries, target):
    """Raise ``NotImplementedError`` if any geometry carries time series.

    Exporters without a time series mapping call this before writing, so
    attached time series are refused instead of being left out of the output.
    """
    for geometry in geometries:
        if getattr(geometry, "time_series", None):
            raise NotImplementedError(
                f"Time series require an explicit {target} mapping; remove them "
                "from the geometry or save the model as .dtcc."
            )
