"""Output-level selection for city reconstruction."""

from ...model import GeometryType


def resolve_output_lod(lod, lod2_method):
    """Validate the requested output LOD and its roof reconstruction method."""
    if lod2_method not in ("catalogue", "ransac"):
        raise ValueError("lod2_method must be 'catalogue' or 'ransac'")
    if isinstance(lod, GeometryType):
        lod = lod.name.lower()
    if not isinstance(lod, str) or lod not in ("lod1", "lod2"):
        raise ValueError(
            "lod selects output and must be 'lod1' or 'lod2'; "
            "stored input geometry is selected automatically"
        )
    if lod != "lod2" and lod2_method != "catalogue":
        raise ValueError("lod2_method='ransac' requires lod='lod2'")
    return lod
