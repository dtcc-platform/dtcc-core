# Copyright(C) 2023 Anders Logg
# Licensed under the MIT License

from abc import abstractmethod
from dataclasses import dataclass, field
from typing import Union

from ..model import Model
from ..values import Field
from .bounds import Bounds
from .semantic_region import SemanticRegion
from .transform import Transform


@dataclass(repr=False)
class Geometry(Model):
    """Base class for all geometry classes.

    Geometry classes represent geometric objects such as point clouds,
    surfaces, polygons, and meshes. They are used to represent the geometry of
    city objects.

    All geometries are stored in a local coordinate system, which may be
    different for each geometry. The transform attribute is used to transform
    the geometry from the local coordinate system to a global coordinate system.

    A geometry may have a list of fields which are values (scalars or vectors)
    defined on the entities of the geometry (vertices, edges, faces, etc.).

    Attributes
    ----------
    bounds : Bounds
        Bounding box of the geometry in the local coordinate system.
    transform : Transform
        Affine transform to a global coordinate system.
    fields: list[Field]
    """

    _bounds: Bounds = field(default_factory=Bounds)
    transform: Transform = field(default_factory=Transform)
    fields: list[Field] = field(default_factory=list)
    regions: list[SemanticRegion] = field(default_factory=list, kw_only=True)

    def _summary_items(self):
        items = [("num_fields", len(self.fields))]
        if self.regions:
            items.append(("num_regions", len(self.regions)))
        return items

    def _info_sections(self):
        from dataclasses import fields

        import numpy as np

        from .._display import field_section

        sections = super()._info_sections()
        sections[0][2].extend(
            [
                ("Bounds (local)", self.bounds.bndstr),
                ("CRS", self.transform.srs or "Not specified"),
                (
                    "Transform",
                    "Identity"
                    if np.array_equal(self.transform.affine, np.eye(4))
                    else str(self.transform.affine),
                ),
            ]
        )
        arrays = [
            (f.name, value.shape, value.dtype)
            for f in fields(self)
            if not f.name.startswith("_")
            and isinstance(value := getattr(self, f.name), np.ndarray)
        ]
        if arrays:
            sections.append(("Arrays", ("Name", "Shape", "Type"), arrays))
        if self.fields:
            sections.append(field_section(self.fields))
        if self.regions:
            sections.append(
                (
                    "Semantic regions",
                    ("Type", "ID", "Elements", "Parent"),
                    [
                        (r.semantic_type, r.id, len(r.indices), r.parent)
                        for r in self.regions
                    ],
                )
            )
        return sections

    @abstractmethod
    def calculate_bounds(self):
        """
        Compute and cache the bounding box for the geometry.

        Returns
        -------
        Bounds
            Calculated bounds for the geometry.
        """
        pass

    @property
    def bounds(self) -> Bounds:
        """
        Bounding box of the geometry in local coordinates.

        Public array edits are not observed. Call ``calculate_bounds()`` to
        refresh derived geometry bounds explicitly; transforms are not applied.
        Grid domain bounds are intrinsic state and should not be recalculated
        as though they were a derived coordinate cache.

        Returns
        -------
        Bounds
            Cached bounds; computed on demand when missing.
        """
        if self._bounds is None or self._bounds.area == 0:
            self.calculate_bounds()
        return self._bounds

    @bounds.setter
    def bounds(self, bounds: Bounds):
        """
        Set the cached bounds for the geometry.

        Parameters
        ----------
        bounds : Bounds
            Bounding box to assign.
        """
        self._bounds = bounds

    def add_field(self, field: Field):
        """Add a field to the geometry.

        Parameters
        ----------
        field : Field
            The field to add to the geometry.
        """
        self.fields.append(field)

    def regions_of(self, semantic_type: str) -> list[SemanticRegion]:
        """Return regions with this exact semantic URI, retaining native arrays."""
        return [
            region for region in self.regions if region.semantic_type == semantic_type
        ]
