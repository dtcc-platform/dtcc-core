"""A semantic entity assigned to elements of an existing geometry."""

from dataclasses import dataclass, field

import numpy as np


@dataclass
class SemanticRegion:
    """Metadata for a set of geometry elements, without copying their coordinates.

    ``indices`` selects surfaces of a MultiSurface/Solid or faces of a Mesh. ``id`` is
    optional and local to the owning geometry. Domain meaning comes from the
    semantic URI and a separately selected profile, not a Python subclass.
    ``parent`` optionally indexes another region in the same geometry's ordered
    ``regions`` list. It does not add the parent's elements to this region or
    vice versa. Parents form an acyclic forest; inverse children are derived.
    Reindex element membership when topology changes, and parent references when
    the region list is reordered. Neither relationship requires an ``id``.
    """

    def __repr__(self):
        from ...common._display import format_repr

        return format_repr(
            type(self).__name__,
            [
                ("semantic_type", self.semantic_type),
                ("id", self.id),
                ("num_elements", len(self.indices)),
                ("num_attributes", len(self.attributes)),
                ("parent", self.parent),
            ],
        )

    semantic_type: str
    indices: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int64))
    id: str | None = None
    attributes: dict = field(default_factory=dict)
    parent: int | None = None
