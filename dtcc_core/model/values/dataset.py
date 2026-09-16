"""Transitional Dataset v2 value containers.

These generic containers are fallback model objects for unresolved dataset
return shapes. They are not the target public Dataset v2 return type when a
domain-specific model can be introduced.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from typing import Any

from ..model import Model


@dataclass(repr=False)
class DatasetCollection(Model):
    """Dataset-aware sequence container for transitional dataset returns."""

    items: list[Any] = field(default_factory=list)

    def _info_sections(self):
        from .._display import sample_section

        sections = super()._info_sections()
        if self.items:
            sections.append(sample_section("Items", self.items))
        return sections

    def _summary_items(self):
        return [("num_items", len(self.items))]

    def __len__(self) -> int:
        return len(self.items)

    def __iter__(self) -> Iterator[Any]:
        return iter(self.items)

    def __getitem__(self, index):
        return self.items[index]

    def to_list(self) -> list[Any]:
        """Return the contained items as a plain list."""
        return list(self.items)


@dataclass(repr=False)
class DatasetValue(Model):
    """Dataset-aware JSON-like value container for transitional fallback returns."""

    value: Any = None

    def _info_sections(self):
        from ...common._display import value_text

        sections = super()._info_sections()
        if isinstance(self.value, Mapping):
            sections.append(("Values", ("Key", "Value"),
                             [(key, value_text(value)) for key, value in self.value.items()]))
        else:
            sections[0][2].append(("Value", value_text(self.value)))
        return sections

    def _summary_items(self):
        return [("value_type", type(self.value).__name__)]

    def to_python(self) -> Any:
        """Return the contained Python value."""
        return self.value

    def __getitem__(self, key):
        return self.value[key]

    def __contains__(self, key) -> bool:
        return key in self.value

    def __len__(self) -> int:
        return len(self.value)

    def __iter__(self):
        return iter(self.value)

    def __eq__(self, other) -> bool:
        if isinstance(other, DatasetValue):
            return self.value == other.value
        return self.value == other

    def get(self, key, default=None):
        """Return ``mapping.get(key, default)`` for mapping values."""
        return self._mapping().get(key, default)

    def keys(self):
        """Return mapping keys for mapping values."""
        return self._mapping().keys()

    def items(self):
        """Return mapping items for mapping values."""
        return self._mapping().items()

    def values(self):
        """Return mapping values for mapping values."""
        return self._mapping().values()

    def _mapping(self) -> Mapping:
        if not isinstance(self.value, Mapping):
            raise TypeError(
                f"DatasetValue contains {type(self.value).__name__}, not a mapping."
            )
        return self.value
