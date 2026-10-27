# Copyright(C) 2024 Dag Wästberg
# Licensed under the MIT License

from dataclasses import dataclass

from .object import Object


@dataclass(repr=False)
class Terrain(Object):
    """Represents a terrain object in a city."""
