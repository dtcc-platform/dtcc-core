from enum import Enum, auto


class RoofType(Enum):
    FLAT = auto()
    GABLED = auto()
    HIPPED = auto()
    UNKNOWN = auto()


class SurfaceSemantic(Enum):
    GROUND = 0
    WALL = 1
    ROOF = 2
