"""File-backed raster tiles: local original files and their metadata."""

from collections.abc import Iterator
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from ..geometry.bounds import Bounds
from ..model import Model


@dataclass(frozen=True)
class RasterTile:
    """One original raster file in a local cache, with its source metadata.

    ``extent`` is the grid cell the file covers, not its imaged area.
    ``resolution`` (metres per pixel) and ``size_bytes`` (download size) are
    None when unknown. A record holds no pixels and never opens its file.
    """

    path: Path
    id: str
    collection: str
    datetime: datetime
    extent: Bounds
    crs: str
    spektraltyp: str
    resolution: float | None
    size_bytes: int | None


@dataclass(repr=False)
class RasterTileCollection(Model):
    """Sequence of file-backed raster tiles.

    The tiles refer to files in a local cache, so the collection cannot be
    exported, published or serialized.
    """

    tiles: list[RasterTile] = field(default_factory=list)

    def _summary_items(self):
        return [("num_tiles", len(self.tiles))]

    def __len__(self) -> int:
        return len(self.tiles)

    def __iter__(self) -> Iterator[RasterTile]:
        return iter(self.tiles)

    def __getitem__(self, index):
        return self.tiles[index]

    def export(self, *args, **kwargs):
        """Refuse export: the tiles are local cache files, not portable artifacts."""
        raise NotImplementedError(
            "Exporting a RasterTileCollection is not supported: its tiles are "
            "local cache files. Use each tile's path directly."
        )
