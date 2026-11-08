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
    None when unknown. ``archive`` is the local archive that holds ``path`` as
    a member, or None. A record holds no pixels and never opens its file.
    A tile restored from canonical exchange has a package-relative reference
    as ``path`` and no package binding, and refuses to load.
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
    archive: Path | None = None
    # Runtime state, never serialized: path is a package reference with no binding.
    _unbound: bool = field(default=False, repr=False, compare=False)

    def load(self, *, max_memory_bytes: int = 2 * 1024**3):
        """Load every band of the local file into an (H, W, C) Raster.

        Reads only the local file, within an estimated working-memory budget.
        A tile restored without a package binding raises ValueError.
        """
        if self._unbound:
            raise ValueError(
                f"Cannot load raster tile {self.id!r}: its reference "
                f"{self.path.as_posix()!r} has no package binding. Load the tile "
                "from the dataset package that contains it."
            )
        from ...io.raster_tiles import load_tile

        return load_tile(self, max_memory_bytes=max_memory_bytes)


@dataclass(repr=False)
class RasterTileCollection(Model):
    """Sequence of file-backed raster tiles.

    Canonical exchange carries the tile metadata and package-relative
    references, never local cache paths. The tiles refer to files in a local
    cache, so the collection cannot be exported or published.
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

    def publish(self, *args, **kwargs):
        """Refuse publishing, before any upload configuration is read."""
        raise NotImplementedError(
            "Publishing a RasterTileCollection is not supported: its tiles are "
            "local cache files."
        )
