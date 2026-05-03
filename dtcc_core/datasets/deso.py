from pathlib import Path
import tempfile
from typing import Literal, Optional

import dtcc_core
from dtcc_core.model import DeSO
from pydantic import Field

from .dataset import DatasetBaseArgs, DatasetDescriptor


class DeSOArgs(DatasetBaseArgs):
    source: Literal["SCB"] = Field("SCB", description="Data source")
    year: Literal[2018, 2025] = Field(2025, description="DeSO geometry vintage")
    format: Optional[Literal["pb", "geojson", "gpkg"]] = Field(
        None,
        description="Output format (pb for protobuf bytes, geojson/gpkg for vector bytes)",
    )


class DeSODataset(DatasetDescriptor):
    name = "deso"
    description = "Swedish DeSO demographic statistical areas from SCB."
    ArgsModel = DeSOArgs
    data_category = "raw"
    result_kind = "administrative_areas"
    python_return_type = "dtcc_core.model.DeSO"

    def build(self, args: DeSOArgs):
        bounds = self.parse_bounds(args.bounds)
        if args.format in ("geojson", "gpkg"):
            gdf = dtcc_core.io.data.deso.download_deso_geodataframe(
                bounds=bounds,
                year=args.year,
                source=args.source,
            )
            return self._export_gdf_to_bytes(gdf, args.format)

        deso: DeSO = dtcc_core.io.data.download_deso(
            bounds=bounds,
            year=args.year,
            source=args.source,
        )
        if args.format == "pb":
            return deso.to_proto().SerializeToString()
        return deso

    @staticmethod
    def _export_gdf_to_bytes(gdf, format: str) -> bytes:
        with tempfile.TemporaryDirectory(delete=True) as tmpdir:
            if format == "geojson":
                path = Path(tmpdir) / "deso.geojson"
                gdf.to_file(path, driver="GeoJSON")
            elif format == "gpkg":
                path = Path(tmpdir) / "deso.gpkg"
                gdf.to_file(path, layer="deso", driver="GPKG")
            else:
                raise ValueError(f"Unsupported DeSO export format: {format}")
            return path.read_bytes()


__all__ = ["DeSOArgs", "DeSODataset"]
