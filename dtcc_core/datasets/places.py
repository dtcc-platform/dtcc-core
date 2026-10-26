from pathlib import Path
import tempfile
from typing import List, Literal, Optional

import dtcc_core
from dtcc_core.model import PlaceCollection
from pydantic import Field

from ..common import warning

from .dataset import DatasetBaseArgs, DatasetDescriptor
from .providers import provider_entry


Category = Literal[
    "accommodation", "craft", "culture", "education", "finance", "food_and_drink",
    "health", "leisure", "office", "other", "public_service", "shopping", "tourism",
]


class PlacesArgs(DatasetBaseArgs):
    """Arguments for the OpenStreetMap places dataset."""

    categories: Optional[List[Category]] = Field(
        None, description="Keep only these place categories; all categories by default."
    )
    roads: bool = Field(
        True,
        description=(
            "Also download OpenStreetMap roads for the same bounds and keep them "
            "as context_roads, drawn underneath the places by plot(). Ignored for "
            "geojson/gpkg output."
        ),
    )
    format: Optional[Literal["geojson", "gpkg"]] = Field(
        None, description="Output format (geojson/gpkg for point vector bytes)"
    )


class PlacesDataset(DatasetDescriptor):
    """Places dataset, registered as ``places``.

    Named places of interest from OpenStreetMap, such as shops, restaurants,
    schools and clinics, for the requested bounds, returned as a
    PlaceCollection in EPSG:3006.

    Call it with keyword arguments defined by ``PlacesArgs``.
    """

    name = "places"
    title = "Places of Interest"
    description = (
        "Named places of interest from OpenStreetMap, such as shops, restaurants, "
        "schools and clinics, for the requested bounds, returned as a "
        "PlaceCollection in EPSG:3006."
    )
    ArgsModel = PlacesArgs
    data_category = "raw"
    result_kind = "place_collection"
    python_return_type = "dtcc_core.model.PlaceCollection"
    default_crs = "EPSG:3006"
    provider = [
        provider_entry("openstreetmap", role="source_provider"),
        provider_entry("dtcc-platform", role="processor"),
    ]
    source = [
        {
            "name": "OpenStreetMap via public Overpass API endpoints",
            "role": "source_provider",
            "service": "overpass-api.de and mirrors",
            "source_terms_status": "requires_review",
        }
    ]
    license = (
        "OpenStreetMap data, © OpenStreetMap contributors, Open Database License "
        "(ODbL) 1.0. Redistributed derived databases must be offered under ODbL "
        "with attribution. Requires review before redistribution."
    )
    collection_period = "Current OpenStreetMap data at query time"
    data_types = ["vector", "points", "places"]
    geographic_coverage = "Global, constrained by requested bounds"
    update_frequency = "Continuous (live OpenStreetMap edits)"
    processing_steps = [
        "Convert the bounds to WGS84 using all four corners",
        "Query Overpass for named features with amenity, shop, tourism, leisure, healthcare, office or craft tags",
        "Drop street furniture such as benches, post boxes and parking spaces",
        "Place ways and relations at their centre and keep features inside the bounds",
        "Assign a coarse category from the OSM key and value",
        "Download OpenStreetMap roads for the same bounds as plot context, unless roads=False",
        "Return a PlaceCollection or serialize GeoJSON/GeoPackage bytes",
    ]
    presentation_headline = "Places of Interest"
    presentation_summary = (
        "Shops, restaurants, schools, clinics and other named places from "
        "OpenStreetMap, grouped into simple categories."
    )
    presentation_narrative = [
        {
            "heading": "What you are seeing",
            "body": (
                "Each point is a named place mapped in OpenStreetMap, with its "
                "category, original OSM tag and, where mapped, opening hours, "
                "cuisine, website and address."
            ),
        },
        {
            "heading": "How to interpret it",
            "body": (
                "Coverage depends on local mappers: city centres are usually well "
                "mapped, outskirts less so. Places without a name are left out."
            ),
        },
        {
            "heading": "Limitations",
            "body": (
                "Closed businesses may remain until someone updates the map, and "
                "categories are a coarse grouping of OSM tags."
            ),
        },
    ]
    key_points = [
        "Named shops, restaurants, schools, clinics, offices and more",
        "Coarse categories plus the original OSM key and value",
        "Opening hours, cuisine, website and address where mapped",
        "ODbL licence; attribution required",
    ]
    presentation_legend = {
        "title": "Place attributes",
        "entries": [
            {"label": "category", "meaning": "Coarse place category, such as food_and_drink or shopping"},
            {"label": "osm_key / osm_value", "meaning": "Original OSM tag, such as amenity=restaurant"},
            {"label": "osm_id / osm_type", "meaning": "OpenStreetMap element reference"},
        ],
    }
    view_hints = {
        "preferred_geometry": "points",
        "default_crs": "EPSG:3006",
        "default_color_attribute": "category",
    }
    presentation_warnings = [
        "OpenStreetMap attribution and ODbL terms must accompany redistributed data.",
        "Public Overpass endpoints can be slow or busy; results are cached per bounds.",
    ]
    presentation_limitations = [
        "Completeness varies with local mapping activity.",
        "Unnamed places and street furniture are excluded.",
    ]

    def build(self, args: PlacesArgs):
        bounds = self.parse_bounds(args.bounds)
        if args.format in ("geojson", "gpkg"):
            gdf = dtcc_core.io.data.osm_places.download_places_geodataframe(
                bounds=bounds, categories=args.categories
            )
            return self._export_gdf_to_bytes(gdf, args.format)

        places: PlaceCollection = dtcc_core.io.data.download_places(
            bounds=bounds, categories=args.categories
        )
        if args.roads:
            try:
                places.context_roads = dtcc_core.io.data.download_roadnetwork(
                    bounds=bounds, provider="OSM", epsg="3006"
                )
            except Exception as e:
                warning(f"Could not download context roads for places; plotting without them: {e}")
        return places

    @staticmethod
    def _export_gdf_to_bytes(gdf, format: str) -> bytes:
        with tempfile.TemporaryDirectory() as tmpdir:
            if format == "geojson":
                path = Path(tmpdir) / "places.geojson"
                gdf.to_file(path, driver="GeoJSON")
            else:
                path = Path(tmpdir) / "places.gpkg"
                gdf.to_file(path, layer="places", driver="GPKG")
            return path.read_bytes()


__all__ = ["PlacesArgs", "PlacesDataset"]
