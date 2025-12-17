#!/usr/bin/env python3

import os
import json
import requests
import pyproj
import geopandas as gpd
from shapely.geometry import box, Polygon, LineString
from platformdirs import user_cache_dir
from .logging import info, warning, debug, error

# ------------------------------------------------------------------------
# 1) Global constants/paths
# ------------------------------------------------------------------------
BASE_CACHE_DIR = user_cache_dir(appname="dtcc-data")
os.makedirs(BASE_CACHE_DIR,exist_ok=True)
CACHE_DIR = os.path.join(BASE_CACHE_DIR,"downloaded_osm")
os.makedirs(CACHE_DIR,exist_ok=True)
# Where we store local cache metadata
CACHE_METADATA_FILE = os.path.join(BASE_CACHE_DIR,"cache_metadata.json")

OVERPASS_URL = "https://overpass-api.de/api/interpreter"

# ------------------------------------------------------------------------
# 2) Utilities for bounding boxes
# ------------------------------------------------------------------------
def is_superset_bbox(bbox_sup, bbox_sub):
    """
    Return True if bbox_sub is fully contained within bbox_sup.
    Each bbox is (xmin, ymin, xmax, ymax).
    """
    xminS, yminS, xmaxS, ymaxS = bbox_sup
    xminT, yminT, xmaxT, ymaxT = bbox_sub
    return (
        xminS <= xminT and
        yminS <= yminT and
        xmaxS >= xmaxT and
        ymaxS >= ymaxT
    )

def filter_gdf_to_bbox(gdf, bbox_3006):
    """
    Filter a GeoDataFrame (already in EPSG:3006) to the specified bounding box by intersection.
    Returns a copy of the subset.
    """
    minx, miny, maxx, maxy = bbox_3006
    bbox_poly = box(minx, miny, maxx, maxy)  # shapely
    return gdf[gdf.geometry.intersects(bbox_poly)].copy()

# ------------------------------------------------------------------------
# 3) Metadata I/O
# ------------------------------------------------------------------------
def load_cache_metadata(meta_path=CACHE_METADATA_FILE):
    """
    Load cached Overpass metadata records from disk.

    Parameters
    ----------
    meta_path : str, default CACHE_METADATA_FILE
        Path to the JSON file containing cached metadata.

    Returns
    -------
    list
        Parsed list of cache records; empty list if the file does not exist.
    """
    if not os.path.exists(meta_path):
        return []
    with open(meta_path, "r", encoding="utf-8") as f:
        return json.load(f)

def save_cache_metadata(records, meta_path=CACHE_METADATA_FILE):
    """
    Persist Overpass cache metadata to disk.

    Parameters
    ----------
    records : list
        Metadata entries to serialize.
    meta_path : str, default CACHE_METADATA_FILE
        Destination JSON file path.

    Returns
    -------
    None
        Metadata is written to ``meta_path``.
    """
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(records, f, indent=2)

def find_superset_record(bbox_3006, records):
    """
    Look for a record in 'records' whose bounding box is a superset of bbox_3006.
    Return the first match, or None if none found.
    """
    for rec in records:
        rbox = rec["bbox"]  # (xmin, ymin, xmax, ymax)
        if is_superset_bbox(rbox, bbox_3006):
            return rec
    return None

# ------------------------------------------------------------------------
# 4) Overpass logic
# ------------------------------------------------------------------------
def download_overpass_buildings(bbox_3006):
    """
    1) Transform bbox_3006 -> EPSG:4326 (lat/lon).
    2) Query Overpass for building footprints in that bounding box.
    3) Return a GeoDataFrame in EPSG:3006.
    """
    transformer = pyproj.Transformer.from_crs("EPSG:3006", "EPSG:4326", always_xy=True)
    xmin, ymin, xmax, maxy = bbox_3006
    min_lon, min_lat = transformer.transform(xmin, ymin)
    max_lon, max_lat = transformer.transform(xmax, maxy)

    query = f"""
    [out:json];
    way["building"]({min_lat},{min_lon},{max_lat},{max_lon});
    (._;>;);
    out body;
    """
    info(f"Querying Overpass for buildings in bbox={bbox_3006}")
    resp = requests.post(OVERPASS_URL, data={"data": query}, timeout=60)
    resp.raise_for_status()
    data = resp.json()

    # Parse
    nodes = {}
    for elem in data.get("elements", []):
        if elem["type"] == "node":
            nid = elem["id"]
            lat = elem["lat"]
            lon = elem["lon"]
            nodes[nid] = (lat, lon)

    footprints_ll = []
    for elem in data.get("elements", []):
        if elem["type"] == "way" and "nodes" in elem:
            refs = elem["nodes"]
            coords = []
            for r in refs:
                if r in nodes:
                    coords.append(nodes[r])  # (lat, lon)
            if len(coords) > 2:
                if coords[0] != coords[-1]:
                    coords.append(coords[0])  # close ring
                footprints_ll.append(coords)

    # Convert lat-lon -> polygons in EPSG:4326
    polygons_4326 = []
    for ring in footprints_ll:
        ring_lonlat = [(lon, lat) for (lat, lon) in ring]
        polygons_4326.append(Polygon(ring_lonlat))

    gdf_4326 = gpd.GeoDataFrame(
        {"osm_id": range(len(polygons_4326))},
        geometry=polygons_4326,
        crs="EPSG:4326"
    )
    gdf_3006 = gdf_4326.to_crs("EPSG:3006")
    return gdf_3006

def download_overpass_roads(bbox_3006):
    """
    1) Transform bbox_3006 -> EPSG:4326 (lat/lon).
    2) Query Overpass for roads (highways) in that bounding box.
    3) Return a GeoDataFrame in EPSG:3006.
    """
    transformer = pyproj.Transformer.from_crs("EPSG:3006", "EPSG:4326", always_xy=True)
    xmin, ymin, xmax, maxy = bbox_3006
    min_lon, min_lat = transformer.transform(xmin, ymin)
    max_lon, max_lat = transformer.transform(xmax, maxy)

    query = f"""
    [out:json];
    way["highway"]({min_lat},{min_lon},{max_lat},{max_lon});
    (._;>;);
    out body;
    """
    info(f"Querying Overpass for roads in bbox={bbox_3006}")
    resp = requests.post(OVERPASS_URL, data={"data": query}, timeout=60)
    resp.raise_for_status()
    data = resp.json()

    # Parse
    nodes = {}
    for elem in data.get("elements", []):
        if elem["type"] == "node":
            nid = elem["id"]
            lat = elem["lat"]
            lon = elem["lon"]
            nodes[nid] = (lat, lon)

    roads_ll = []
    for elem in data.get("elements", []):
        if elem["type"] == "way" and "nodes" in elem:
            refs = elem["nodes"]
            coords = []
            for r in refs:
                if r in nodes:
                    coords.append(nodes[r])  # (lat, lon)
            if len(coords) > 1:
                roads_ll.append(coords)

    # lat-lon -> lines in EPSG:4326
    lines_4326 = []
    for line_coords in roads_ll:
        line_lonlat = [(lon, lat) for (lat, lon) in line_coords]
        lines_4326.append(LineString(line_lonlat))

    gdf_4326 = gpd.GeoDataFrame(
        {"osm_id": range(len(lines_4326))},
        geometry=lines_4326,
        crs="EPSG:4326"
    )
    gdf_3006 = gdf_4326.to_crs("EPSG:3006")
    return gdf_3006

# ------------------------------------------------------------------------
# 5) Superset-based caching logic for Buildings
# ------------------------------------------------------------------------
def get_buildings_for_bbox(bbox_3006):
    """
    1) load metadata
    2) look for superset => filter local
    3) otherwise => Overpass => store => add to metadata
    4) return GDF in EPSG:3006
    """
    records = load_cache_metadata()
    sup_rec = find_superset_record(bbox_3006, [r for r in records if r["type"] == "buildings"])
    if sup_rec:
        debug("Found superset bounding box for buildings:", sup_rec["bbox"])
        gdf_all = gpd.read_file(sup_rec["filepath"], layer=sup_rec["layer"])
        subset_gdf = filter_gdf_to_bbox(gdf_all, bbox_3006)
        saved_filename = sup_rec['filepath']
        info(f"Subset size: {len(subset_gdf)} features for buildings in bbox={bbox_3006}")
        return subset_gdf, saved_filename
    else:
        # Overpass
        new_gdf = download_overpass_buildings(bbox_3006)
        info(f"Downloaded {len(new_gdf)} building footprints from Overpass.")
        # store
        out_filename = os.path.join(CACHE_DIR,f"buildings_{bbox_3006[0]}_{bbox_3006[1]}_{bbox_3006[2]}_{bbox_3006[3]}.gpkg")
        new_gdf.to_file(out_filename, layer="buildings", driver="GPKG")
        # update metadata
        records.append({
            "type": "buildings",
            "bbox": list(bbox_3006),
            "filepath": out_filename,
            "layer": "buildings"
        })
        save_cache_metadata(records)
        return new_gdf, out_filename

# ------------------------------------------------------------------------
# 6) Superset-based caching logic for Roads
# ------------------------------------------------------------------------
def get_roads_for_bbox(bbox_3006):
    """
    1) load metadata
    2) look for superset => filter local
    3) otherwise => Overpass => store => add to metadata
    4) return GDF in EPSG:3006
    """
    records = load_cache_metadata()
    sup_rec = find_superset_record(bbox_3006, [r for r in records if r["type"] == "roads"])
    if sup_rec:
        debug("Found superset bounding box for roads:", sup_rec["bbox"])
        gdf_all = gpd.read_file(sup_rec["filepath"], layer=sup_rec["layer"])
        saved_filename = sup_rec['filepath']
        subset_gdf = filter_gdf_to_bbox(gdf_all, bbox_3006)
        info(f"Subset size: {len(subset_gdf)} features for roads in bbox={bbox_3006}")
        return subset_gdf, saved_filename
    else:
        # Overpass
        new_gdf = download_overpass_roads(bbox_3006)
        info(f"Downloaded {len(new_gdf)} road features from Overpass.")
        # store
        out_filename = os.path.join(CACHE_DIR, f"roads_{bbox_3006[0]}_{bbox_3006[1]}_{bbox_3006[2]}_{bbox_3006[3]}.gpkg")
        new_gdf.to_file(out_filename, layer="roads", driver="GPKG")
        # update metadata
        records.append({
            "type": "roads",
            "bbox": list(bbox_3006),
            "filepath": out_filename,
            "layer": "roads"
        })
        save_cache_metadata(records)
        return new_gdf, out_filename


# ------------------------------------------------------------------------
# Example usage
# ------------------------------------------------------------------------
'''
if __name__ == "__main__":
    # bigger bounding box in EPSG:3006
    bboxA = (267000, 6519000, 270000, 6521000)
    # smaller bounding box, fully inside A
    bboxB = (268000, 6519500, 269000, 6520000)

    print("=== Buildings: BBox A ===")
    bldgA = get_buildings_for_bbox(bboxA)
    print("Buildings A size:", len(bldgA))

    print("=== Buildings: BBox B (subset) ===")
    bldgB = get_buildings_for_bbox(bboxB)
    print("Buildings B size:", len(bldgB))

    print("=== Roads: BBox A ===")
    roadsA = get_roads_for_bbox(bboxA)
    print("Roads A size:", len(roadsA))

    print("=== Roads: BBox B (subset) ===")
    roadsB = get_roads_for_bbox(bboxB)
    print("Roads B size:", len(roadsB))
'''
