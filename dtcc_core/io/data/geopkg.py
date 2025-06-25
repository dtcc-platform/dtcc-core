#!/usr/bin/env python3

import asyncio
import os
import json
import aiohttp
import requests
from platformdirs import user_cache_dir
from .logging import info, warning, debug, error

CACHE_DIR = user_cache_dir(appname="dtcc-data")
os.makedirs(CACHE_DIR,exist_ok=True)
# Where we store local cache metadata
CACHE_FILE = os.path.join(CACHE_DIR,"tile_cache_superset.json")

# The FastAPI endpoint
DEFAULT_SERVER_URL = "http://127.0.0.1:8000/tiles"

try:
    import nest_asyncio
    nest_asyncio.apply()
except ImportError:
    pass

def load_cache():
    """
    Load or create an empty cache metadata from tile_cache_superset.json.
    We'll store a list of records like:
    [
      {
        "bbox": [268234.462, 6473567.915, 278234.462, 6483567.915],
        "zipfile": "tiles_268234.462_6473567.915_278234.462_6483567.915.zip"
      },
      ...
    ]
    """
    if not os.path.exists(CACHE_FILE):
        return []
    with open(CACHE_FILE, "r", encoding="utf-8") as f:
        return json.load(f)

def save_cache(cache_data):
    """Save the cache metadata to tile_cache_superset.json."""
    with open(CACHE_FILE, "w", encoding="utf-8") as f:
        json.dump(cache_data, f, indent=2)

def is_superset_bbox(bbox_sup, bbox_sub):
    """
    Return True if bbox_sub is fully contained within bbox_sup.
    Each bbox is (minx, miny, maxx, maxy).
    """
    minxS, minyS, maxxS, maxyS = bbox_sup
    minxT, minyT, maxxT, maxyT = bbox_sub
    return (
        minxS <= minxT and
        minyS <= minyT and
        maxxS >= maxxT and
        maxyS >= maxyT
    )

def find_superset_in_cache(bbox, cache_data):
    """
    Given 'bbox' and a list of records in cache_data,
    check if any record's bbox is a superset of 'bbox'.
    Return that record if found, else None.
    """
    for rec in cache_data:
        cached_bbox = tuple(rec["bbox"])
        if is_superset_bbox(cached_bbox, bbox):
            return rec
    return None

def post_gpkg_request(url, session, xmin, ymin, xmax, ymax, buffer_value=0):
    """
    Sends a POST request to the FastAPI endpoint with the given bounding box & buffer.
    Returns the parsed JSON response.
    Example: url = "http://127.0.0.1:8000/tiles"
    """
    payload = {
        "minx": xmin,
        "miny": ymin,
        "maxx": xmax,
        "maxy": ymax,
    }
    debug(f"[POST] to {url} with payload={payload}")
    resp = session.post(f'{url}/tiles', json=payload, timeout=30)
    debug(resp)
    if resp.status_code != 200:
        raise RuntimeError(
            f"Request failed with status {resp.status_code}:\n{resp.text}"
        )
    data = resp.json()
    return data

async def download_gpkg_file(session, base_url, filename, output_dir):
    """
    Download a single .laz file asynchronously with aiohttp if not already cached.
    The endpoint is assumed to be: f"{base_url}/get/lidar/{filename}"
    We'll store the downloaded file in output_dir/filename.
    """
    url = f"{base_url}/get/gpkg/{filename}"
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, filename)

    # 1) Check local cache
    if os.path.exists(out_path):
        info(f"File {filename} already in cache, skipping download.")
        return  # skip

    # 2) If not cached, download
    info(f"Downloading {filename} from {url}")
    async with session.get(url) as resp:
        if resp.status == 200:
            content = await resp.read()
            with open(out_path, "wb") as f:
                f.write(content)
            info(f"Saved {filename} to {out_path}")
        else:
            warning(f"Failed to download {filename}, status code={resp.status}")

async def download_all_gpkg_files(base_url, filenames, output_dir="downloaded_gpkg"):
    """
    Given a list of filenames, downloads them all asynchronously from
    base_url/get/lidar/<filename> using aiohttp, skipping any local cache hits.
    """
        
    async with aiohttp.ClientSession() as session:
        tasks = []
        for fname in filenames:
            tasks.append(download_gpkg_file(session, base_url, fname, output_dir))
        # Run all downloads concurrently
        await asyncio.gather(*tasks)

def run_download_files(base_url, filenames, output_dir="downloaded_gpkg"):
    """
    Entry point to run the async download with asyncio, skipping already cached files.
    """
    if not filenames:
        info("No files to download.")
        return
    debug(f"Downloading {len(filenames)} files in parallel (with cache check)...")
    asyncio.run(download_all_gpkg_files(base_url, filenames, output_dir))
    info("All downloads finished.")

def download_tiles(user_bbox, session, server_url=DEFAULT_SERVER_URL):
    """
    1) Check if any cached bounding box is a superset of 'bbox'. If so, skip request.
    2) Otherwise, POST the bounding box to 'server_url' to get a ZIP file.
    3) Save the ZIP to a local file named e.g. 'tiles_minx_miny_maxx_maxy.zip'.
    4) Add a new cache record with the bounding box and zip filename.
    5) If later bounding boxes are subsets of this one, we skip new requests.
    """
    try:
        response_data = post_gpkg_request(
            server_url,
            session,
            xmin=user_bbox[0],
            ymin=user_bbox[1],
            xmax=user_bbox[2],
            ymax=user_bbox[3],
            buffer_value=2000
        )
    except Exception as e:
        warning(f"Error occurred: {e}")
        return
    returned_tiles = response_data["tiles"]
    output_dir = os.path.join(CACHE_DIR,'downloaded-gpkg')
    # D) Download files in parallel (with local cache)
    # filenames_to_download = [tile["filename"] for tile in returned_tiles]
    run_download_files(server_url, returned_tiles, output_dir=output_dir)
    return [os.path.join(output_dir, filename) for filename in returned_tiles]


'''
def main():
    """
    Demonstration of superset-based caching for bounding boxes.
    1) We request a bigger bounding box => triggers server request (cache miss).
    2) We then request a smaller bounding box => sees superset in cache => skip.
    """
    bigger_bbox = (268000, 6473500, 278000, 6483500)  # EPSG:3006
    smaller_bbox = (269000, 6474000, 270000, 6475000) # fully inside bigger_bbox

    # First call: bigger => no superset => fetch from server
    download_tiles(bigger_bbox)

    # Second call: smaller => superset found => skip
    download_tiles(smaller_bbox)


if __name__ == "__main__":
    main()
'''
