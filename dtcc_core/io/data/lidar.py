# Copyright (C) 2024 Vasilis Naserentin
# Licensed under the MIT License

import os
import time
import requests
import folium
import pyproj
import asyncio
import aiohttp
from platformdirs import user_cache_dir
from .logging import info, debug, warning, error
from dtcc_core.common.progress import report_progress

try:
    import nest_asyncio
    nest_asyncio.apply()
except ImportError:
    pass
# -----------------------------------------------------------------------
# Reusable Helper Functions
# -----------------------------------------------------------------------

_REQUEST_TIMEOUT_SECONDS = 30
_REQUEST_MAX_ATTEMPTS = 4
_REQUEST_RETRY_BACKOFF_SECONDS = 2.0


def _env_int(name, default, minimum=None):
    value = os.environ.get(name)
    if value is None:
        return default

    try:
        parsed = int(value)
    except ValueError:
        warning(f"Ignoring invalid {name}={value!r}; expected an integer.")
        return default

    if minimum is not None and parsed < minimum:
        warning(
            f"Ignoring invalid {name}={value!r}; expected an integer >= {minimum}."
        )
        return default

    return parsed


_DOWNLOAD_TOTAL_TIMEOUT_SECONDS = _env_int("DTCC_LIDAR_DOWNLOAD_TOTAL_TIMEOUT", 120, minimum=1)
_DOWNLOAD_CONNECT_TIMEOUT_SECONDS = _env_int("DTCC_LIDAR_DOWNLOAD_CONNECT_TIMEOUT", 30, minimum=1)
_DOWNLOAD_SOCK_READ_TIMEOUT_SECONDS = _env_int("DTCC_LIDAR_DOWNLOAD_SOCK_READ_TIMEOUT", 120, minimum=1)
_DOWNLOAD_MAX_ATTEMPTS = _env_int("DTCC_LIDAR_DOWNLOAD_MAX_ATTEMPTS", 4, minimum=1)
_DOWNLOAD_RETRY_BACKOFF_SECONDS = _env_int("DTCC_LIDAR_DOWNLOAD_RETRY_BACKOFF", 5, minimum=0)
_DOWNLOAD_RATE_LIMIT_BACKOFF_SECONDS = _env_int(
    "DTCC_LIDAR_DOWNLOAD_RATE_LIMIT_BACKOFF", 30, minimum=0
)
_DOWNLOAD_MAX_CONCURRENCY = _env_int("DTCC_LIDAR_DOWNLOAD_MAX_CONCURRENCY", 3, minimum=1)


class LidarDownloadError(RuntimeError):
    def __init__(self, filename, status, retry_after=None):
        self.filename = filename
        self.status = status
        self.retry_after = retry_after
        super().__init__(f"Download failed for {filename} with status {status}")


def post_lidar_request(url, session, xmin, ymin, xmax, ymax, buffer_value=0):
    """
    Sends a POST request to the FastAPI endpoint with the given bounding box & buffer.
    Returns the parsed JSON response.
    Example: url = "http://127.0.0.1:8000/get_lidar"
    """
    payload = {
        "xmin": round(xmin),
        "ymin": round(ymin),
        "xmax": round(xmax),
        "ymax": round(ymax),
        "buffer": buffer_value
    }
    debug(f"[POST] to {url} with payload={payload}")
    last_error = None
    for attempt in range(1, _REQUEST_MAX_ATTEMPTS + 1):
        try:
            resp = session.post(
                url,
                json=payload,
                timeout=_REQUEST_TIMEOUT_SECONDS,
            )
            debug(resp)
            if resp.status_code != 200:
                raise RuntimeError(
                    f"Request failed with status {resp.status_code}:\n{resp.text}"
                )
            return resp.json()
        except (requests.RequestException, RuntimeError) as exc:
            last_error = exc
            if attempt >= _REQUEST_MAX_ATTEMPTS:
                raise
            backoff = _REQUEST_RETRY_BACKOFF_SECONDS * attempt
            warning(
                "Lidar tile lookup attempt %d/%d failed: %s. Retrying in %.1fs.",
                attempt,
                _REQUEST_MAX_ATTEMPTS,
                exc,
                backoff,
            )
            time.sleep(backoff)

    raise RuntimeError(f"Lidar tile lookup failed: {last_error}")


def plot_bboxes_folium(user_bbox, tiles, out_html="client_map.html", crs_from="EPSG:3006"):
    """
    Plots:
      - The user_bbox (in green) in one layer
      - Each tile bounding box from 'tiles' (in blue) in another layer
    on a Folium map, with a LayerControl so you can toggle them on/off.

    'user_bbox' is a tuple (xmin, ymin, xmax, ymax) in EPSG:3006.
    'tiles' is a list of dicts, each like:
        { "filename": "...", "xmin": ..., "ymin": ..., "xmax": ..., "ymax": ... }
    """
    xmin_u, ymin_u, xmax_u, ymax_u = user_bbox

    # 1) Create a transformer if your data is in EPSG:3006 -> EPSG:4326
    transformer = pyproj.Transformer.from_crs(crs_from, "EPSG:4326", always_xy=True)

    # 2) We'll track lat/lon extents to center the map properly
    all_lons = []
    all_lats = []

    def add_bbox_coords(x1, y1, x2, y2):
        """
        Transform bounding box corners to WGS84 and record extent.

        Parameters
        ----------
        x1, y1, x2, y2 : float
            Corner coordinates in the source CRS (defaults to EPSG:3006).

        Returns
        -------
        tuple[float, float, float, float]
            Bounding box in longitude/latitude order suitable for Folium.
        """
        # Convert corners to lat/lon, add them to extents
        min_lon, min_lat = transformer.transform(x1, y1)
        max_lon, max_lat = transformer.transform(x2, y2)
        all_lons.extend([min_lon, max_lon])
        all_lats.extend([min_lat, max_lat])
        return (min_lon, min_lat, max_lon, max_lat)

    # 3) Convert user bbox to lat/lon & store extents
    user_min_lon, user_min_lat, user_max_lon, user_max_lat = add_bbox_coords(
        xmin_u, ymin_u, xmax_u, ymax_u
    )

    # 4) Convert each tile bbox to lat/lon & store extents
    converted_tiles = []
    for tile in tiles:
        txmin = tile["xmin"]
        tymin = tile["ymin"]
        txmax = tile["xmax"]
        tymax = tile["ymax"]
        tmin_lon, tmin_lat, tmax_lon, tmax_lat = add_bbox_coords(txmin, tymin, txmax, tymax)
        converted_tiles.append({
            "filename": tile["filename"],
            "min_lon": tmin_lon,
            "min_lat": tmin_lat,
            "max_lon": tmax_lon,
            "max_lat": tmax_lat
        })

    # 5) Determine center of the map
    if not all_lons or not all_lats:
        # fallback if empty
        map_center = [59.0, 18.0]  # some default
    else:
        center_lon = (min(all_lons) + max(all_lons)) / 2.0
        center_lat = (min(all_lats) + max(all_lats)) / 2.0
        map_center = [center_lat, center_lon]

    # 6) Create Folium map
    m = folium.Map(location=map_center, zoom_start=8)

    # 7) Create two FeatureGroups:
    user_fg = folium.FeatureGroup(name="User BBox", show=True)
    tile_fg = folium.FeatureGroup(name="Lidar Tiles", show=True)

    # 8) Add rectangle for user bbox (green) to user_fg
    folium.Rectangle(
        bounds=[(user_min_lat, user_min_lon), (user_max_lat, user_max_lon)],
        color="green",
        fill=False,
        tooltip="User BBox"
    ).add_to(user_fg)

    # 9) Add rectangles for each tile (blue) to tile_fg
    for t in converted_tiles:
        tooltip_text = f"{t['filename']} ({t['min_lon']:.5f},{t['min_lat']:.5f}) -> ({t['max_lon']:.5f},{t['max_lat']:.5f})"
        folium.Rectangle(
            bounds=[(t["min_lat"], t["min_lon"]), (t["max_lat"], t["max_lon"])],
            color="blue",
            fill=False,
            tooltip=tooltip_text
        ).add_to(tile_fg)

    # 10) Add the FeatureGroups to the map
    user_fg.add_to(m)
    tile_fg.add_to(m)

    # 11) Add a LayerControl so we can toggle them on/off
    folium.LayerControl(collapsed=False).add_to(m)

    # 12) Save
    m.save(out_html)
    debug(f"Map saved to {out_html}")


# ------------------------------------------------------------------------
# Async download with caching
# ------------------------------------------------------------------------
async def download_laz_file(session, base_url, filename, output_dir, semaphore):
    """
    Download a single .laz file asynchronously with aiohttp if not already cached.
    The endpoint is assumed to be: f"{base_url}/get/lidar/{filename}"
    We'll store the downloaded file in output_dir/filename.
    """
    url = f"{base_url}/get/lidar/{filename}"
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, filename)

    # 1) Check local cache
    if os.path.exists(out_path):
        info(f"File {filename} already in cache, skipping download.")
        return  # skip

    tmp_path = f"{out_path}.part"
    timeout = aiohttp.ClientTimeout(
        total=_DOWNLOAD_TOTAL_TIMEOUT_SECONDS,
        connect=_DOWNLOAD_CONNECT_TIMEOUT_SECONDS,
        sock_read=_DOWNLOAD_SOCK_READ_TIMEOUT_SECONDS,
    )

    async with semaphore:
        for attempt in range(1, _DOWNLOAD_MAX_ATTEMPTS + 1):
            try:
                info(
                    f"Downloading {filename} from {url} "
                    f"(attempt {attempt}/{_DOWNLOAD_MAX_ATTEMPTS})"
                )
                async with session.get(url, timeout=timeout) as resp:
                    if resp.status != 200:
                        retry_after = resp.headers.get("Retry-After")
                        raise LidarDownloadError(filename, resp.status, retry_after)

                    with open(tmp_path, "wb") as f:
                        async for chunk in resp.content.iter_chunked(1024 * 1024):
                            f.write(chunk)

                os.replace(tmp_path, out_path)
                info(f"Saved {filename} to {out_path}")
                return
            except (aiohttp.ClientError, asyncio.TimeoutError, OSError, RuntimeError) as exc:
                try:
                    os.remove(tmp_path)
                except OSError:
                    pass

                if attempt >= _DOWNLOAD_MAX_ATTEMPTS:
                    raise

                if isinstance(exc, LidarDownloadError) and exc.status == 429:
                    if exc.retry_after is not None:
                        try:
                            backoff = int(exc.retry_after)
                        except ValueError:
                            backoff = _DOWNLOAD_RATE_LIMIT_BACKOFF_SECONDS * attempt
                    else:
                        backoff = _DOWNLOAD_RATE_LIMIT_BACKOFF_SECONDS * attempt
                else:
                    backoff = _DOWNLOAD_RETRY_BACKOFF_SECONDS * attempt

                warning(
                    f"Lidar file download attempt {attempt}/{_DOWNLOAD_MAX_ATTEMPTS} "
                    f"failed for {filename}: {type(exc).__name__}: {exc}. "
                    f"Retrying in {backoff:.1f}s."
                )
                await asyncio.sleep(backoff)

async def download_all_lidar_files(base_url, filenames, output_dir="downloaded_laz"):
    """
    Given a list of filenames, downloads them all asynchronously from
    base_url/get/lidar/<filename> using aiohttp, skipping any local cache hits.
    """
    total_files = len(filenames)
    completed = [0]  # Use list to allow modification in nested function

    async def download_with_progress(session, base_url, filename, output_dir):
        await download_laz_file(session, base_url, filename, output_dir, semaphore)
        completed[0] += 1
        report_progress(
            current=completed[0],
            total=total_files,
            message=f"Downloading files ({completed[0]}/{total_files})..."
        )

    report_progress(percent=0, message=f"Downloading {total_files} files...")
    info(
        "Lidar download settings: "
        f"max_concurrency={_DOWNLOAD_MAX_CONCURRENCY}, "
        f"total_timeout={_DOWNLOAD_TOTAL_TIMEOUT_SECONDS}s, "
        f"connect_timeout={_DOWNLOAD_CONNECT_TIMEOUT_SECONDS}s, "
        f"sock_read_timeout={_DOWNLOAD_SOCK_READ_TIMEOUT_SECONDS}s, "
        f"max_attempts={_DOWNLOAD_MAX_ATTEMPTS}, "
        f"rate_limit_backoff={_DOWNLOAD_RATE_LIMIT_BACKOFF_SECONDS}s"
    )

    connector = aiohttp.TCPConnector(limit=_DOWNLOAD_MAX_CONCURRENCY)
    semaphore = asyncio.Semaphore(_DOWNLOAD_MAX_CONCURRENCY)
    async with aiohttp.ClientSession(connector=connector) as session:
        tasks = []
        for fname in filenames:
            tasks.append(download_with_progress(session, base_url, fname, output_dir))
        # Run all downloads concurrently
        await asyncio.gather(*tasks)

def run_download_files(base_url, filenames, output_dir="downloaded_laz"):
    """
    Entry point to run the async download with asyncio, skipping already cached files.
    """
    if not filenames:
        info("No files to download.")
        return
    debug(f"Downloading {len(filenames)} files in parallel (with cache check)...")
    asyncio.run(download_all_lidar_files(base_url, filenames, output_dir))
    info("All downloads finished.")


# ------------------------------------------------------------------------
# The single function that does everything for the user
# ------------------------------------------------------------------------
def download_lidar(user_bbox, session, buffer_val=0, base_url="http://127.0.0.1:8000",
                   output_map="client_map.html", output_dir= None):
    """
    1) POST the bounding box + buffer to the server -> get intersecting tiles
    2) Plot user bbox + tile bboxes in a Folium map
    3) Download .laz files in parallel, skipping any that exist locally
    """
    if output_dir is None: 
        cache_dir = user_cache_dir(appname="dtcc-data")
        os.makedirs(cache_dir, exist_ok=True)
        output_dir = os.path.join(cache_dir,'downloaded_laz')
    # A) Prepare endpoint
    endpoint_post = f"{base_url}/get_lidar"

    # B) Call the server
    try:
        response_data = post_lidar_request(
            endpoint_post,
            session,
            xmin=user_bbox[0],
            ymin=user_bbox[1],
            xmax=user_bbox[2],
            ymax=user_bbox[3],
            buffer_value=buffer_val
        )
    except Exception as e:
        warning(f"Error occurred: {e}")
        raise

    debug("Response from server:", response_data)

    # C) Plot bboxes
    returned_tiles = response_data["tiles"]
    output_map = os.path.join(cache_dir,output_map)
    plot_bboxes_folium(user_bbox, returned_tiles, out_html=output_map, crs_from="EPSG:3006")

    # D) Download files in parallel (with local cache)
    filenames_to_download = [tile["filename"] for tile in returned_tiles]
    run_download_files(base_url, filenames_to_download, output_dir=output_dir)
    return [os.path.join(output_dir, filename) for filename in filenames_to_download]


# ------------------------------------------------------------------------
# Example usage: only 3 lines needed
# ------------------------------------------------------------------------
#if __name__ == "__main__":
#    user_bbox = (267000, 6519000, 268000, 6521000)
#    buffer_val = 100
#    download_lidar(user_bbox, buffer_val)
