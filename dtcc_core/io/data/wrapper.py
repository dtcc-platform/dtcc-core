
#!/usr/bin/env python3
import requests
import getpass
import sys
import os
from pathlib import Path
from .overpass import get_roads_for_bbox, get_buildings_for_bbox
from .geopkg import CACHE_DIR as GPKG_CACHE_DIR, download_tiles
from .lidar import download_lidar
from dtcc_core import io
from dtcc_core.model import Bounds
from .logging import info, warning, debug, error

# We'll allow "lidar" or "roads" or "footprints" for data_type, and "dtcc" or "OSM" for provider.
valid_types = ["lidar", "roads", "footprints"]
valid_providers = ["dtcc", "OSM"]

# Env-overridable backend URLs (defaults preserve current behavior).
_DTCC_BASE = os.environ.get("DTCC_DATA_URL", "http://compute.dtcc.chalmers.se")
DTCC_LIDAR_URL = os.environ.get("DTCC_LIDAR_URL", f"{_DTCC_BASE}:8000")
DTCC_GPKG_URL  = os.environ.get("DTCC_GPKG_URL",  f"{_DTCC_BASE}:8001")

# We'll keep a single global SSH client in memory
SSH_CLIENT = None
SSH_CREDS = {
    "username": None,
    "password": None
}
sessions = []


def _bounds_overlap(lhs: Bounds, rhs: Bounds) -> bool:
    return not (
        lhs.xmax < rhs.xmin
        or lhs.xmin > rhs.xmax
        or lhs.ymax < rhs.ymin
        or lhs.ymin > rhs.ymax
    )


def _find_cached_footprint_files(bounds: Bounds) -> list[str]:
    cache_dir = Path(GPKG_CACHE_DIR) / "downloaded-gpkg"
    if not cache_dir.is_dir():
        return []

    matching_files: list[str] = []
    for path in sorted(cache_dir.glob("*.gpkg")):
        try:
            file_bounds = io.footprints.building_bounds(path)
        except Exception as exc:
            warning(f"Skipping cached footprint tile {path.name}: {exc}")
            continue
        if _bounds_overlap(file_bounds, bounds):
            matching_files.append(str(path))

    return matching_files


def _load_cached_footprints(bounds: Bounds):
    cached_files = _find_cached_footprint_files(bounds)
    if not cached_files:
        return None

    info(
        "Using %d cached footprint tile(s) from %s",
        len(cached_files),
        os.path.dirname(cached_files[0]),
    )
    buildings = io.load_footprints(cached_files, bounds=bounds)
    return buildings if buildings else None

def get_authenticated_session(base_url: str, username: str, password: str) -> requests.Session:
    """
    1. POST to /auth/token to obtain a bearer token.
    2. Create a requests.Session that automatically sends the token for future requests during runtime.
    """
    # 1) Obtain the token
    token_url = f"{base_url.rstrip('/')}/auth/token"
    payload = {"username": username, "password": password}

    response = requests.post(token_url, json=payload)
    if response.status_code != 200:
        error(f"Token request failed. Status code: {response.status_code}")
        return

    data = response.json()
    if "token" not in data:
        raise RuntimeError(f"No token found in response: {data}")

    token = data["token"]

    # 2) Create and return a Session with the token in headers
    session = requests.Session()
    session.headers.update({"Authorization": f"Bearer {token}"})
    return session

class SSHAuthenticationError(Exception):
    """Raised if SSH authentication fails."""
    pass

def _ssh_connect_if_needed():
    """
    Ensures we're authenticated via SSH to data.dtcc.chalmers.se.
    If not connected, prompts user for username/password, tries to connect.
    On success, we store the SSH client in memory for future calls.

    In non-interactive mode (e.g., when running as a server), credentials
    can be provided via DTCC_SSH_USERNAME and DTCC_SSH_PASSWORD environment
    variables.
    """
    global SSH_CLIENT, SSH_CREDS
    global sessions
    # If no credentials, prompt user
    if not sessions:
        # Check for environment variables first
        USERNAME = os.environ.get("DTCC_SSH_USERNAME")
        PASSWORD = os.environ.get("DTCC_SSH_PASSWORD")

        if not USERNAME or not PASSWORD:
            # Only prompt if in an interactive terminal
            if sys.stdin.isatty():
                info("SSH Authentication required for dtcc provider.")
                USERNAME = input("Enter SSH username: ")
                PASSWORD = getpass.getpass("Enter SSH password: ")
            else:
                warning("Non-interactive mode: SSH authentication skipped. "
                       "Set DTCC_SSH_USERNAME and DTCC_SSH_PASSWORD environment variables.")
                return None

        lidar_session = get_authenticated_session(DTCC_LIDAR_URL, USERNAME, PASSWORD)
        gpkg_session  = get_authenticated_session(DTCC_GPKG_URL,  USERNAME, PASSWORD)
        return lidar_session, gpkg_session
    return sessions

    # # Create a new SSH client
    # SSH_CLIENT = paramiko.SSHClient()
    # SSH_CLIENT.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    # try:
    #     SSH_CLIENT.connect(
    #         hostname="data.dtcc.chalmers.se",
    #         username=SSH_CREDS["username"],
    #         password=SSH_CREDS["password"]
    #     )
    # except paramiko.AuthenticationException as e:
    #     # If auth fails, raise an error and reset SSH_CLIENT
    #     SSH_CLIENT = None
    #     raise SSHAuthenticationError(f"SSH authentication failed: {e}")

    # print("SSH authenticated with data.dtcc.chalmers.se (no SFTP).")

def download_data(data_type: str, provider: str, bounds: Bounds, epsg = '3006', url = None):
    """
    A wrapper for downloading data, but with a dummy step for actual file transfer.
    If provider='dtcc', we do an SSH-based authentication check and then simulate a download.
    If provider='OSM', we just do a dummy download with no SSH.

    :param data_type: 'lidar' or 'roads' or 'footprints'
    :param provider: 'dtcc' or 'OSM'
    :return: dict with info about the (dummy) download
    """
    # Resolve per-service URLs: explicit `url` overrides env; otherwise use env-backed defaults.
    if url is not None:
        lidar_url, gpkg_url = f"{url}:8000", f"{url}:8001"
    else:
        lidar_url, gpkg_url = DTCC_LIDAR_URL, DTCC_GPKG_URL
    # Ensure user provided bounding box is a dtcc.Bounds object.
    if isinstance(bounds,(tuple | list)):
        bounds = Bounds(xmin=bounds[0],ymin=bounds[1],xmax=bounds[2],ymax=bounds[3])
    if not isinstance(bounds,Bounds):
        raise TypeError("user_bbox parameter must be of dtcc.Bounds type.")
    
    # user_bbox = user_bbox.tuple
    if not epsg == '3006':
        warning('Please enter the coordinates in EPSG:3006')
        return
    # Validate
    if data_type not in valid_types:
        raise ValueError(f"Invalid data_type '{data_type}'. Must be one of {valid_types}.")
    if provider not in valid_providers:
        raise ValueError(f"Invalid provider '{provider}'. Must be one of {valid_providers}.")

    if provider == "dtcc":

        global sessions
        session = requests.Session()
        if data_type == 'lidar':
            info('Starting the Lidar files download from dtcc source')
            files = download_lidar(bounds.tuple, session, base_url=lidar_url)
            if not files:
                raise RuntimeError("No lidar data available for the requested bounding box.")
            debug(files)
            pc = io.load_pointcloud(files,bounds=bounds)
            return pc
        elif data_type == 'footprints':
            cached_footprints = _load_cached_footprints(bounds)
            if cached_footprints is not None:
                return cached_footprints
            info("Starting the footprints download from dtcc source")
            files = download_tiles(bounds.tuple, session, server_url=gpkg_url)
            if not files:
                raise RuntimeError(
                    f"Footprint download failed for bounds {bounds.tuple}."
                )
            foots = io.load_footprints(files,bounds= bounds)
            return foots 
        else:
            error("Incorrect data type.")
        return

    else:  
        if data_type == 'footprints':
            info("Starting footprints files download from OSM source")
            gdf, filename = get_buildings_for_bbox(bounds.tuple)
            footprints = io.load_footprints(filename, bounds=bounds)
            return footprints
        elif data_type == 'roads':
            info('Start the roads files download from OSM source')
            gdf, filename = get_roads_for_bbox(bounds.tuple)
            roads = io.load_roadnetwork(filename)
            return roads
        else:
            error('Please enter a valid data type')
        return
   
def download_pointcloud(bounds: Bounds, provider = 'dtcc', epsg = '3006'):
    """
    Download a point cloud from the specified provider within the given bounds.

    Args:
        bounds (Bounds): The geographic bounds to download the point cloud data for.
        provider (str, optional): The data provider, defaults to 'dtcc'.
        epsg (str, optional): The EPSG code for the coordinate reference system, defaults to '3006'.

    Returns:
        Result of the download_data function call for 'lidar' data type if provider is 'dtcc'.

    Raises:
        Error if an invalid provider is specified.
    """

    if not provider or provider.lower() == 'dtcc':
        return download_data('lidar', 'dtcc', bounds, epsg=epsg)
    else:
        error("Please enter a valid provider")

def download_footprints(bounds: Bounds, provider = 'dtcc', epsg = '3006'):
    """
    Download building footprints from the specified provider within the given bounds.

    Args:
        bounds (Bounds): The geographic bounds to download the building footprints data for.
        provider (str, optional): The data provider, defaults to 'dtcc'.
        epsg (str, optional): The EPSG code for the coordinate reference system, defaults to '3006'.

    Returns:
        Result of the download_data function call for 'footprints' data type if provider is 'dtcc'.
        Result of the download_data function call for 'footprints' data type if provider is 'OSM'.

    Raises:
        Error if an invalid provider is specified.
    """
    if not provider or provider.lower() == 'dtcc':
        return download_data('footprints', 'dtcc', bounds, epsg=epsg)
    elif provider.upper() == 'OSM':
        return download_data('footprints', "OSM", bounds, epsg = epsg)
    else:
        error("Please enter a valid provider")

def download_roadnetwork(bounds: Bounds, provider = 'dtcc', epsg='3006'):
    """
    Download road network data from the specified provider within the given bounds.

    Args:
        bounds (Bounds): The geographic bounds to download the road network data for.
        provider (str, optional): The data provider, defaults to 'dtcc'.
        epsg (str, optional): The EPSG code for the coordinate reference system, defaults to '3006'.

    Returns:
        Result of the download_data function call for 'roads' data type if provider is 'OSM'.

    Raises:
        Error if an invalid provider is specified.
    """
    if provider and provider.upper() == 'OSM':
        download_data('roads', "OSM", bounds, epsg=epsg)
    else:
        error("Please enter a valid provider")
