
#!/usr/bin/env python3
import requests
import getpass
from dtcc_data.overpass import get_roads_for_bbox, get_buildings_for_bbox
from dtcc_data.geopkg import download_tiles
from .geopkg import download_tiles_dataset
from .lidar import download_lidar
from dtcc_core import io
from dtcc_core.model import Bounds
from .dtcc_logging import info, warning, debug, error
from .auth import create_authenticated_session
import os
# We'll allow "lidar" or "roads" or "footprints" for data_type, and "dtcc" or "OSM" for provider.
valid_types = ["lidar", "roads", "footprints"]
valid_providers = ["dtcc", "OSM"]

BASE_URL = os.getenv("BASE_URL",'http://127.0.0.1:8002')

# We'll keep a single global SSH client in memory
SSH_CLIENT = None
SSH_CREDS = {
    "username": None,
    "password": None
}
def download_data(data_type: str, provider: str, bounds: Bounds, epsg = '3006', url = BASE_URL):
    """
    A wrapper for downloading data, but with a dummy step for actual file transfer.
    If provider='dtcc', we do an SSH-based authentication check and then simulate a download.
    If provider='OSM', we just do a dummy download with no SSH.

    :param data_type: 'lidar' or 'roads' or 'footprints'
    :param provider: 'dtcc' or 'OSM'
    :return: dict with info about the (dummy) download
    """
    
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

        session = create_authenticated_session()
        if not session:
            return
        if data_type == 'lidar':
            info('Starting the Lidar files download from dtcc source')
            files = download_lidar(bounds.tuple, session, base_url=f'{url}')
            debug(files)
            pc = io.load_pointcloud(files,bounds=bounds)
            return pc
        elif data_type == 'footprints':
            info("Starting the footprints download from dtcc source")
            files = download_tiles(bounds.tuple, session, server_url=f"{url}")
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
    if not provider or provider.lower() == 'dtcc':
        return download_data('lidar', 'dtcc', bounds, epsg=epsg)
    else:
        error("Please enter a valid provider")

def download_footprints(bounds: Bounds, provider = 'dtcc', epsg = '3006'):
    if not provider or provider.lower() == 'dtcc':
        return download_data('footprints', 'dtcc', bounds, epsg=epsg)
    elif provider.upper() == 'OSM':
        return download_data('footprints', "OSM", bounds, epsg = epsg)
    else:
        error("Please enter a valid provider")

def download_roadnetwork(bounds: Bounds, provider = 'dtcc', epsg='3006'):
    if provider and provider.upper() == 'OSM':
        download_data('roads', "OSM", bounds, epsg=epsg)
    else:
        error("Please enter a valid provider")


def download_footprints_dataset(bounds: Bounds, dataset: str, provider = 'dtcc', epsg = '3006', url = BASE_URL):
    if isinstance(bounds,(tuple | list)):
        bounds = Bounds(xmin=bounds[0],ymin=bounds[1],xmax=bounds[2],ymax=bounds[3])
    if not isinstance(bounds,Bounds):
        raise TypeError("bounds must be of dtcc.Bounds type.")
    if epsg != '3006':
        warning('Please enter the coordinates in EPSG:3006')
        return
    if provider.lower() != 'dtcc':
        error("Only 'dtcc' provider supported for dataset-aware footprints")
        return
    session = create_authenticated_session()
    info(f"Starting footprints download for dataset '{dataset}' from dtcc source")
    files = download_tiles_dataset(bounds.tuple, session, dataset=dataset, server_url=f"{url}")
    if not files:
        return None
    foots = io.load_footprints(files, bounds=bounds)
    return foots
