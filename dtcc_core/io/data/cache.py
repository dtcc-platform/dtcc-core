from pathlib import Path
from platformdirs import user_cache_dir
import shutil

cache_dir = Path(user_cache_dir("dtcc-data"))  # Replace with your app name


def empty_cache(cache_type = None):
    """
    Clears the contents of the cache directory.

    Parameters
    ----------
    cache_type : optional
        Currently unused, reserved for future functionality to clear specific
        types of cache.
    """

    if cache_dir.exists():
        for item in cache_dir.iterdir():
            if item.is_file() or item.is_symlink():
                item.unlink()  
            elif item.is_dir():
                shutil.rmtree(item)  
