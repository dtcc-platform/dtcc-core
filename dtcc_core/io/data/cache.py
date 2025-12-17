from pathlib import Path
from platformdirs import user_cache_dir
import shutil

cache_dir = Path(user_cache_dir("dtcc-data"))  # Replace with your app name


def _is_within_cache(path: Path, base: Path) -> bool:
    """
    Ensure the resolved path stays within the cache root to avoid following
    symlinks outside of the cache directory.
    """
    try:
        path.resolve(strict=False).relative_to(base)
        return True
    except (RuntimeError, ValueError):
        return False


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
        cache_root = cache_dir.resolve()
        for item in cache_root.iterdir():
            # Skip anything resolving outside the cache dir
            if not _is_within_cache(item, cache_root):
                continue
            if item.is_symlink():
                item.unlink()
            elif item.is_dir(follow_symlinks=False):
                shutil.rmtree(item)
            elif item.is_file():
                item.unlink()
