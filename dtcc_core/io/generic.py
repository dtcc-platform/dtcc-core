# Copyright(C) 2023 Anders Logg
# Licensed under the MIT License

import pathlib
from .logging import info, warning, error


def save(object, path, name, formats, *args, **kwargs):
    """
    Save an object using a registered format handler.

    Parameters
    ----------
    object : Any
        Object to serialize.
    path : str or pathlib.Path
        Output path including extension.
    name : str
        Logical resource name (used in messages).
    formats : dict
        Mapping of supported types to extension-handler dictionaries.
    *args, **kwargs :
        Passed through to the selected saver.

    Raises
    ------
    SystemExit
        If the object type or file extension is unsupported (via ``error``).
    """
    if not type(object) in formats:
        error(f'Unable to save {name}; type "{type(object)}" not supported')
    path = pathlib.Path(path)

    path_suffix = path.suffix
    two_level_suffix = ''.join(path.suffixes[-2:]) if len(path.suffixes) >= 2 else path_suffix

    if path_suffix not in formats[type(object)] and two_level_suffix not in formats[type(object)]:
        error(
            f"Unable to save {name} ({type(object).__name__}); format {path.suffix} not supported"
        )

    info(f"Saving {name} ({type(object).__name__}) to {path}")
    if two_level_suffix in formats[type(object)]:
        saver = formats[type(object)][two_level_suffix]
    elif path_suffix in formats[type(object)]:
        saver = formats[type(object)][path_suffix]
    else:
        error(f"Unable to save {name}; format {path.suffix} not supported")

    saver(object, path, *args, **kwargs)


def load(path, name, type, formats, *args, **kwargs):
    """
    Load an object using a registered format handler.

    Parameters
    ----------
    path : str or pathlib.Path or list[str | pathlib.Path]
        Path(s) to read from.
    name : str
        Logical resource name (used in messages).
    type : type
        Expected target class key in ``formats``.
    formats : dict
        Mapping of supported types to extension-handler dictionaries.
    *args, **kwargs :
        Passed through to the selected loader.

    Returns
    -------
    Any
        Object returned by the loader.

    Raises
    ------
    SystemExit
        If the type or file extension is unsupported (via ``error``).
    """
    if not type in formats:
        error(f'Unable to load {name}; type "{type.__name__}" not supported')
    if isinstance(path, (list, tuple)):
        path = [pathlib.Path(p) for p in path]
        path_suffix = path[0].suffix
        two_level_suffix = ''.join(path[0].suffixes[-2:]) # e.g. ['.json.zip]

    else:
        path = pathlib.Path(path)
        path_suffix = path.suffix
        two_level_suffix = ''.join(path.suffixes[-2:]) # e.g. ['.json.zip]
    if path_suffix not in formats[type] and two_level_suffix not in formats[type]:
        error(f"Unable to load {name}; format {path.suffix} not supported")
    info(f"Loading {name} ({type.__name__}) from {path}")
    if two_level_suffix in formats[type]:
        loader = formats[type][two_level_suffix]
    elif path_suffix in formats[type]:
        loader = formats[type][path_suffix]
    else:
        error(f"Unable to load {name}; format {path.suffix} not supported")
    return loader(path, *args, **kwargs)


def list_io(name, load_formats, save_formats):
    """
    Summarize supported load/save formats for a resource.

    Parameters
    ----------
    name : str
        Logical resource name.
    load_formats : dict
        Supported loaders keyed by type.
    save_formats : dict
        Supported savers keyed by type.

    Returns
    -------
    dict
        Mapping with ``load_formats`` and ``save_formats`` keys.
    """
    return {
        "load_formats": load_formats.keys(),
        "save_formats": save_formats.keys(),
    }


def print_io(name, load_formats, save_formats):
    """
    Print supported formats for loading and saving a resource.

    Parameters
    ----------
    name : str
        Logical resource name.
    load_formats : dict
        Supported loaders keyed by type.
    save_formats : dict
        Supported savers keyed by type.

    Returns
    -------
    None
        Output is printed to stdout.
    """
    print(f"load_{name}() supports the following data types and formats:")
    print("")
    N = max([len(t.__name__) for t in load_formats])
    for t in load_formats:
        n = N - len(t.__name__)
        formats = ", ".join([f for f in load_formats[t]])
        print(f"  {t.__name__}: {n*' '}{formats}")
    print("")
    N = max([len(t.__name__) for t in save_formats])
    print(f"save_{name}() supports the following data types and formats:")
    print("")
    for t in save_formats:
        n = N - len(t.__name__)
        formats = ", ".join([f for f in save_formats[t]])
        print(f"  {t.__name__}: {n*' '} {formats}")
