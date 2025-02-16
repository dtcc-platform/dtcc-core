# Copyright(C) 2023 Anders Logg
# Licensed under the MIT License

import pathlib
from .logging import info, warning, error


def save(object, path, name, formats, *args, **kwargs):
    if not type(object) in formats:
        error(f'Unable to save {name}; type "{type(object)}" not supported')
    path = pathlib.Path(path)
    if path.suffix not in formats[type(object)]:
        error(
            f"Unable to save {name} ({type(object).__name__}); format {path.suffix} not supported"
        )
    info(f"Saving {name} ({type(object).__name__}) to {path}")
    formats[type(object)][path.suffix](object, path, *args, **kwargs)


def load(path, name, type, formats, *args, **kwargs):
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
    return {
        "load_formats": load_formats.keys(),
        "save_formats": save_formats.keys(),
    }


def print_io(name, load_formats, save_formats):
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
