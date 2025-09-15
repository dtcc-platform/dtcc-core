#!/usr/bin/env python3
"""
Check that every public API function exported via __all__ in the package
has at least one executed line in coverage.

Usage:
  python scripts/check_public_api_calls.py \
      --package dtcc_core \
      --coverage-file tests/coverage.json

Exit codes:
  0: All public functions covered by at least one executed line
  1: One or more public functions missed
  2: Invalid input or unexpected error
"""

from __future__ import annotations

import argparse
import importlib
import inspect
import json
import os
import pkgutil
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple


@dataclass
class FuncInfo:
    module: str
    name: str
    qualname: str
    file: str
    start: int
    end: int

    @property
    def key(self) -> Tuple[str, int, int, str]:
        return (self.file, self.start, self.end, self.qualname)


def iter_modules(package_name: str) -> Iterable[str]:
    pkg = importlib.import_module(package_name)
    # Include the package itself
    yield package_name
    if hasattr(pkg, "__path__"):
        for m in pkgutil.walk_packages(pkg.__path__, pkg.__name__ + "."):
            yield m.name


def get_public_functions_from_loaded_module(module) -> List[FuncInfo]:
    """Collect public functions from a loaded module (uses its __all__)."""
    modname = module.__name__
    names = getattr(module, "__all__", None)
    if not names:
        return []

    found: List[FuncInfo] = []
    for name in names:
        try:
            obj = getattr(module, name)
        except Exception:
            # Name not present as attribute in module (may be a submodule hint or typo);
            # skip gracefully, since we walk all submodules separately.
            continue

        # Only consider plain functions for this check
        if not inspect.isfunction(obj) and not inspect.isbuiltin(obj):
            continue

        # Try to resolve source file and lines
        try:
            src_file = inspect.getsourcefile(obj) or inspect.getfile(obj)
        except TypeError:
            # Builtin or C-accelerated without Python source
            src_file = None
        if not src_file:
            # No Python source to attribute coverage; treat as miss
            # but capture as best-effort with synthetic range
            found.append(
                FuncInfo(
                    module=modname,
                    name=name,
                    qualname=getattr(obj, "__qualname__", name),
                    file="<built-in>",
                    start=0,
                    end=0,
                )
            )
            continue

        try:
            lines, start = inspect.getsourcelines(obj)
        except (OSError, IOError):
            # Could not get source lines; mark as unknown range
            lines, start = [], 0

        end = start + max(len(lines), 1) - 1
        found.append(
            FuncInfo(
                module=modname,
                name=name,
                qualname=getattr(obj, "__qualname__", name),
                file=os.path.realpath(src_file),
                start=int(start),
                end=int(end),
            )
        )
    return found


def load_coverage_executed_lines(coverage_json_path: str) -> Dict[str, Set[int]]:
    with open(coverage_json_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    files = data.get("files") or {}
    executed: Dict[str, Set[int]] = {}
    for path, entry in files.items():
        # Normalize to real absolute path
        try:
            real = os.path.realpath(path)
        except Exception:
            real = path
        lines = set(entry.get("executed_lines") or [])
        executed[real] = lines
    return executed


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Check public API function calls using coverage")
    parser.add_argument("--package", required=True, help="Root package to scan, e.g., dtcc_core")
    parser.add_argument("--coverage-file", required=True, help="Path to coverage JSON file")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail if any submodule under the package fails to import",
    )
    args = parser.parse_args(argv)

    # Discover public functions
    func_map: Dict[Tuple[str, int, int, str], FuncInfo] = {}
    total_modules = 0
    import_errors: List[Tuple[str, str]] = []
    loaded_modules: List[object] = []
    for modname in iter_modules(args.package):
        total_modules += 1
        try:
            module = importlib.import_module(modname)
            loaded_modules.append(module)
        except Exception as e:
            import_errors.append((modname, f"{e}"))
            continue

    if args.strict and import_errors:
        print("ERROR: One or more submodules failed to import in strict mode:")
        for name, msg in import_errors:
            print(f" - {name}: {msg}")
        return 2

    for module in loaded_modules:
        funcs = get_public_functions_from_loaded_module(module)
        for fi in funcs:
            func_map[fi.key] = fi

    functions: List[FuncInfo] = list(func_map.values())

    # Load coverage executed lines
    try:
        executed_by_file = load_coverage_executed_lines(args.coverage_file)
    except Exception as e:
        print(f"ERROR: Failed to read coverage file '{args.coverage_file}': {e}")
        return 2

    # Evaluate coverage for functions
    misses: List[FuncInfo] = []
    covered_count = 0
    for fi in sorted(functions, key=lambda x: (x.file, x.start, x.name)):
        if fi.file == "<built-in>":
            # Consider built-ins as missed (no python source to track)
            misses.append(fi)
            continue

        executed = executed_by_file.get(fi.file)
        if not executed:
            # No executed lines recorded for file -> missed
            misses.append(fi)
            continue

        # Check any executed line within [start, end]
        if any(fi.start <= ln <= fi.end for ln in executed):
            covered_count += 1
        else:
            misses.append(fi)

    total = len(functions)
    print(f"Scanned modules: {total_modules}")
    print(f"Public functions found: {total}")
    print(f"Functions covered: {covered_count}")
    print(f"Functions missed: {len(misses)}")

    if misses:
        print("\nMissed public functions (module:name @ file:start-end):")
        for fi in misses:
            loc = f"{fi.file}:{fi.start}-{fi.end}" if fi.file else "<unknown>"
            print(f" - {fi.module}:{fi.name} ({fi.qualname}) @ {loc}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
