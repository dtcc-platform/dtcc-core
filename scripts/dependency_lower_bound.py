#!/usr/bin/env python3
"""
Print the lowest version pyproject.toml allows for one runtime dependency,
taken from its ">=" or "==" specifier.

Usage:
  python scripts/dependency_lower_bound.py aiohttp

CI installs that version to test against it, so the check follows the bound
in pyproject.toml instead of repeating it.

Exit codes:
  0: Version printed
  2: The dependency is missing or has no lower bound
"""

from __future__ import annotations

import sys
import tomllib
from pathlib import Path

from packaging.requirements import Requirement
from packaging.utils import canonicalize_name

PYPROJECT = Path(__file__).resolve().parent.parent / "pyproject.toml"


def lower_bound(name: str) -> str | None:
    project = tomllib.loads(PYPROJECT.read_text(encoding="utf-8"))["project"]
    for line in project["dependencies"]:
        requirement = Requirement(line)
        if canonicalize_name(requirement.name) != canonicalize_name(name):
            continue
        for specifier in requirement.specifier:
            if specifier.operator in {">=", "=="}:
                return specifier.version
    return None


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    version = lower_bound(sys.argv[1])
    if version is None:
        print(f"{sys.argv[1]} has no lower bound in {PYPROJECT}", file=sys.stderr)
        return 2
    print(version)
    return 0


if __name__ == "__main__":
    sys.exit(main())
