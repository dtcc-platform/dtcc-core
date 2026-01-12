#!/usr/bin/env python3
"""
Script to find all functions and methods missing docstrings in dtcc_core directory.
"""

import os
import re
import ast
from typing import List, Tuple, Dict
from pathlib import Path


def has_docstring(node):
    """Check if a function/method has a docstring."""
    if not node.body:
        return False

    first_stmt = node.body[0]
    return (
        isinstance(first_stmt, ast.Expr)
        and isinstance(first_stmt.value, ast.Constant)
        and isinstance(first_stmt.value.value, str)
    )


def analyze_file(file_path: str) -> List[Tuple[str, str, int]]:
    """Analyze a Python file and return functions/methods missing docstrings."""
    missing_docstrings = []

    try:
        with open(file_path, "r", encoding="utf-8") as f:
            content = f.read()

        # Parse the AST
        tree = ast.parse(content)

        # Find all function and method definitions
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef):
                # Skip private methods starting with underscore (except __init__ and other special methods)
                if node.name.startswith("_") and not node.name.startswith("__"):
                    continue

                if not has_docstring(node):
                    missing_docstrings.append((file_path, node.name, node.lineno))

            elif isinstance(node, ast.AsyncFunctionDef):
                # Skip private async methods
                if node.name.startswith("_") and not node.name.startswith("__"):
                    continue

                if not has_docstring(node):
                    missing_docstrings.append(
                        (file_path, f"async {node.name}", node.lineno)
                    )

    except Exception as e:
        print(f"Error analyzing {file_path}: {e}")

    return missing_docstrings


def find_missing_docstrings(root_dir: str) -> List[Tuple[str, str, int]]:
    """Find all functions/methods missing docstrings in the given directory."""
    all_missing = []

    for root, dirs, files in os.walk(root_dir):
        # Skip external dependencies and cache directories
        if "external" in root or "__pycache__" in root:
            continue
        for file in files:
            if file.endswith(".py"):
                file_path = os.path.join(root, file)
                missing = analyze_file(file_path)
                all_missing.extend(missing)

    return all_missing


def main():
    root_dir = Path(__file__).parent.parent  # Adjust to your dtcc_core directory
    src_dir = root_dir / "dtcc_core"
    missing_docstrings = find_missing_docstrings(src_dir)

    if missing_docstrings:
        with open(root_dir / "MISSING_DOCSTRINGS.md", "w", encoding="utf-8") as f:

            print(
                f"Found {len(missing_docstrings)} functions/methods missing docstrings:\n"
            )
            f.write(
                f"# Missing Docstrings\n\n"
            )
            f.write(
                f"Found {len(missing_docstrings)} functions/methods missing docstrings.\n\n"
            )
            f.write(
                f"Run `python scripts/find_missing_docstrings.py` to regenerate this file.\n\n"
            )

            # Group by file for better readability
            by_file = {}
            for file_path, func_name, line_no in missing_docstrings:
                # Convert absolute path to relative path from repo root
                rel_path = os.path.relpath(file_path, root_dir)
                if rel_path not in by_file:
                    by_file[rel_path] = []
                by_file[rel_path].append((func_name, line_no))

            for file_path in sorted(by_file.keys()):
                print(f"File: {file_path}")
                f.write(f"## {file_path}\n\n")
                for func_name, line_no in sorted(
                    by_file[file_path], key=lambda x: x[1]
                ):
                    print(f"  - {func_name} (line {line_no})")
                    f.write(f"- `{func_name}` (line {line_no})\n")
                print()
                f.write("\n")
    else:
        print("No functions/methods missing docstrings found!")


if __name__ == "__main__":
    main()
