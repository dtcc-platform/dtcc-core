"""Capture lossless focused fixtures from the frozen 2026-09-19 survey.

The full benchmark artifacts stay below ``benchmarks/runs`` (which is ignored).
This script extracts only the interaction groups named by the active footprint
cleaning plan into tracked test data, retaining world-coordinate WKB and source
indices.  It deliberately does not run either cleaner.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from dtcc_core.builder.cleaning.contract import _interpret_input
from dtcc_core.io.model import load_model
from dtcc_core.model import City, GeometryType
from sandbox.cleaning_graph_prototype import interaction_groups, polygon_parts


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "benchmarks/runs/2026-09-19_full-survey-audit"
OUTPUT = ROOT / "tests/data/cleaning/full-survey-focused-groups.json"

# The legacy generation containing the valid saved input, and the exact groups
# identified by the audit.  ``None`` retains every group for the two tail cases
# whose performance boundary is the whole tile rather than one correctness bug.
CASES = {
    "gothenburg:006": ("legacy2", [84]),
    "helsingborg:007": ("legacy2", [42]),
    "helsingborg:015": ("legacy2", [2]),
    "linkoping:066": ("legacy2", [45]),
    "lund:016": ("legacy3", [8]),
    "norrkoping:003": ("legacy3", None),
    "norrkoping:039": ("legacy3", [23]),
    "stockholm:083": ("legacy2", [7, 8]),
    "vasteras:057": ("legacy3", None),
}


def _source_path(generation: str, city: str, tile: str) -> Path:
    task = f"survey_city_flat_mesh_city_grid_{city}_{tile}_baseline"
    return AUDIT / generation / city / "tasks" / task / "cleaning/raw.dtcc"


def _atoms(city: City):
    values = []
    for source_index, building in enumerate(city.buildings):
        geometry = building.flatten_geometry(GeometryType.LOD0)
        if geometry is None:
            continue
        polygon = geometry.to_polygon(simplify=0.0)
        if polygon is None or polygon.is_empty:
            continue
        interpreted, _ = _interpret_input([polygon])
        values.extend((part, source_index) for part in polygon_parts(interpreted))
    return values


def capture() -> dict:
    cases = []
    for case_id, (generation, selected_groups) in CASES.items():
        city_name, tile = case_id.split(":")
        source = _source_path(generation, city_name, tile)
        city = load_model(source, expected_type=City)
        atoms = _atoms(city)
        groups = interaction_groups(
            [polygon for polygon, _ in atoms],
            1.0000005,  # 2 * epsilon + delta + the checker's delta tolerance
        )
        indices = range(len(groups)) if selected_groups is None else selected_groups
        captured_groups = []
        for group_index in indices:
            members = groups[group_index]
            captured_groups.append(
                {
                    "group": group_index,
                    "polygon_wkb_hex": [atoms[index][0].wkb_hex for index in members],
                    "source_indices": [atoms[index][1] for index in members],
                }
            )
        cases.append(
            {
                "case_id": f"city_grid:{case_id}",
                "task_id": (
                    f"survey:city_flat_mesh:city_grid:{city_name}:{tile}:baseline"
                ),
                "source_artifact": str(source.relative_to(ROOT)),
                "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                "source_building_count": len(city.buildings),
                "interaction_group_count": len(groups),
                "groups": captured_groups,
            }
        )
    return {
        "version": 1,
        "crs": "EPSG:3006",
        "coordinate_units": "metres",
        "delta": 0.5,
        "epsilon": 0.25,
        "group_distance": 1.0000005,
        "provenance": (
            "Exact interpreted polygon atoms from the valid saved raw artifacts "
            "of the 2026-09-19 1,000-tile paired audit."
        ),
        "cases": cases,
    }


def main() -> None:
    if not AUDIT.is_dir():
        raise SystemExit(f"missing frozen audit directory: {AUDIT}")
    payload = capture()
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {len(payload['cases'])} cases to {OUTPUT}")


if __name__ == "__main__":
    main()
