"""Report measured feature costs; research only, not a derivation of delta.

The former nearest-feature model charged entire edges for a short contact at
one endpoint. It therefore confused boundary sampling with narrow passages.
Its city delta recommendations are withdrawn. This report uses actual meshes
from the corrected precondition sweep; it makes no cross-family extrapolation
and cannot choose a city-wide cleaning floor.

  .venv/bin/python sandbox/cleaning_delta_budget_probe.py \
      --preconditions /tmp/preconditions/preconditions.json --output /tmp/costs.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def measured_costs(data):
    """Compare face counts within each family at the same requested spacing.

    A square hole's perimeter shrinks with its side, unlike the fixed facing
    walls. Keeping families separate avoids assigning a common contact length.
    Failed meshes stay visible; they are never interpreted as free geometry.
    """
    references = {}
    for row in data["rows"]:
        if row["separation"] == 1.0:
            for attempt in row["attempts"]:
                if attempt["outcome"] == "meshed":
                    references[row["family"], attempt["requested_spacing"]] = attempt[
                        "faces"
                    ]
    rows = []
    for row in data["rows"]:
        for attempt in row["attempts"]:
            base = references.get((row["family"], attempt["requested_spacing"]))
            faces = attempt.get("faces") if attempt["outcome"] == "meshed" else None
            rows.append(
                dict(
                    family=row["family"],
                    separation=row["separation"],
                    requested_spacing=attempt["requested_spacing"],
                    outcome=attempt["outcome"],
                    faces=faces,
                    reference_faces=base,
                    extra_faces=(
                        faces - base if faces is not None and base is not None else None
                    ),
                )
            )
    return rows


def main():
    import hashlib

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preconditions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    data = json.loads(args.preconditions.read_text())
    # Old measurements used a passage of width 1+s. Refuse that calibration.
    source = Path(__file__).with_name("cleaning_mesher_precondition_probe.py")
    expected = hashlib.sha256(source.read_bytes()).hexdigest()
    if (
        data.get("scripts", {}).get("sandbox/cleaning_mesher_precondition_probe.py")
        != expected
    ):
        parser.error(
            "preconditions do not match the current corrected probe; rerun the sweep"
        )
    result = {
        "scope": "Observed mesh costs only; no derived delta or city-cost prediction.",
        "preconditions": str(args.preconditions),
        "preconditions_sha256": hashlib.sha256(
            args.preconditions.read_bytes()
        ).hexdigest(),
        "rows": measured_costs(data),
    }
    with args.output.open("x") as handle:
        handle.write(json.dumps(result, indent=2) + "\n")
    print(f"Wrote {len(result['rows'])} measured costs to {args.output}")


if __name__ == "__main__":
    main()
