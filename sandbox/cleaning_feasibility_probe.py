"""Feasibility and baseline measurement; never called by production cleaning.

Milestones A and B of
`.agent/plans/2026-09-18-footprint-cleaning-feasibility-and-architecture-decision.md`.
A disposable global construction alters the interpreted occupancy by a swept
parameter, and the independent checker alone decides what conforms. Nothing
here is a candidate cleaner. One fixed candidate ladder is applied to every
group, so no result is per-city tuned. A sweep that finds nothing reports
"none found"; it never asserts that nothing exists.

Milestone C re-runs the recorded construction unchanged, only to add the
per-group canonical vertex count, accepted edits and time that its manifests
never stored. It imports `sandbox/cleaning_patch_probe.py` and does not edit it.

  .venv/bin/python sandbox/cleaning_feasibility_probe.py --output /tmp/feasibility \
      --run benchmarks/runs/2026-09-17_100049_quick
  .venv/bin/python sandbox/cleaning_feasibility_probe.py --output /tmp/feasibility \
      --run benchmarks/runs/2026-09-17_100049_quick --scaling

Run the two on an otherwise idle machine and never at the same time: both
report seconds.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from shapely import STRtree, distance, linestrings, make_valid, points, set_precision
from shapely import simplify as shapely_simplify
from shapely.errors import GEOSException
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dtcc_core.builder.cleaning.contract import (  # noqa: E402
    FidelityBudget,
    check_cleaning_contract,
    _SEPARATION_RELATIVE_TOLERANCE,
    _canonical_graph,
    _interpret_input,
    admissibility,
)
from sandbox.cleaning_graph_prototype import polygon_parts  # noqa: E402

# Zero comes first so an already admissible group is witnessed at no alteration
# at all. The ladder brackets the answer; a short bisection then narrows it.
EPSILON_LADDER = (
    0.0,
    0.005,
    0.01,
    0.02,
    0.04,
    0.0625,
    0.09,
    0.125,
    0.175,
    0.25,
    0.35,
    0.5,
    0.75,
    1.0,
    1.5,
    2.0,
    3.0,
    5.0,
)
EPSILON_REFINEMENTS = 5

# The declared control of milestone B: topology-preserving simplification alone,
# swept over one tolerance ladder and canonicalized. Kept separate from the
# wider milestone A ladder so the baseline is never credited with snapping.
BASELINE_TOLERANCES = (0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0)


def rounded(value, digits=9):
    """Round floats before writing, so a saved result is readable and diffable."""
    if isinstance(value, float):
        return round(value, digits)
    if isinstance(value, dict):
        return {k: rounded(v, digits) for k, v in value.items()}
    if isinstance(value, list):
        return [rounded(v, digits) for v in value]
    return value


def write_json(destination, payload):
    destination.write_text(json.dumps(rounded(payload), indent=1) + "\n")


def unresolved_groups(city):
    """Groups the recorded neighbouring-clearance construction left unresolved."""
    manifest = json.loads(
        (ROOT / "docs/design/footprint-cleaning-patch-results.json").read_text()
    )
    return {
        entry["group"]
        for entry in manifest["neighbor_clearance_experiment"]["unresolved"]
        if entry["city"] == city
    }


def provenance(run=None):
    """Record what produced a number, in the form the earlier manifests use."""
    import hashlib
    import platform

    import numpy
    import shapely
    import scipy

    def digest(path):
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()

    record = {
        "scripts": {
            name: digest(ROOT / name)
            for name in (
                "sandbox/cleaning_feasibility_probe.py",
                "dtcc_core/builder/cleaning/contract.py",
                "sandbox/cleaning_graph_prototype.py",
                "sandbox/cleaning_pipeline_probe.py",
            )
        },
        "versions": {
            "shapely": shapely.__version__,
            "geos": ".".join(str(n) for n in shapely.geos_version),
            "numpy": numpy.__version__,
            "scipy": scipy.__version__,
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "measurement_environment": (
            "Linux aarch64 virtual machine on the author's workstation, not the "
            "repository's macOS .venv. Library and mesher revisions match; "
            "absolute timings are therefore not comparable with earlier "
            "manifests, and a calibration sample is recorded separately."
        ),
    }
    if run is not None:
        run = Path(run)
        results = json.loads((run / "results.json").read_text())["results"]
        record["source_artifacts"] = [
            {
                "city": r["case"]["city"],
                "path": str(run / r["artifacts"]["cleaning_input"]["path"]),
                "sha256": digest(run / r["artifacts"]["cleaning_input"]["path"]),
                "raw_sha256": digest(
                    (run / r["artifacts"]["cleaning_input"]["path"]).parent / "raw.dtcc"
                ),
            }
            for r in results
            if r["dataset"] == "city_flat_mesh"
        ]
    return record


def candidate_ladder():
    """One fixed, ordered ladder of global alterations, applied to every group.

    Snap rounding at a pixel width and topology-preserving simplification at a
    tolerance are the two operations the plan names. Morphological closing at a
    radius is added because neither of those ever fills a sub-delta gap between
    two distinct components, which is the defect the corpus is full of.
    Compositions are included because no single operation removes both dense
    boundary sampling and sub-delta gaps. The ladder is fixed here, before any
    corpus result is seen, and is applied unchanged to every group.
    """
    ladder = [()]
    ladder += [(("snap", w),) for w in (0.01, 0.05, 0.1, 0.25, 0.5, 1.0)]
    ladder += [(("simplify", t),) for t in BASELINE_TOLERANCES]
    ladder += [(("close", r),) for r in (0.0625, 0.125, 0.25, 0.5)]
    ladder += [
        (("close", r), ("simplify", t))
        for r in (0.0625, 0.125, 0.25, 0.5)
        for t in (0.01, 0.1, 0.5)
    ]
    ladder += [
        (("close", r), ("snap", w))
        for r in (0.125, 0.25, 0.5)
        for w in (0.125, 0.25, 0.5)
    ]
    ladder += [
        (("close", r), ("simplify", 0.05), ("snap", w))
        for r in (0.125, 0.25)
        for w in (0.125, 0.25, 0.5)
    ]
    ladder += [
        (("simplify", t), ("snap", w)) for t in (0.05, 0.25) for w in (0.125, 0.25, 0.5)
    ]
    ladder += [
        (("snap", w), ("simplify", t)) for w in (0.125, 0.25, 0.5) for t in (0.05, 0.25)
    ]
    return ladder


def candidate_name(operations):
    return "+".join(f"{kind}:{value:g}" for kind, value in operations) or "identity"


def reference_output(occupied, operations):
    """Alter the interpreted occupancy globally, then canonicalize.

    Disposable by intent: every operation is single shot and terminates, none
    carries a defect model, and none is applied case by case. Mitre joins keep
    a closed geometry polygonal; round joins would sample every corner far
    below delta and make the result inadmissible for a reason the experiment
    introduced. Keeping the polygonal parts and re-noding them is the
    canonicalization step, since a simplified boundary can fold onto itself and
    the contract's subdivision is built from an interior-disjoint coverage.
    """
    geometry = occupied
    for kind, value in operations:
        if kind == "close":
            geometry = geometry.buffer(value, join_style="mitre").buffer(
                -value, join_style="mitre"
            )
        elif kind == "snap":
            geometry = set_precision(geometry, value, mode="valid_output")
        elif kind == "simplify":
            geometry = shapely_simplify(geometry, value, preserve_topology=True)
        else:
            raise ValueError(f"unknown reference operation {kind}")
    pieces = polygon_parts(make_valid(geometry))
    return polygon_parts(unary_union(pieces)) if pieces else []


def defect_classes(polygons, delta):
    """Split the checker's subscale pairs by kind and by boundary curve.

    Descriptive only. The split repeats the checker's own queries on the
    checker's own canonical graph and its totals are asserted against
    `admissibility` by the caller, so this breakdown cannot quietly disagree
    with the conformance decision. Two features on the same connected component
    of the graph lie on one boundary curve, which is dense sampling or a thin
    neck; different components are two distinct curves too close together.
    """
    counts = {
        "vertex_vertex_same_curve": 0,
        "vertex_vertex_other_curve": 0,
        "vertex_edge_same_curve": 0,
        "vertex_edge_other_curve": 0,
    }
    vertices, edges = _canonical_graph([p.boundary for p in polygons])
    if not vertices:
        return counts
    vertex_ids = {xy: i for i, xy in enumerate(vertices)}
    component = list(range(len(vertices)))

    def root(i):
        while component[i] != i:
            component[i] = component[component[i]]
            i = component[i]
        return i

    endpoints = np.array([[vertex_ids[a], vertex_ids[b]] for a, b in edges])
    for a, b in endpoints:
        ra, rb = root(int(a)), root(int(b))
        if ra != rb:
            component[rb] = ra
    curve = np.array([root(i) for i in range(len(vertices))])
    point_geometries = points(vertices)
    segments = linestrings(edges)
    tolerance = delta * _SEPARATION_RELATIVE_TOLERANCE
    for kind, tree, targets in (
        (0, STRtree(point_geometries), point_geometries),
        (1, STRtree(segments), segments),
    ):
        for start in range(0, len(vertices), 256):
            batch = point_geometries[start : start + 256]
            source, target = tree.query(batch, predicate="dwithin", distance=delta)
            source = source + start
            keep = (
                target > source
                if kind == 0
                else np.all(endpoints[target] != source[:, None], axis=1)
            )
            source, target = source[keep], target[keep]
            if not len(source):
                continue
            short = (
                distance(point_geometries[source], targets[target]) < delta - tolerance
            )
            source, target = source[short], target[short]
            other = (
                curve[source] != curve[target if kind == 0 else endpoints[target][:, 0]]
            )
            name = "vertex_vertex" if kind == 0 else "vertex_edge"
            counts[f"{name}_other_curve"] += int(np.count_nonzero(other))
            counts[f"{name}_same_curve"] += int(len(source) - np.count_nonzero(other))
    return counts


def measure_candidate(occupied, operations, delta, *, with_classes=True):
    """Build one candidate output and record what the checker says about it."""
    name = candidate_name(operations)
    started = time.perf_counter()
    try:
        output = reference_output(occupied, operations)
        state = admissibility(output, delta)
    except GEOSException as error:
        return None, {
            "candidate": name,
            "outcome": "numerical_failure",
            "error": str(error),
            "seconds": time.perf_counter() - started,
        }
    record = {
        "candidate": name,
        "outcome": "admissible" if state["resolved"] else "unresolved",
        "topology_ok": state["topology_ok"],
        "canonical_vertex_count": state.get("canonical_vertex_count"),
        "subscale_pairs": state.get("subscale_pairs"),
        "interior_overlap_pairs": state.get("interior_overlap_pairs"),
        "nonmanifold_boundary_vertices": state.get("nonmanifold_boundary_vertices"),
        "separation_capped_at_delta": state.get("separation_capped_at_delta"),
        "polygon_count": len(output),
    }
    if with_classes and "canonical_vertex_count" in state:
        classes = defect_classes(output, delta)
        total = sum(classes.values())
        if total != state["subscale_pairs"]:
            raise AssertionError(
                f"defect class total {total} disagrees with the checker's "
                f"{state['subscale_pairs']} for candidate {name}"
            )
        record["defect_classes"] = classes
    record["seconds"] = time.perf_counter() - started
    return output, record


def witnessed_epsilon(raw, admissible, delta):
    """Smallest epsilon at which some admissible candidate passes fidelity.

    Admissibility does not depend on epsilon, so the candidate set is fixed
    before the search. Each ladder value prepares one budget and every
    candidate is offered to it, which makes the reported value the minimum over
    candidates rather than the value of whichever candidate was tried first.
    Candidates are offered in order of increasing symmetric difference only to
    reach a pass sooner; the decision is the checker's in either case.
    """
    if not admissible:
        return {"found": False, "reason": "no admissible candidate in the ladder"}
    original = unary_union(polygon_parts(unary_union([p for p in raw])))
    ordered = sorted(
        admissible,
        key=lambda item: original.symmetric_difference(unary_union(item[1])).area,
    )
    budgets = 0

    def passes(epsilon):
        nonlocal budgets
        budgets += 1
        budget = FidelityBudget(raw, epsilon)
        for name, output in ordered:
            if budget.measure(output)["status"] == "pass":
                return name
        return None

    lower = None
    for epsilon in EPSILON_LADDER:
        name = passes(epsilon)
        if name is not None:
            break
        lower = epsilon
    else:
        return {
            "found": False,
            "reason": f"no pass up to epsilon {EPSILON_LADDER[-1]}",
            "budgets": budgets,
        }
    upper, witness = epsilon, name
    if lower is not None:
        for _ in range(EPSILON_REFINEMENTS):
            middle = (lower + upper) / 2
            name = passes(middle)
            if name is None:
                lower = middle
            else:
                upper, witness = middle, name
    return {
        "found": True,
        "epsilon": upper,
        "excluded_below": lower,
        "candidate": witness,
        "budgets": budgets,
    }


def survey_group(raw, delta, *, ladder, detail=False):
    """Measure one interaction group: candidate states and a sufficient epsilon.

    The full candidate table is kept only where it is read, because 57
    candidates on 538 groups is evidence nobody can open. The untreated state,
    which candidates were admissible, the least unresolved candidate and the
    witnessed epsilon are kept for every group.
    """
    started = time.perf_counter()
    occupied, interpretation = _interpret_input(raw)
    records, admissible = [], []
    for index, operations in enumerate(ladder):
        output, record = measure_candidate(
            occupied, operations, delta, with_classes=detail or index == 0
        )
        records.append(record)
        if record["outcome"] == "admissible":
            admissible.append((record["candidate"], output))
    witness = witnessed_epsilon(raw, admissible, delta)
    unresolved = [r for r in records if r["outcome"] == "unresolved"]
    row = {
        "polygon_count": len(raw),
        "interpretation": interpretation,
        "raw_state": records[0],
        "admissible_candidates": [name for name, _ in admissible],
        "least_unresolved": (
            min(
                (
                    {
                        "candidate": r["candidate"],
                        "subscale_pairs": r["subscale_pairs"],
                        "canonical_vertex_count": r["canonical_vertex_count"],
                        "nonmanifold_boundary_vertices": r[
                            "nonmanifold_boundary_vertices"
                        ],
                    }
                    for r in unresolved
                ),
                key=lambda r: (r["nonmanifold_boundary_vertices"], r["subscale_pairs"]),
            )
            if unresolved
            else None
        ),
        "witness": witness,
        "seconds": time.perf_counter() - started,
    }
    if detail:
        row["candidates"] = records
    return row


def baseline_group(raw, delta, epsilon, tolerances):
    """Milestone B control: simplify then canonicalize, judged at the contract.

    One global tolerance, the same for every group, and the contract's own
    delta and epsilon. Nothing is retried, tuned or repaired per case. Tolerance
    zero records the untreated group so the residual has a fixed reference.
    """
    occupied, interpretation = _interpret_input(raw)
    records = []
    for tolerance in (0.0, *tolerances):
        started = time.perf_counter()
        operations = () if tolerance == 0 else (("simplify", tolerance),)
        try:
            output = reference_output(occupied, operations)
            report = check_cleaning_contract(raw, output, delta=delta, epsilon=epsilon)
        except GEOSException as error:
            records.append(
                {
                    "tolerance": tolerance,
                    "status": "numerical_failure",
                    "error": str(error),
                    "seconds": time.perf_counter() - started,
                }
            )
            continue
        state = report["admissibility"]
        record = {
            "tolerance": tolerance,
            "status": report["status"],
            "resolved": state["resolved"],
            "fidelity_status": report["fidelity"]["status"],
            "canonical_vertex_count": state.get("canonical_vertex_count"),
            "subscale_pairs": state.get("subscale_pairs"),
            "nonmanifold_boundary_vertices": state.get("nonmanifold_boundary_vertices"),
            "polygon_count": len(output),
        }
        # Only a fidelity result that is not a clean pass needs its areas on
        # record; writing two zeros for every conforming group buries the ones
        # that matter.
        if report["fidelity"]["status"] != "pass":
            record["lost_protected_area"] = report["fidelity"]["lost_protected_area"]
            record["added_outside_budget_area"] = report["fidelity"][
                "added_outside_budget_area"
            ]
        if state.get("interior_overlap_pairs"):
            record["interior_overlap_pairs"] = state["interior_overlap_pairs"]
        if "canonical_vertex_count" in state:
            classes = defect_classes(output, delta)
            if sum(classes.values()) != state["subscale_pairs"]:
                raise AssertionError(
                    "defect class total disagrees with the checker at tolerance "
                    f"{tolerance}"
                )
            record["defect_classes"] = classes
        record["seconds"] = time.perf_counter() - started
        records.append(record)
    return {
        "polygon_count": len(raw),
        "interpretation": interpretation,
        "tolerances": records,
    }


def survey(groups, delta, *, ladder, label, detailed=(), progress=True):
    rows = []
    for index, raw in enumerate(groups):
        row = survey_group(raw, delta, ladder=ladder, detail=index in detailed)
        row["group"] = index
        rows.append(row)
        if progress and (row["seconds"] > 2 or index % 25 == 0):
            print(
                f"{label} group {index}: "
                f"{'witnessed at ' + format(row['witness']['epsilon'], '.4g') if row['witness']['found'] else 'none found'}"
                f" ({row['seconds']:.2f}s)",
                flush=True,
            )
    return rows


def recorded_point_contact_witness():
    """The analytic corner-cut witness already recorded for `point_contact`.

    Measured here only to separate two different meanings of a failed sweep. It
    is an inherited result, not a ladder candidate and not a tuned one: the
    sampled leg length 0.375 is copied unchanged from the automatic partition
    probe. The ladder never sees it.
    """
    from shapely.geometry import Polygon, box

    raw = [box(0, 0, 2, 2), box(2, 2, 4, 4)]
    t = 0.375
    cuts = unary_union(
        [
            Polygon([(2, 2), (2 - t, 2), (2, 2 - t)]),
            Polygon([(2, 2), (2 + t, 2), (2, 2 + t)]),
        ]
    )
    return raw, polygon_parts(unary_union(raw).difference(cuts))


def examples_survey(ladder):
    from sandbox.cleaning_pipeline_probe import topology_examples

    rows = {
        name: survey_group(raw, 0.5, ladder=ladder, detail=True)
        for name, (raw, _) in topology_examples().items()
    }
    raw, witness = recorded_point_contact_witness()
    state = admissibility(witness, 0.5)
    rows["point_contact"]["recorded_witness"] = {
        "source": "sampled corner cuts, leg 0.375, from sandbox/cleaning_pipeline_probe.py",
        "outcome": "admissible" if state["resolved"] else "unresolved",
        "subscale_pairs": state.get("subscale_pairs"),
        "separation_capped_at_delta": state.get("separation_capped_at_delta"),
        "witness": witnessed_epsilon(raw, [("recorded_witness", witness)], 0.5),
    }
    return rows


def scaling_survey(run, destination, *, budget=None, skip=()):
    """Milestone C: per-group work against canonical vertex count.

    The manifests record per-city totals only, so the per-group field the fit
    needs does not exist yet and the recorded construction is re-run unchanged
    to add it. `sandbox/cleaning_patch_probe.py` is imported, not edited: its
    hash is part of the earlier manifests' provenance. The construction's own
    progress history already carries the canonical vertex count before the
    first edit, so the vertex count and the work are measured on the same run.

    Results append to a JSON-lines file and finished groups are never redone,
    because the measurement harness cannot hold a process for the whole corpus.
    A budget stops the pass between groups rather than inside one, so no
    partial group is ever recorded. `skip` names groups this harness cannot
    time at all; they are reported as omitted rather than silently missing.
    """
    from sandbox.cleaning_patch_probe import construct
    from sandbox.cleaning_pipeline_probe import saved_city_groups

    started = time.perf_counter()
    done = set()
    if destination.exists():
        for line in destination.read_text().splitlines():
            row = json.loads(line)
            done.add((row["city"], row["group"]))
    for record, delta, epsilon, groups in saved_city_groups(run):
        city = record["case"]["city"]
        for index, raw in enumerate(groups):
            if (city, index) in done or f"{city}:{index}" in skip:
                continue
            if budget is not None and time.perf_counter() - started > budget:
                print(f"budget reached before {city} group {index}", flush=True)
                return False
            group_started = time.perf_counter()
            output, report = construct(
                raw, delta=delta, epsilon=epsilon, merge_buildings=True
            )
            progress = report["progress"]
            row = {
                "city": city,
                "group": index,
                "polygon_count": len(raw),
                "initial_vertices": progress[0][2],
                "final_vertices": progress[-1][2],
                "initial_subscale_pairs": progress[0][1],
                "final_subscale_pairs": progress[-1][1],
                "edits": report["edits"],
                "checks": report["checks"],
                "outcome": report["outcome"],
                "reason": report["reason"],
                "seconds": round(time.perf_counter() - group_started, 6),
            }
            with destination.open("a") as handle:
                handle.write(json.dumps(row) + "\n")
            print(
                f"{city} {index}: {row['outcome']} v={row['initial_vertices']} "
                f"edits={row['edits']} {row['seconds']:.2f}s",
                flush=True,
            )
    return True


def power_fit(pairs, label):
    """Least squares on log-log, with the fit quality that decides how to read it."""
    pairs = [(x, y) for x, y in pairs if x > 0 and y > 0]
    if len(pairs) < 3:
        return {"label": label, "count": len(pairs), "exponent": None}
    x = np.log(np.array([a for a, _ in pairs], dtype=float))
    y = np.log(np.array([b for _, b in pairs], dtype=float))
    design = np.vstack([x, np.ones_like(x)]).T
    (slope, intercept), *_ = np.linalg.lstsq(design, y, rcond=None)
    predicted = design @ [slope, intercept]
    residual = float(((y - predicted) ** 2).sum())
    total = float(((y - y.mean()) ** 2).sum())
    return {
        "label": label,
        "count": len(pairs),
        "exponent": float(slope),
        "coefficient": float(np.exp(intercept)),
        "r_squared": 1 - residual / total if total else None,
        "log_rmse": float(np.sqrt(residual / len(pairs))),
    }


def scaling_fits(rows):
    """The fits milestone C asks for, plus the ones that say why they differ.

    Work-limited groups are censored: the 2000-proposal cap truncates both
    their time and their accepted edits, and every one of the largest groups is
    censored. Both the all-groups and the uncensored fits therefore understate
    growth, in opposite ways, and the per-proposal cost is reported because it
    is the one quantity the cap does not truncate.
    """
    free = [r for r in rows if r["reason"] != "work_limit"]
    return {
        "censoring": {
            "work_limited_groups": len(rows) - len(free),
            "largest_group_is_censored": all(
                r["reason"] == "work_limit"
                for r in sorted(rows, key=lambda r: -r["initial_vertices"])[:10]
            ),
            "reading": (
                "Every fit below is a lower bound on growth: the cap removes "
                "exactly the work the largest groups would otherwise do."
            ),
        },
        "fits": [
            power_fit(
                [(r["initial_vertices"], r["seconds"]) for r in rows],
                "runtime against canonical vertex count, all groups",
            ),
            power_fit(
                [(r["initial_vertices"], r["edits"]) for r in rows],
                "accepted edits against canonical vertex count, all groups",
            ),
            power_fit(
                [(r["initial_vertices"], r["seconds"]) for r in free],
                "runtime against canonical vertex count, uncensored",
            ),
            power_fit(
                [(r["initial_vertices"], r["edits"]) for r in free],
                "accepted edits against canonical vertex count, uncensored",
            ),
            power_fit(
                [(r["initial_vertices"], r["checks"]) for r in rows],
                "proposals checked against canonical vertex count",
            ),
            power_fit(
                [
                    (r["initial_vertices"], r["seconds"] / r["checks"])
                    for r in rows
                    if r["checks"]
                ],
                "seconds per proposal against canonical vertex count",
            ),
            power_fit(
                [(r["initial_subscale_pairs"], r["edits"]) for r in rows],
                "accepted edits against initial short pairs",
            ),
            power_fit(
                [(r["edits"], r["seconds"]) for r in rows],
                "runtime against accepted edits",
            ),
        ],
        "concentration": {
            "measured_seconds": sum(r["seconds"] for r in rows),
            "share_in_slowest": {
                str(k): sum(
                    r["seconds"] for r in sorted(rows, key=lambda r: -r["seconds"])[:k]
                )
                / sum(r["seconds"] for r in rows)
                for k in (1, 5, 10, 20)
            },
        },
    }


# Cumulative time under one named frame each. The first two are disjoint: the
# descent calls them separately per proposal, and neither appears inside the
# other. The last two are components that fall inside those, not terms beside
# them, because the solver runs inside a proposal and the tree is queried from
# both. The one-off final `check_cleaning_contract` is left out; it is one call
# against thousands.
PROFILE_TERMS = {
    "conformance check": (
        ("cleaning/contract.py", "admissibility"),
        ("cleaning/contract.py", "accepts_union"),
    ),
    "proposal construction": (
        ("cleaning_patch_probe.py", "patch_batches"),
        ("cleaning_patch_probe.py", "local_edits"),
        ("cleaning_patch_probe.py", "coupled_move"),
        ("cleaning_patch_probe.py", "move_boundary"),
    ),
    "of which linear program": (("_minimize.py", "minimize"),),
    "of which neighbourhood query": (("strtree.py", "query"),),
}


def profile_group(run, city, index):
    """Attribute one group's construction time to the terms milestone C names.

    Each term is the cumulative time under named frames that the descent enters
    directly, so nothing inside a term is counted twice. The two disjoint terms
    are the conformance check and proposal construction; the linear program and
    the neighbourhood query are reported as components inside them, which is
    what they are.
    """
    import cProfile
    import pstats

    from sandbox.cleaning_patch_probe import construct
    from sandbox.cleaning_pipeline_probe import saved_city_groups

    for record, delta, epsilon, groups in saved_city_groups(run):
        if record["case"]["city"] != city:
            continue
        raw = groups[index]
        profiler = cProfile.Profile()
        started = time.perf_counter()
        profiler.enable()
        _, report = construct(raw, delta=delta, epsilon=epsilon, merge_buildings=True)
        profiler.disable()
        elapsed = time.perf_counter() - started
        stats = pstats.Stats(profiler).stats
        terms, matched = {}, {}
        for term, patterns in PROFILE_TERMS.items():
            for file, function in patterns:
                for (path, _, name), entry in stats.items():
                    if file in str(path).replace("\\", "/") and function == name:
                        terms[term] = terms.get(term, 0.0) + entry[3]
                        matched.setdefault(term, []).append(f"{Path(path).name}:{name}")
        return {
            "city": city,
            "group": index,
            "profiled_seconds": elapsed,
            "outcome": report["outcome"],
            "reason": report["reason"],
            "checks": report["checks"],
            "edits": report["edits"],
            "cumulative_seconds_by_term": terms,
            "entry_points": {k: sorted(v) for k, v in matched.items()},
            "note": (
                "cProfile inflates absolute time; read the shares. The two "
                "terms without an 'of which' prefix are disjoint and together "
                "cover the descent; the other two lie inside them."
            ),
        }
    raise ValueError(f"city {city} not found in {run}")


def run_profile(args):
    """Milestone C: where one group's construction time goes."""
    profiles = []
    for item in args.profile.split(","):
        city, index = item.split(":")
        profiles.append(profile_group(args.run, city, int(index)))
        print(
            city,
            index,
            {
                k: round(v, 2)
                for k, v in profiles[-1]["cumulative_seconds_by_term"].items()
            },
            flush=True,
        )
    write_json(
        args.output / "profile.json",
        {"run": str(args.run), **provenance(), "profiles": profiles},
    )


def run_scaling(args):
    """Milestone C: per-group work, resumable, assembled once it is complete."""
    skip = {s for s in args.skip.split(",") if s}
    lines = args.output / "scaling.jsonl"
    complete = scaling_survey(args.run, lines, budget=args.budget, skip=skip)
    rows = [json.loads(line) for line in lines.read_text().splitlines()]
    print(f"{len(rows)} groups recorded; complete={complete}", flush=True)
    if not complete:
        return
    write_json(
        args.output / "scaling.json",
        {
            "scope": (
                "Milestone C: the recorded neighbouring-clearance construction "
                "re-run unchanged to add per-group canonical vertex counts, "
                "accepted edits and time. Conforming counts are expected to "
                "reproduce the manifest exactly, and the per-city totals double "
                "as the calibration between this machine and the one the "
                "manifests were measured on."
            ),
            "run": str(args.run),
            "omitted_groups": sorted(skip),
            **provenance(args.run),
            **scaling_fits(rows),
            "profiles": [
                profile
                for path in sorted(args.output.glob("profile*.json"))
                for profile in json.loads(path.read_text())["profiles"]
            ],
            "groups": rows,
        },
    )


def run_sweep(args, ladder):
    """Milestones A and B: the feasibility sweep and the declared control."""
    from sandbox.cleaning_pipeline_probe import saved_city_groups

    wanted = {c for c in args.cities.split(",") if c}
    examples = examples_survey(ladder)
    write_json(
        args.output / "examples.json",
        {"ladder": [candidate_name(c) for c in ladder], "examples": examples},
    )
    for name, row in examples.items():
        print(
            name,
            "witnessed at" if row["witness"]["found"] else "none found",
            row["witness"].get("epsilon", ""),
            row["witness"].get("candidate", ""),
            flush=True,
        )
    if not args.run:
        return
    cities, baseline = [], []
    for record, delta, epsilon, groups in saved_city_groups(args.run):
        city = record["case"]["city"]
        if wanted and city not in wanted:
            continue
        rows = survey(
            groups,
            delta,
            ladder=ladder,
            label=city,
            detailed=unresolved_groups(city),
        )
        for row in list(rows):
            # A sweep that found nothing is the one case where the whole
            # candidate table has to be readable afterwards.
            if row["witness"]["found"] or "candidates" in row:
                continue
            index = row["group"]
            rows[index] = survey_group(groups[index], delta, ladder=ladder, detail=True)
            rows[index]["group"] = index
        cities.append(
            {
                "case": record["case"],
                "delta": delta,
                "contract_epsilon": epsilon,
                "groups": rows,
                "seconds": sum(r["seconds"] for r in rows),
            }
        )
        witnessed = sum(r["witness"]["found"] for r in rows)
        within = sum(
            r["witness"]["found"] and r["witness"]["epsilon"] <= epsilon for r in rows
        )
        print(
            f"{city}: {witnessed}/{len(rows)} witnessed, "
            f"{within} at epsilon <= {epsilon} ({cities[-1]['seconds']:.1f}s)",
            flush=True,
        )
        control = []
        for index, raw in enumerate(groups):
            row = baseline_group(raw, delta, epsilon, BASELINE_TOLERANCES)
            row["group"] = index
            control.append(row)
        baseline.append(
            {
                "case": record["case"],
                "delta": delta,
                "epsilon": epsilon,
                "groups": control,
                "seconds": sum(t["seconds"] for r in control for t in r["tolerances"]),
            }
        )
        passes = {
            t: sum(
                any(
                    x["tolerance"] == t and x["status"] == "pass"
                    for x in r["tolerances"]
                )
                for r in control
            )
            for t in BASELINE_TOLERANCES
        }
        print(f"{city} baseline pass counts by tolerance: {passes}", flush=True)
        write_json(
            args.output / "cities.json",
            {
                "run": str(args.run),
                "scope": (
                    "Milestone A feasibility sweep. A disposable global "
                    "construction witnesses a sufficient epsilon; it is not "
                    "a candidate cleaner and its epsilon is an upper bound."
                ),
                "delta": cities[0]["delta"],
                "contract_epsilon": cities[0]["contract_epsilon"],
                "epsilon_ladder": list(EPSILON_LADDER),
                "epsilon_refinements": EPSILON_REFINEMENTS,
                "ladder": [candidate_name(c) for c in ladder],
                **provenance(args.run),
                "analytic_examples": examples,
                "cities": cities,
            },
        )
        write_json(
            args.output / "baseline.json",
            {
                "run": str(args.run),
                "control": (
                    "Milestone B declared control: GEOS topology-preserving "
                    "simplification at one global tolerance, then "
                    "canonicalization, judged by the independent checker at "
                    "the contract's delta and epsilon. Tolerance 0 is the "
                    "untreated group."
                ),
                "tolerances": [0.0, *BASELINE_TOLERANCES],
                **provenance(args.run),
                "cities": baseline,
            },
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--run", type=Path)
    parser.add_argument("--cities", type=str, default="")
    parser.add_argument("--scaling", action="store_true")
    parser.add_argument("--budget", type=float, default=None)
    parser.add_argument("--skip", type=str, default="")
    parser.add_argument("--profile", type=str, default="")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.profile:
        run_profile(args)
    elif args.scaling:
        run_scaling(args)
    else:
        run_sweep(args, candidate_ladder())


if __name__ == "__main__":
    main()
