"""Proof examples for finite occupancy selection, not a production cleaner."""

import pytest
from shapely.geometry import Point, Polygon, box
from shapely.ops import unary_union

from dtcc_core.builder.cleaning.contract import check_cleaning_contract
from sandbox.cleaning_pipeline_probe import (
    automatic_carrier,
    automatic_probe,
    example_carrier,
    occupancy_search,
    parts,
    topology_examples,
)


def test_one_occupancy_rule_handles_contacts_passages_and_holes():
    for name, (raw, cells) in topology_examples().items():
        original = [p.wkb for p in raw]
        output, report = occupancy_search(raw, cells)
        assert report["outcome"] == "conforming", (name, report)
        assert report["contract"]["status"] == "pass"
        assert report["checked"] <= report["assignments"] == 2 ** report["free"]
        assert [p.wkb for p in raw] == original
        if name == "courtyard_passage":
            assert len(output) == 1 and len(output[0].interiors) == 1
            assert not output[0].covers(Point(6, 5))
        elif name == "tiny_hole":
            assert len(output) == 1 and not output[0].interiors
        elif name == "point_contact":
            assert len(output) == 2 and output[0].distance(output[1]) >= 0.5
        else:
            assert len(output) == 1


def test_contact_requires_a_joint_choice_and_insufficient_carrier_is_unresolved():
    raw, cells = topology_examples()["point_contact"]
    corner = Polygon([(2, 2), (1.6, 2), (2, 1.6)])
    one_change = parts(unary_union(raw).difference(corner))
    partial = check_cleaning_contract(raw, one_change, delta=0.5, epsilon=0.25)
    assert partial["fidelity"]["status"] == "pass"
    assert not partial["admissibility"]["resolved"]
    output, limited = occupancy_search(raw, example_carrier(raw))
    assert output is None and limited["reason"] == "carrier_search_exhausted"
    output, refined = occupancy_search(raw, cells)
    assert refined["outcome"] == "conforming"
    # Permuting inputs or cells does not select a different geometric solution.
    reversed_output, _ = occupancy_search(raw[::-1], cells[::-1])
    assert unary_union(output).equals(unary_union(reversed_output))


def test_fixed_band_rejects_bad_carriers_and_does_not_relax_to_find_a_solution():
    raw, cells = topology_examples()["point_contact"]
    output, strict = occupancy_search(raw, cells, epsilon=0)
    assert output is None and strict["reason"] == "carrier_search_exhausted"
    _, mixed = occupancy_search(raw, [unary_union(raw).convex_hull])
    assert mixed["reason"] == "carrier_cell_conflict"
    _, missing = occupancy_search(raw, [raw[0]])
    assert missing["reason"] == "carrier_omits_protected_space"


def test_resolved_input_is_unchanged_when_represented_in_carrier():
    raw = [box(0, 0, 3, 3)]
    output, report = occupancy_search(raw, example_carrier(raw))
    assert report["outcome"] == "conforming"
    assert unary_union(output).equals(raw[0])
    assert report["changed_area"] == 0


def test_work_limit_and_malformed_carrier_are_explicit():
    raw, cells = topology_examples()["point_contact"]
    output, report = occupancy_search(raw, cells, max_free_cells=1)
    assert output is None and report["reason"] == "work_limit"
    assert report["checked"] == 0
    with pytest.raises(ValueError, match="interior-disjoint"):
        occupancy_search(raw, [box(0, 0, 3, 3), box(1, 1, 4, 4)])
    with pytest.raises(ValueError, match="valid"):
        occupancy_search(raw, [Polygon([(0, 0), (1, 1), (0, 1), (1, 0)])])
    with pytest.raises(ValueError, match="epsilon"):
        occupancy_search(raw, cells, epsilon=-1)
    with pytest.raises(ValueError, match="max_free_cells"):
        occupancy_search(raw, cells, max_free_cells=13)


def test_automatic_partitions_represent_input_and_conforming_alternatives():
    # Witnesses are chosen after construction, never supplied to the generator.
    for name, report in automatic_probe().items():
        assert report["reason"] == "constructed", (name, report)
        assert report["original"]["reconstruction_difference_area"] == 0
        assert report["witness"]["reconstruction_difference_area"] == 0
        assert report["witness"]["contract"]["status"] == "pass"


def test_automatic_partition_is_independent_of_input_order():
    raw = topology_examples()["point_contact"][0]
    original = [p.wkb for p in raw]
    forward, _ = automatic_carrier(raw)
    reverse, _ = automatic_carrier(raw[::-1])
    assert [c.normalize().wkb for c in forward] == [c.normalize().wkb for c in reverse]
    assert [p.wkb for p in raw] == original


def test_automatic_partition_preserves_resolved_input_and_returns_no_partial_carrier():
    raw = [box(0, 0, 3, 3)]
    cells, report = automatic_carrier(raw)
    assert report["initially_resolved"] and report["chords"] == 0
    assert unary_union(cells).equals(raw[0])
    touching = topology_examples()["point_contact"][0]
    cells, report = automatic_carrier(touching, max_chords=1)
    assert cells is None and report["reason"] == "chord_limit"
    with pytest.raises(ValueError, match="epsilon"):
        automatic_carrier(raw, epsilon=-1)
