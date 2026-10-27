"""Synthetic checks for the sandbox catalogue model and exact patch transfer."""

import sys
from pathlib import Path

import numpy as np
import pytest
from shapely import affinity
from shapely.geometry import Polygon, box
from shapely.ops import unary_union

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "sandbox"))
from roof_catalogue import RoofModel, fit_catalogue, planes


@pytest.mark.parametrize(
    "kind,parameters",
    [
        ("flat", [8.0]),
        ("shed", [8.0, 0.15, -0.25]),
        ("gable", [0.0, 8.0, 0.55]),
        ("hip", [0.0, 8.0, 0.55, 0.6]),
    ],
)
def test_recovers_noisy_catalogue_roofs(kind, parameters):
    rng = np.random.default_rng(2026)
    footprint = box(-12, -6, 12, 6)
    xy = rng.uniform([-12, -6], [12, 6], size=(900, 2))
    pp = planes(kind, np.array(parameters), np.array(footprint.exterior.coords))
    z = np.min(xy @ pp[:, :2].T + pp[:, 2], axis=1) + rng.normal(0, 0.035, len(xy))
    z[:30] += 3
    model = fit_catalogue(footprint, np.c_[xy, z], 0)
    truth = np.min(xy @ pp[:, :2].T + pp[:, 2], axis=1)
    assert model.kind == kind
    assert np.median(np.abs(model.predict(xy) - truth)) < 0.08
    assert model.diagnostics["supported"]


def test_partition_preserves_courtyard_and_covers_changed_footprint():
    original = box(-10, -5, 10, 5)
    pp = planes("hip", [0, 8, 0.5, 0.7], np.array(original.exterior.coords))
    model = RoofModel("hip", np.array([0, 8, 0.5, 0.7]), np.zeros(2), pp, {})
    target = Polygon(
        [(-11, -6), (11, -6), (11, 6), (-11, 6)],
        holes=[[(-1, -1), (1, -1), (1, 1), (-1, 1)]],
    )
    patches = model.patches(target)
    union = unary_union([p["polygon"] for p in patches])
    assert union.symmetric_difference(target).area < 1e-8
    assert sum(p["polygon"].area for p in patches) - union.area < 1e-8
    for p in patches:
        xy = np.array(p["polygon"].exterior.coords)
        assert np.allclose(xy @ p["coeff"][:2] + p["coeff"][2], model.predict(xy))


def test_translation_and_rotation_invariance():
    rng = np.random.default_rng(21)
    foot = box(-10, -5, 10, 5)
    xy = rng.uniform([-10, -5], [10, 5], size=(500, 2))
    z = 8 + 0.6 * (5 - np.abs(xy[:, 1]))
    angle = 0.63
    rotation = np.array(
        [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    )
    offset = np.array([386000, 6174000])
    world = xy @ rotation.T + offset
    footprint = affinity.translate(
        affinity.rotate(foot, np.rad2deg(angle), origin=(0, 0)), *offset
    )
    model = fit_catalogue(footprint, np.c_[world, z], 0)
    assert model.kind == "gable"
    assert np.max(np.abs(model.predict(world) - z)) < 0.03


def test_sparse_cloud_is_not_fitted():
    assert fit_catalogue(box(0, 0, 5, 5), np.zeros((11, 3)), 0) is None


def test_compound_height_step_is_flagged_for_review():
    rng = np.random.default_rng(93)
    xy = rng.uniform([-10, -5], [10, 5], size=(600, 2))
    z = np.where(xy[:, 0] < 0, 6.0, 12.0) + rng.normal(0, 0.03, len(xy))
    model = fit_catalogue(box(-10, -5, 10, 5), np.c_[xy, z], 0)
    assert not model.diagnostics["supported"]
    assert model.diagnostics["review_reasons"]
