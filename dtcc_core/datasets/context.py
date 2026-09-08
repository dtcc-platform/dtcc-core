"""Helpers for attaching Dataset v2 context to native DTCC objects."""

from __future__ import annotations

from typing import Any

from dtcc_core.model.model import Model as DTCCModel

from .schema import DatasetContext


def attach_dataset_context(obj: Any, context: DatasetContext) -> Any:
    """Attach dataset context to native DTCC model objects when possible.

    Bare Python containers are returned unchanged until their datasets migrate
    to semantic model objects or transitional Dataset v2 containers.
    """
    if isinstance(obj, DTCCModel):
        obj.dataset_context = context
        return obj

    # TODO(Dataset v2): migrate any remaining bare list results to
    # dataset-aware model objects.
    if isinstance(obj, list):
        return obj

    # TODO(Dataset v2): migrate any remaining bare dict results to
    # dataset-aware model objects.
    if isinstance(obj, dict):
        return obj

    return obj
