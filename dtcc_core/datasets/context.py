"""Helpers for attaching Dataset v2 context to native DTCC objects."""

from __future__ import annotations

from typing import Any

from dtcc_core.model.model import Model as DTCCModel

from .schema import DatasetContext


def attach_dataset_context(obj: Any, context: DatasetContext) -> Any:
    """Attach dataset context to native DTCC model objects when possible.

    Bare Python containers are returned unchanged until their datasets migrate
    to typed DatasetCollection or DatasetValue results.
    """
    if isinstance(obj, DTCCModel):
        obj.dataset_context = context
        return obj

    # TODO(Phase 1C): migrate remaining bare list results to typed dataset-aware
    # collection models.
    if isinstance(obj, list):
        return obj

    # TODO(Phase 1C): migrate remaining bare dict results to typed dataset-aware
    # value/layer models.
    if isinstance(obj, dict):
        return obj

    return obj
