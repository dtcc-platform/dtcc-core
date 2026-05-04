from __future__ import annotations

import importlib

from dtcc_core.common.dtcc_logging import (
    init_logging,
    set_log_level,
    get_logger,
    get_python_logger,
)
from dtcc_core import common
from dtcc_core import model
from dtcc_core import builder

register_model_method = builder.register_model_method
Bounds = model.Bounds

_LAZY_SUBMODULES = {
    "io": "dtcc_core.io",
    "datasets": "dtcc_core.datasets",
    "plotting": "dtcc_core.plotting",
    "reproject": "dtcc_core.reproject",
}


def __getattr__(name: str):
    if name in _LAZY_SUBMODULES:
        module = importlib.import_module(_LAZY_SUBMODULES[name])
        globals()[name] = module
        return module
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


__all__ = [
    "init_logging",
    "set_log_level",
    "get_logger",
    "get_python_logger",
    "common",
    "model",
    "builder",
    "Bounds",
    "register_model_method",
    "io",
    "datasets",
    "plotting",
    "reproject",
]
