from __future__ import annotations

import importlib
import importlib.util
import os
from typing import Literal

MesherName = Literal["auto", "dtcc_mesher", "triangle", "spade"]

_DEFAULT_2D_MESHER_ENV = "DTCC_2D_MESHER"
_default_2d_mesher_override: MesherName | None = None


def _normalize_mesher_name(mesher: str) -> MesherName:
    normalized = mesher.strip().lower()
    if normalized not in {"auto", "dtcc_mesher", "triangle", "spade"}:
        raise ValueError(
            f"Unsupported 2D mesher '{mesher}'. "
            f"Expected one of: auto, dtcc_mesher, triangle, spade."
        )
    return normalized  # type: ignore[return-value]


def _builder_backend_available(name: str) -> bool:
    spec = importlib.util.find_spec("dtcc_core.builder._dtcc_builder")
    if spec is None:
        return False

    try:
        builder_module = importlib.import_module("dtcc_core.builder._dtcc_builder")
    except Exception:
        return False

    if hasattr(builder_module, "triangulation_backends"):
        try:
            return name in builder_module.triangulation_backends()
        except Exception:
            return False

    if name == "triangle":
        return bool(getattr(builder_module, "HAVE_TRIANGLE", False))
    if name == "spade":
        return bool(getattr(builder_module, "HAVE_SPADE", False))
    return False


def available_2d_meshers() -> list[str]:
    meshers: list[str] = []
    if importlib.util.find_spec("dtcc_mesher") is not None:
        meshers.append("dtcc_mesher")
    if _builder_backend_available("triangle"):
        meshers.append("triangle")
    if _builder_backend_available("spade"):
        meshers.append("spade")
    return meshers


def resolve_2d_mesher(mesher: str | None = None) -> str:
    requested = mesher
    if requested is None:
        requested = _default_2d_mesher_override or os.getenv(
            _DEFAULT_2D_MESHER_ENV, "auto"
        )

    normalized = _normalize_mesher_name(requested)
    if normalized == "auto":
        available = available_2d_meshers()
        if "dtcc_mesher" in available:
            return "dtcc_mesher"
        if "triangle" in available:
            return "triangle"
        if "spade" in available:
            return "spade"
        raise RuntimeError("No supported 2D mesher backend is available.")

    available = available_2d_meshers()
    if normalized not in available:
        raise RuntimeError(
            f"Requested 2D mesher '{normalized}' is not available. "
            f"Available meshers: {', '.join(available) or 'none'}."
        )

    return normalized


def get_default_2d_mesher() -> str:
    return resolve_2d_mesher()


def set_default_2d_mesher(mesher: str | None = None) -> str:
    global _default_2d_mesher_override

    if mesher is None:
        _default_2d_mesher_override = None
        return resolve_2d_mesher()

    normalized = _normalize_mesher_name(mesher)
    if normalized == "auto":
        _default_2d_mesher_override = "auto"
        return resolve_2d_mesher()

    resolved = resolve_2d_mesher(normalized)
    _default_2d_mesher_override = normalized
    return resolved
