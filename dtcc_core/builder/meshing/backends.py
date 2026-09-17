from __future__ import annotations

import importlib
import importlib.util
import os
from typing import Literal, cast

MesherName = Literal["auto", "dtcc_mesher", "triangle"]

_DEFAULT_2D_MESHER_ENV = "DTCC_2D_MESHER"
_default_2d_mesher_override: MesherName | None = None
_MESHER_CHOICES: tuple[str, ...] = ("auto", "dtcc_mesher", "triangle")
_BUILDER_MESHER_CHOICES: tuple[str, ...] = ("triangle",)


def _normalize_mesher_name(mesher: str) -> MesherName:
    normalized = mesher.strip().lower()
    if normalized not in _MESHER_CHOICES:
        raise ValueError(
            f"Unsupported 2D mesher '{mesher}'. "
            f"Expected one of: {', '.join(_MESHER_CHOICES)}."
        )
    return cast(MesherName, normalized)


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

    if name in _BUILDER_MESHER_CHOICES:
        return bool(getattr(builder_module, "HAVE_TRIANGLE", False))
    return False


def available_2d_meshers() -> list[str]:
    """Return the 2D mesher backends usable in this installation.

    ``"dtcc_mesher"`` is listed when the ``dtcc_mesher`` package is
    installed, and ``"triangle"`` when dtcc-core was built with Triangle
    support.

    Returns
    -------
    list[str]
        Available mesher names, in order of preference.
    """
    meshers: list[str] = []
    if importlib.util.find_spec("dtcc_mesher") is not None:
        meshers.append("dtcc_mesher")
    for name in _BUILDER_MESHER_CHOICES:
        if _builder_backend_available(name):
            meshers.append(name)
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
        raise RuntimeError("No supported 2D mesher backend is available.")

    available = available_2d_meshers()
    if normalized not in available:
        raise RuntimeError(
            f"Requested 2D mesher '{normalized}' is not available. "
            f"Available meshers: {', '.join(available) or 'none'}."
        )

    return normalized


def get_default_2d_mesher() -> str:
    """Return the 2D mesher used when a meshing function gets ``mesher=None``.

    The default is the mesher set with ``set_default_2d_mesher``, otherwise
    the ``DTCC_2D_MESHER`` environment variable, otherwise ``"auto"``.
    ``"auto"`` picks ``"dtcc_mesher"`` when available, then ``"triangle"``.

    Returns
    -------
    str
        The concrete mesher name, never ``"auto"``.

    Raises
    ------
    ValueError
        If the configured name is not a supported mesher.
    RuntimeError
        If the configured mesher is not available, or no mesher is available
        for ``"auto"``.
    """
    return resolve_2d_mesher()


def set_default_2d_mesher(mesher: str | None = None) -> str:
    """Set the 2D mesher used when a meshing function gets ``mesher=None``.

    The setting lasts for the current Python process and takes precedence
    over the ``DTCC_2D_MESHER`` environment variable.

    Parameters
    ----------
    mesher : {"auto", "dtcc_mesher", "triangle"}, optional
        Mesher to use by default. ``None`` clears the setting, so the
        environment variable or ``"auto"`` applies again.

    Returns
    -------
    str
        The concrete mesher the new default resolves to.

    Raises
    ------
    ValueError
        If ``mesher`` is not a supported name.
    RuntimeError
        If ``mesher`` is not available in this installation. The previous
        default is then kept.
    """
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
