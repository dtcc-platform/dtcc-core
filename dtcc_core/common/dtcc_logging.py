# Copyright(C) 2023 Anders Logg
# Licensed under the MIT License

import logging as _logging
from typing import Dict, Tuple, Callable

from rich.console import Console

from .handler import LoggingHandler

# Import the shared console from progress module for coordination
# This ensures log messages appear above the progress bar
try:
    from .progress import get_console
    _console = get_console()
except ImportError:
    # Fallback if progress module not available
    _console = Console(stderr=True,
                       soft_wrap=False)

# Global logger dictionary (name -> function tuple)
loggers: Dict[str, Tuple[Callable, Callable, Callable, Callable, Callable]] = {}

# Global logger objects (name -> logging.Logger)
_logger_objects: Dict[str, _logging.Logger] = {}

# Global logger object
_logger: _logging.Logger = None


def _coerce_level(level: str | int) -> int:
    """Convert a string/int log level to logging level int."""
    if isinstance(level, int):
        return level
    if isinstance(level, str):
        value = _logging.getLevelName(level.upper())
        if isinstance(value, int):
            return value
    raise ValueError(
        f"Invalid log level {level!r}. Valid levels are DEBUG, INFO, WARNING, ERROR, CRITICAL."
    )


_global_log_level: int = _logging.INFO


def _ensure_logger(name: str) -> _logging.Logger:
    """Return a configured logger with DTCC formatting for this source name."""
    logger = _logger_objects.get(name)
    if logger is not None:
        return logger

    logger = _logging.getLogger(name)
    logger.setLevel(_global_log_level)
    logger.propagate = False

    # Remove stale DTCC handlers for this source before attaching one.
    logger.handlers = [
        handler
        for handler in logger.handlers
        if not (isinstance(handler, LoggingHandler) and handler.source_name == name)
    ]
    logger.addHandler(
        LoggingHandler(
            source_name=name,
            console=_console,
        )
    )

    _logger_objects[name] = logger
    return logger


def _callbacks_for(logger: _logging.Logger):
    """Build helper callbacks with legacy DTCC API."""
    def error(message):
        logger.error(message)
        raise RuntimeError(message)

    def critical(message):
        logger.critical(message)
        raise RuntimeError(message)

    return (logger.debug, logger.info, logger.warning, error, critical)


def _init_logging(name):
    "Internal function for initializing logging."

    global _logger

    _logger = _ensure_logger(name)
    return _callbacks_for(_logger)


debug, info, warning, error, critical = _init_logging("dtcc-core")


def init_logging(name="dtcc-core"):
    "Initialize logging for given package."
    callbacks = _init_logging(name)
    loggers[name] = callbacks
    return callbacks


def get_logger(name="dtcc-core"):
    "Get logger callback tuple for given package."
    if name not in loggers:
        loggers[name] = init_logging(name)
    return loggers[name]


def get_python_logger(name: str = "dtcc-core") -> _logging.Logger:
    """Get configured Python logger object for a given package source name."""
    return _ensure_logger(name)


def set_log_level(level):
    """Set log level for all configured DTCC loggers."""
    global _logger
    global _global_log_level

    _global_log_level = _coerce_level(level)

    if _logger is None:
        _logger = _ensure_logger("dtcc-core")
    _logger.setLevel(_global_log_level)

    for logger in _logger_objects.values():
        logger.setLevel(_global_log_level)
