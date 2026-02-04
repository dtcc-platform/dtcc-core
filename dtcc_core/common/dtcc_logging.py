# Copyright(C) 2023 Anders Logg
# Licensed under the MIT License

import logging as _logging
import sys

from rich.logging import RichHandler
from rich.console import Console

# Import the shared console from progress module for coordination
# This ensures log messages appear above the progress bar
try:
    from .progress import get_console
    _console = get_console()
except ImportError:
    # Fallback if progress module not available
    _console = Console(stderr=True)

# Global logger dictionary
loggers = {}

# Global logger object
_logger = None


def _init_logging(name):
    "Internal function for initializing logging"

    global _logger

    # Initialize logger
    _logger = _logging.getLogger(name)
    _logger.setLevel(_logging.INFO)

    # Remove all existing handlers
    _logger.handlers.clear()

    # Create a RichHandler that uses the shared console
    # This ensures log messages coordinate with progress bar display
    handler = RichHandler(
        console=_console,
        show_time=True,
        show_path=False,
        rich_tracebacks=True,
        markup=True,
    )

    # Add handler to logger
    _logger.addHandler(handler)

    # Only log at first logger
    _logger.propagate = False

    # Also set the root logger's handlers to use rich
    _logging.root.handlers.clear()
    _logging.root.addHandler(handler)
    _logging.root.setLevel(_logging.INFO)

    # Define error and critical as print + exit
    def error(message):
        """
        Log an error and terminate the process.

        Parameters
        ----------
        message : str
            Error message to log before exiting.

        Raises
        ------
        SystemExit
            Always raised after logging the message.
        """
        _logger.error(message)
        raise RuntimeError(message)

    def critical(message):
        """
        Log a critical error and terminate the process.

        Parameters
        ----------
        message : str
            Critical message to log before exiting.

        Raises
        ------
        SystemExit
            Always raised after logging the message.
        """
        _logger.critical(message)
        raise RuntimeError(message)

    return (_logger.debug, _logger.info, _logger.warning, error, critical)


debug, info, warning, error, critical = _init_logging("dtcc-common")


def init_logging(name="dtcc-core"):
    "Initialize logging for given package"
    return _init_logging(name)


def get_logger(name="dtcc-core"):
    "Get logger for given package"
    if name not in loggers:
        loggers[name] = _init_logging(name)
    return loggers[name]


def set_log_level(level):
    """Set log level. Valid levels are:

    "DEBUG"
    "INFO"
    "WARNING"
    "ERROR"
    "CRITICAL"

    """
    global _logger
    if _logger is None:
        _init_logging("dtcc-core")
    _logger.setLevel(level)
