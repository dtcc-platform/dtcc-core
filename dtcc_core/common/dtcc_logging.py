# Copyright(C) 2023 Anders Logg
# Licensed under the MIT License

import logging as _logging
import sys

# Global logger dictionary
loggers = {}

# Global logger object
_logger = None


def _init_logging(name):
    "Internal function for initializing logging"

    global _logger

    # Set log format
    format = "%(asctime)s [%(name)s] [%(levelname)s] %(message)s"

    # Initialize logger
    _logger = _logging.getLogger(name)
    _logger.setLevel(_logging.INFO)

    # Remove all existing handlers
    _logger.handlers.clear()

    # Create a new handler explicitly using stdout
    handler = _logging.StreamHandler(sys.stdout)
    handler.setFormatter(_logging.Formatter(format))

    # Add handler to logger
    _logger.addHandler(handler)

    # Only log at first logger
    _logger.propagate = False

    # Also set the root logger's handlers to stdout to override any previous settings
    _logging.root.handlers.clear()
    _logging.root.addHandler(handler)
    _logging.root.setLevel(_logging.INFO)

    # Define error and critical as print + exit
    def error(message):
        _logger.error(message)
        exit(1)

    def critical(message):
        _logger.critical(message)
        exit(1)

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
