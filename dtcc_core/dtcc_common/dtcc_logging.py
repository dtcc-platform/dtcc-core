# Copyright(C) 2023 Anders Logg
# Licensed under the MIT License

import logging as _logging

# Global logger object
_logger = None

def _init_logging(name):
    "Internal function for initializing logging"

    global _logger

    # Initialize logger
    format = "%(asctime)s [%(name)s] [%(levelname)s] %(message)s"
    _logging.basicConfig(format=format)
    _logger = _logging.getLogger(name)
    _logger.setLevel(_logging.INFO)

    # Define error and critical as print + exit
    def error(message):
        _logger.error(message)
        exit(1)

    def critical(message):
        _logger.critical(message)
        exit(1)

    return (_logger.debug, _logger.info, _logger.warning, error, critical)


debug, info, warning, error, critical = _init_logging("dtcc-common")


def init_logging(name):
    "Initialize logging for given package"
    debug(f"Initializing logging for {name}")
    return _init_logging(name)


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
        _init_logging("dtcc-common")
    _logger.setLevel(level)
