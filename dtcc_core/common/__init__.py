from .dtcc_logging import init_logging, get_logger

debug, info, warning, error, critical = get_logger()

__all__ = [
    "init_logging",
    "get_logger",
    "debug",
    "info",
    "warning",
    "error",
    "critical",
]
