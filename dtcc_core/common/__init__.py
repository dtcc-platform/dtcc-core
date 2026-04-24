from .dtcc_logging import (
    init_logging,
    get_logger,
    get_python_logger,
    make_table,
    log_renderable,
    log_table,
)

debug, info, warning, error, critical = get_logger()

__all__ = [
    "init_logging",
    "get_logger",
    "get_python_logger",
    "make_table",
    "log_renderable",
    "log_table",
    "debug",
    "info",
    "warning",
    "error",
    "critical",
]
