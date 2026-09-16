"""Public logging settings and critical-error behavior."""

import logging

import pytest

import dtcc_core
from dtcc_core import common
from dtcc_core.common import dtcc_logging


def test_global_log_level_and_critical_error(monkeypatch, caplog):
    original = dtcc_logging._global_log_level
    logger = common.get_python_logger()
    monkeypatch.setattr(logger, "handlers", [caplog.handler])
    try:
        dtcc_core.set_log_level("ERROR")
        assert common.get_python_logger() is logger
        assert logger.level == logging.ERROR
        with pytest.raises(ValueError, match="Invalid log level"):
            dtcc_core.set_log_level("invalid")
        assert logger.level == logging.ERROR
        with pytest.raises(RuntimeError, match="critical test message"):
            common.critical("critical test message")
        assert [(r.levelno, r.message) for r in caplog.records] == [
            (logging.CRITICAL, "critical test message")
        ]
    finally:
        dtcc_core.set_log_level(original)
