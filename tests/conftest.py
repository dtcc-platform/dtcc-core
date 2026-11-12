"""Shared pytest fixtures."""

from __future__ import annotations

import socket

import pytest


@pytest.fixture
def no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    """Make any attempt to open a socket connection raise."""

    def deny(*_args, **_kwargs):
        raise RuntimeError("network disabled")

    monkeypatch.setattr(socket.socket, "connect", deny)
    monkeypatch.setattr(socket, "create_connection", deny)
