"""Recognise read timeouts that requests reports as connection errors."""

from urllib3.exceptions import MaxRetryError, ReadTimeoutError


def is_read_timeout(exc: BaseException) -> bool:
    """Whether ``exc`` is a read timeout that requests raised as ``ConnectionError``.

    requests raises ``ReadTimeout`` only for some read timeouts. With no read
    retries left, urllib3 raises ``MaxRetryError`` with a ``ReadTimeoutError``
    reason, and a server that stalls while sending the body raises a bare
    ``ReadTimeoutError``; requests wraps both in ``ConnectionError``.
    """
    reason = exc.args[0] if exc.args else None
    if isinstance(reason, MaxRetryError):
        reason = reason.reason
    return isinstance(reason, ReadTimeoutError)
