import requests
from urllib3.exceptions import MaxRetryError, NewConnectionError, ReadTimeoutError

from dtcc_core.io.data.timeouts import is_read_timeout


def test_read_timeouts_wrapped_in_connection_errors_are_recognised():
    stalled = ReadTimeoutError(None, "/", "Read timed out.")
    assert is_read_timeout(requests.ConnectionError(stalled))
    assert is_read_timeout(
        requests.ConnectionError(MaxRetryError(None, "/", reason=stalled))
    )


def test_other_connection_errors_are_not_read_timeouts():
    refused = NewConnectionError(None, "Connection refused")
    assert not is_read_timeout(
        requests.ConnectionError(MaxRetryError(None, "/", reason=refused))
    )
    assert not is_read_timeout(requests.ConnectionError("down"))
    assert not is_read_timeout(requests.ConnectionError())
    assert not is_read_timeout(requests.ConnectionError(None))
