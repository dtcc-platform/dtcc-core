"""Scratch probe: fails on macOS and Windows on the first run attempt only."""

import os
import platform


def test_rerun_probe():
    first_attempt = os.environ.get("GITHUB_RUN_ATTEMPT") == "1"
    assert not (first_attempt and platform.system() != "Linux")
