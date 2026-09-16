"""Exercise the checker CLI against real coverage from a small public API."""

import os
from pathlib import Path
import subprocess
import sys
import textwrap

import pytest


CHECKER = Path(__file__).resolve().parents[1] / "scripts/check_public_api_calls.py"


def check_api(tmp_path, source, invocation):
    package = tmp_path / "sample_api"
    package.mkdir()
    (package / "__init__.py").write_text(textwrap.dedent(source))
    runner = tmp_path / "runner.py"
    runner.write_text("import sample_api\n" + invocation)
    env = dict(os.environ, COVERAGE_FILE=str(tmp_path / ".coverage"))
    # Override any pytest-cov subprocess configuration from the parent suite.
    for name in list(env):
        if name.startswith("COV_CORE_") or name == "COVERAGE_PROCESS_START":
            del env[name]
    for args in (
        ["-m", "coverage", "run", "--source=sample_api", str(runner)],
        ["-m", "coverage", "json", "-o", "coverage.json"],
    ):
        subprocess.run(
            [sys.executable, *args], cwd=tmp_path, env=env,
            check=True, capture_output=True, text=True,
        )
    env["PYTHONPATH"] = str(tmp_path)
    return subprocess.run(
        [sys.executable, str(CHECKER), "--strict", "--package", "sample_api",
         "--coverage-file", "coverage.json"],
        cwd=tmp_path, env=env, capture_output=True, text=True,
    )


@pytest.mark.parametrize("call_bodies", [False, True])
def test_import_does_not_count_as_body_execution(tmp_path, call_bodies):
    result = check_api(tmp_path, '''
        from functools import wraps

        def decorate(fn):
            @wraps(fn)
            def wrapper(*args, **kwargs):
                return fn(*args, **kwargs)
            return wrapper

        __all__ = ["plain", "decorated", "multiline"]

        def plain():
            """A docstring is not evidence of a call."""
            return 1

        @decorate
        def decorated():
            return 2

        def multiline(
            value=int(
                "3"
            ),
        ):
            return value
    ''', "sample_api.plain(); sample_api.decorated(); sample_api.multiline()\n"
        if call_bodies else "")
    assert result.returncode == (0 if call_bodies else 1), result.stdout + result.stderr
    assert f"Functions covered: {3 if call_bodies else 0}" in result.stdout


@pytest.mark.parametrize("call_bodies", [False, True])
def test_ambiguous_definition_lines_are_not_call_evidence(tmp_path, call_bodies):
    result = check_api(tmp_path, '''
        __all__ = ["one_line", "docstring_only"]
        def one_line(): return 1
        def docstring_only():
            """This function has no separate executable body line."""
    ''', "sample_api.one_line(); sample_api.docstring_only()\n" if call_bodies else "")
    assert result.returncode == 1, result.stdout + result.stderr
    assert "Functions covered: 0" in result.stdout
    assert "Functions missed: 2" in result.stdout
