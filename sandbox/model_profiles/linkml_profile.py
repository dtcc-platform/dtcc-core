"""Load the shared backend without importing Core in the isolated tooling process."""

import importlib.util
from pathlib import Path

_path = Path(__file__).resolve().parents[2] / 'dtcc_core/model/_profile_backend.py'
_spec = importlib.util.spec_from_file_location('_dtcc_profile_backend', _path)
_backend = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_backend)

LinkMLProfile = _backend.LinkMLProfile
graph_checks = _backend.graph_checks
issue = _backend.issue
payload = _backend.payload
