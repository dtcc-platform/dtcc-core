"""Local, versioned schema selection for the canonical exchange boundary."""

from functools import lru_cache
from importlib.resources import files
import re

SCHEMA_ID = 'https://github.com/dtcc-platform/dtcc-core/schemas/model'
SEMANTIC_NAMESPACE = SCHEMA_ID + '#'
DEFAULT_VERSION = '0.9.0'


@lru_cache(maxsize=8)
def _profile(schema_id, version):
    from .profiles import SemanticProfile

    if schema_id != SCHEMA_ID or re.fullmatch(r'[0-9]+\.[0-9]+\.[0-9]+', version) is None:
        raise ValueError(f'Unsupported semantic schema {schema_id!r} version {version!r}')
    path = files('dtcc_core').joinpath('schemas', 'model', version, 'schema.yaml')
    if not path.is_file():
        raise ValueError(f'Unsupported semantic schema {schema_id!r} version {version!r}')
    profile = SemanticProfile(str(path))
    if (profile.profile_id, profile.profile_version) != (schema_id, version):
        raise ValueError('Bundled schema declaration does not match its selected identity/version')
    return profile


def validate_admitted(model, schema_id, version):
    """Evaluate semantics after canonical admission, never instead of admission."""
    report = _profile(schema_id, version)._validate_admitted(model)
    if not report.valid:
        details = '\n'.join(f'{issue.path}: {issue.message}' for issue in report.issues[:10])
        if len(report.issues) > 10:
            details += f'\n... and {len(report.issues) - 10} more issues'
        raise ValueError(f'Schema validation failed ({schema_id}, {version}):\n{details}')
