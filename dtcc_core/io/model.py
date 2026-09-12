"""File I/O for DTCC Protobuf model artifacts (.dtcc)."""

import os
from pathlib import Path
import tempfile

from ..model import exchange


def save_model(model, path, *, validate_schema=True):
    """Atomically save the supported canonical model subset to a .dtcc file.

    The message structure is defined by dtcc_core/proto/dtcc.proto.
    Invalid or unsupported state fails before replacing an existing file.
    Semantic schema validation runs by default; False bypasses only that stage.
    """
    path = Path(path)
    if path.suffix.lower() != '.dtcc':
        raise ValueError("Canonical model files use the .dtcc extension")
    data = exchange.dumps(model, validate_schema=validate_schema)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def load_model(path, *, expected_type=None, validate_schema=True):
    """Load a canonical artifact using its own type/version discriminator."""
    with Path(path).open('rb') as stream:
        data = stream.read(exchange.MAX_BYTES + 1)
    model = exchange.loads(data, validate_schema=validate_schema)
    if expected_type is not None and type(model) is not expected_type:
        raise ValueError(f"Expected {expected_type.__name__}, artifact contains {type(model).__name__}")
    return model
