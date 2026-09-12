"""Validate a canonical native model with an explicitly selected local schema.

Exit 0 for valid, 1 for semantic violations, 2 for configuration/native failures.
Requires DTCC and the optional experiment requirements in the same environment.
"""

import argparse
from dataclasses import asdict
import json
from pathlib import Path

from dtcc_core import io
from dtcc_core.model.profiles import SemanticProfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('model', type=Path)
    parser.add_argument('schema', type=Path)
    args = parser.parse_args()
    try:
        profile = SemanticProfile(args.schema)
        report = profile.validate(io.load_model(args.model))
    except (ImportError, OSError, ValueError, TypeError, NotImplementedError) as exc:
        parser.exit(2, f'Validation could not run: {exc}\n')
    print(json.dumps({'valid': report.valid, **asdict(report)}, indent=2))
    return 0 if report.valid else 1


if __name__ == '__main__':
    raise SystemExit(main())
