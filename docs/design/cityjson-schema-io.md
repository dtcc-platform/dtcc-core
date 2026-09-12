# Standard schema at strict CityJSON boundaries

Strict CityJSON import and export now evaluate the selected standard DTCC schema,
using the same semantic evaluator as canonical .dtcc I/O. The boundary was introduced with schema 0.5.0; the current schema is
[0.9.0](../../dtcc_core/schemas/model/0.9.0/schema.yaml), including
[qualified records](qualified-values.md). The wire layout is unchanged.
The strict subset now also includes [triangular TINRelief](strict-tin-relief.md),
mapped explicitly to the existing Terrain and Mesh classes.

```python
from dtcc_core import io

city = io.load_city("buildings.city.json", strict=True)
city.save("buildings.dtcc")
city.save("buildings-out.city.json", strict=True)

unfinished = io.load_city("unfinished.city.json", strict=True, validate_schema=False)
unfinished.save("unfinished-out.city.json", strict=True, validate_schema=False)
```

The same options work for .json.zip, City.save_cityjson, io.load_cityjson and the
direct cityjson.load dictionary API, write_cityjson.to_cityjson and
write_cityjson.save. Public dispatch delegates to
those same implementations rather than maintaining a second JSON/ZIP reader.

| CityJSON mode | Omitted option (or None) | validate_schema=True | validate_schema=False |
| --- | --- | --- | --- |
| strict=True | Evaluate standard schema | Evaluate standard schema | Skip semantics only |
| Permissive default | Existing permissive behavior | Error: requires strict=True | Existing permissive behavior |

A string such as "false" and numeric flags such as 0 are rejected. The None default
is specific to CityJSON's mode-dependent behavior; canonical I/O still requires a
literal boolean. Schema evaluation never runs on ordinary attribute access.

## Admission and selection

On strict import, source-format and native integrity checks complete first. The
semantic evaluator then checks declared building attributes, counts, surface owners,
opening hosts and other known rules. The returned City records the selected schema
ID/version even when semantic evaluation was bypassed. A declaration identifies
the intended contract; it is not a validation certificate.

The supported CityJSON mapping contains no DTCC schema declaration. Import selects
the bundled default, currently 0.9.0. Strict export uses the root City's stored
schema ID/version, or the bundled default if both are absent, through the same
selection logic as canonical saving. A partial/malformed declaration fails even
with bypass; an unavailable schema fails evaluation unless explicitly bypassed.
Export does not mutate the root's selection. Stored Object profile labels do not
select another schema or bypass the existing format restrictions.

CityJSON output does not encode DTCC schema metadata. Reimport selects the default
again; .dtcc is the format that preserves the declaration. In particular, bypassing
an unavailable old schema during CityJSON export does not preserve that old schema
selection on reimport or promise that the result satisfies the current default.
There is no new extension field, schema download, migration or silent value coercion.

## Failures and bypass

Semantic failures raise ValueError with object IDs and native attribute/relationship
paths before returning an imported model or opening an export destination. Tests
cover preservation of existing plain JSON and ZIP files on semantic failures.
This is not a new crash-atomic writer: the existing physical JSON/ZIP writing path
is unchanged after conversion succeeds.

The explicit bypass still enforces strict source/native checks and unsupported
mapping failures, including duplicate JSON keys, array and region indices,
containment, dangling/cyclic region hosts, representation restrictions and invalid
schema identities. For example, a Window referring to a GroundSurface is a semantic
host-type failure that can be bypassed; an out-of-range host index is a structural
failure that cannot. Unknown metadata remains open under the standard schema.

Each strict import/conversion invokes native admission once and semantic admission
once when enabled. It calls the existing admitted-model evaluator directly, without
serializing the model for validation or re-entering standalone SemanticProfile.validate.
The compiled schema is reused; validity is evaluated against the current mutable data.

## Remaining scope

The existing faithful CityJSON subset is unchanged. Native builder representation
IDs/roles, unsupported feature kinds and other unmapped state still require explicit
mappings. Permissive loading can omit unsupported content and is not certified by
this change. Other external format adapters and default artifact-only packages
remain separate adoption work. Qualified measurements, code spaces and additional
urban themes remain schema work; this does not claim complete CityGML conformance.

See the [implementation plan](../../.agent/plans/2026-09-12-cityjson-schema.md) for
observed verification and installed-package evidence.
