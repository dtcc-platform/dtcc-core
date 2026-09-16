# 3DBAG attributes as qualified DTCC values

Current schema 0.8.0 retains this mapping and adds [strict triangular TINRelief](strict-tin-relief.md).
The 0.7.0 schema and performance evidence below are archived milestone evidence.

Implemented 12 September 2026. Standard [schema 0.7.0](https://github.com/dtcc-platform/dtcc-core/blob/349cdfa8ae3eabceb65a7131aa6a1e8c847c0e46/schemas/archive/model/0.7.0/schema.yaml)
adds `ElevationMeasurement` and optional building/part `elevation_measurements`.
The explicit [load_3dbag loader](../../dtcc_core/io/three_d_bag.py) implements the
audited `b3_h_dak` naming convention. Ordinary `load_city` preserves supplier
attributes without applying this interpretation. No Python model class, wire
revision, dependency or second schema validator is introduced.

## Source and interpretation boundary

The evidence is the unchanged [3DBAG tutorial tile](https://3d.bk.tudelft.nl/opendata/cityjson/3dcities/v2.0/9-284-556.city.json),
SHA256 `2bb5d22ae2cbfe2096041e3a79b3e826f43d3a8c9824eeb019feab4e0a2742ba`:
1,110 Buildings, 1,111 BuildingParts and 4,443 representations. Its metadata gives
EPSG:7415 but no 3DBAG release. Source attributes occur on the Buildings; they are
not copied to Parts or individual roof regions. Attribute ownership alone does not
establish which roof patches contributed to a statistic.

[3DBAG's layer definitions](https://docs.3dbag.nl/en/schema/layers/#calculation-of-height-values)
describe the old roof names as statistics from reconstructed LoD2.2 roofs and the
ground value as the fifth percentile of nearby AHN ground points (4 m radius).
These values are NAP elevations. A roof elevation is not a roof-to-ground distance.
The loader preserves source reference names rather than inventing a roof-point ID
or deriving a scalar building height.

The layer page's EPSG:7408 citation is not used: the local PROJ/EPSG database resolves
EPSG:7415 into EPSG:28992 and [NAP height EPSG:5709](https://www.opengis.net/def/crs/EPSG/0/5709).
The loader requires EPSG:7415 (including equivalent CRS spellings), records the
vertical CRS explicitly and performs no coordinate, datum or unit conversion.
Negative values relative to NAP are valid elevations.

Current [attribute documentation](https://docs.3dbag.nl/en/schema/attributes/)
uses newer `b3_h_50p`, `b3_h_70p`, `b3_h_min`, `b3_h_max` names. This loader rejects
those names rather than treating them as equivalent without another mapping audit.
The [release notes](https://docs.3dbag.nl/en/overview/release_notes/) document changing
attributes and a generated attribute schema introduced in 2025. The mapping ID
below versions DTCC's interpretation, not an inferred 3DBAG release.

## Mapping decisions

| Source attribute | Meaning used by this mapping | DTCC result |
| --- | --- | --- |
| `b3_h_maaiveld` | Ground elevation statistic, metres relative to NAP | Signed elevation record, reference code equal to source attribute name |
| `b3_h_dak_min` | Minimum reconstructed roof elevation | Same carrier; source-specific reference retained |
| `b3_h_dak_50p` | Median reconstructed roof elevation | Same carrier; no selection as representative building height |
| `b3_h_dak_70p` | 70th percentile reconstructed roof elevation | Same carrier; no implicit extrusion |
| `b3_h_dak_max` | Maximum reconstructed roof elevation | Same carrier; no fabricated ridge measurement |
| `b3_dak_type` | Source roof categorization | `horizontal`, `multiple horizontal`, `slanted` become qualified `roof_type` codes in the mapping namespace |
| `b3_dak_type` failure/unknown/other codes | Includes `unknown`, `no points`, `no planes` | Retain raw source value only; do not claim a physical roof class |
| `b3_pw_bron`, `b3_pw_datum`, `b3_pw_selectie_reden` | Point-cloud lineage | Retain exactly on source owner, without assigning capture status or date to a derived statistic |
| `b3_val3dity_lod12/13/22`, `b3_reconstructie_onvolledig` | Supplier reconstruction/validation evidence | Preserve; not DTCC geometric certification |
| `b3_rmse_lod12/13/22`, `b3_volume_lod12/13/22` | Supplier fit/volume quantities | Preserve; qualified error/volume semantics not introduced in this slice |
| `b3_nodata_fractie_ahn3/4`, `b3_nodata_radius_ahn3/4`, `b3_puntdichtheid_ahn3/4`, `b3_mutatie_ahn3_ahn4`, `b3_kas_warenhuis` | Source coverage/selection/reconstruction metadata | Preserve without interpreting them as building function, confidence or measured height |
| `identificatie`, `oorspronkelijkbouwjaar`, `status`, `geconstateerd`, document/validity/registration fields | BAG identity and administrative history | Preserve spelling, precision, nulls and values; no replacement of object ID, invented dates or state conversion |

The cached `b3_pw_datum` values are integer years. Current documentation describes
year precision, using point timestamps or a file-date fallback. Neither an exact
acquisition instant nor a universal survey method can be inferred from it. No
`status: measured` is assigned to reconstructed roof statistics. The original
metadata is the retained provenance evidence; the new record's `source` points
to the source attribute on that same object.

All mapped source numbers must be finite numbers, excluding booleans. Missing or
null numbers add no record; nothing becomes zero. Unrecognized roof strings stay
raw. A non-string roof code fails. Existing destination keys fail before returning
a result, even if equal, so an enriched file must be loaded with ordinary
`load_city`. Source-only imports can be repeated independently. There is no implicit
merge or synchronization after editing canonical attributes.

## Versioned record and Python use

`ElevationMeasurement` is a closed, identifier-free dictionary requiring `value`,
`unit`, `vertical_reference`, and `reference`; optional `source` carries provenance.
The value is signed; units reuse `m`, `cm`, `mm`, `km`. `reference` accepts a nonblank
string or existing `QualifiedCode`. A vertical reference is a nonblank vertical CRS
identifier; schema validation checks its presence/type, not geodetic correctness
or remote registration. Authoring other vertical references is possible without
claiming that DTCC can transform them.

The mapping identifier is
`https://github.com/dtcc-platform/dtcc-core/mappings/3dbag/b3_h_dak/1`.
Its `/elevation-references` code space contains exactly the five source names above;
its `/roof-types` code space contains the three physical codes above. These are
local logical identifiers defined by this document, not endpoints fetched during
I/O or an official CityGML code list. The loader enforces its mapping vocabulary;
the generic schema enforces qualified record shape.

```python
from dtcc_core import io
from dtcc_core.datasets import load_model_package

city = io.load_3dbag("3dbag.city.json", extent_policy="recompute")
building = city.buildings[0]
ground = next(r for r in building.attributes["elevation_measurements"]
              if r["reference"]["value"] == "b3_h_maaiveld")
print(ground["value"], ground["unit"], ground["vertical_reference"])
print(building.height)  # None in the audited tile
print(building.attributes["roof_type"]["value"])
city.save("mapped.dtcc")
restored = io.load_city("mapped.dtcc")
package = city.export("mapped.dtccpkg", canonical=True)
with_context = load_model_package(package.path)
```

`extent_policy="recompute"` is explicit for the tile's 1,111 inconsistent extent
summaries; unchanged geometry and source summaries remain evidence. Loader defaults
still reject discrepancies. Canonical packages retain Dataset Context, including
the mapping step and extent evidence. Standalone `.dtcc` preserves record provenance
and supplier attributes, but does not carry root Dataset Context. Strict CityJSON
export preserves enriched attributes; this is DTCC-enriched CityJSON, not an assertion
of conformance to a particular 3DBAG release schema. Reimport selects the current
DTCC schema; CityJSON does not carry the DTCC schema selection.

Default schema validation applies after mapping and at native/package/strict
CityJSON boundaries. `validate_schema=False` skips semantic evaluation only;
source admission, source CRS and mapping interpretation checks still apply.
Schema 0.6.0 is archived unchanged, outside the runtime bundle. No wire migration
is required: current `.dtcc` remains an ordinary `DTCC.ModelFile` Protobuf v6.

## Evidence

Run `python -m sandbox.model_profiles.three_d_bag_example SOURCE OUTPUT` from the
repository with the pinned source. The example compares every native source fact
after removing only enrichment, exercises native/package/strict CityJSON, preserves
package context, and checks that removing a required vertical reference cannot
overwrite an existing file. It measures warm schema evaluation separately from
default serialization. See the active plan's completion section for observed
results and timings; no cross-platform or full geometric validity claim is made.

Observed on the development Mac, with schema initialization excluded:

| Operation | Warm median |
| --- | --- |
| Schema evaluation, unchanged tile | 0.289 s |
| Schema evaluation, tile plus 5,550 elevations and 1,110 roof codes | 0.693 s |
| Enriched native encode, default validation | 1.932 s |
| Enriched native decode, default validation | 2.916 s |

Schema comparisons use five alternating evaluations; codecs use three runs.
The explicit loader took 3.418 s with schema already initialized (one run).
Enriched native size is 26,854,163 bytes, versus 24,504,572 before enrichment.
The increase carries the new qualified metadata while retaining the original facts.
No material schema-cost change is claimed from adding an unused class; the extra
0.403 s here is evaluation of additional records on this workload.

Native and canonical package comparisons are exact. CityJSON preserves IDs,
attributes, representations, shells and region assignments; coordinates differ
by at most 0.000375 m at the writer's 0.001 m quantization, within the existing
0.00050001 m comparison tolerance. CityJSON is not byte-exact native persistence.
Artifacts/report: `/private/tmp/dtcc-3dbag-mapped`; regression log:
`/private/tmp/dtcc-3dbag-regression.log` (686 passed, 1 skipped).
