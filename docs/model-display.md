# Model display and inspection

DTCC provides two presentations: a compact representation for ordinary Python
display, and an explicit detailed report for inspection.

```python
from dtcc_core.model import Point, Mesh

point = Point(x=1.0, y=2.0, z=3.0)
print(point)       # Point(x=1.0, y=2.0, z=3.0)
print(Mesh())      # <Mesh(num_vertices=0, num_faces=0, num_fields=0)>
```

`repr(obj)`, `str(obj)`, `print(obj)` and ordinary string interpolation use the
same compact representation. It includes the actual class name and named scalar
facts such as identity, counts, shape or dtype. It does not expand arrays or
children, calculate bounds or statistics, or access providers. Long identifiers
are abbreviated. A list of models therefore stays readable too.

Angle brackets mark a **summary**, not an executable expression. Only complete
finite `Point` and `Bounds` values opt into constructor expressions without
brackets. A point with fields, regions, a nondefault transform, dataset context or
schema declarations is a summary too. Derived bounds caches are not part of a
point's value. Constructor expressions assume the relevant DTCC classes are in
scope; they preserve coordinate values rather than NumPy scalar storage types.
Use `.copy()` for copying and the native I/O APIs for persistence, rather than
using `eval()` as a general reconstruction mechanism.

## Detailed reports

```python
mesh.info()                         # Print a detailed report; return None
text = mesh.info(print=False)        # Return the same plain text, without printing
mesh.info(presentation=False)        # Omit attached dataset context
with open("mesh-info.txt", "w") as f:
    mesh.print_info(file=f)          # Existing file-output convenience
```

Reports use the shared DTCC table styling for properties, fields, attributes,
geometry representations, grouped children and domain statistics. Explanations
and configuration guidance use ordinary text. Model tables show up to 20 rows
and state how many were omitted; sensor and vehicle samples show the first three.
Dataset parameter help, catalogue listings and context tables remain complete.
Inspection summarizes array shapes and types rather than dumping their values.
Use the named data attributes directly when you need the full arrays or tables.

Detailed inspection may calculate bounds and existing domain statistics. Bounds
follow each model's existing coordinate and cache semantics: public array edits
may require `calculate_bounds()` to refresh them. Geometry bounds are local;
object envelopes do not apply transforms. Inspection does not validate geometry
or certify compatible coordinate frames.

Dataset-produced models append metadata, presentation and provenance tables by
default. All reports are plain text without ANSI codes, treat values as literal
text rather than Rich markup, and remain visible regardless of logging level.
`.tree()` remains the explicit recursive hierarchy inspection method. File
inspection functions returning metadata dictionaries and logging's `info()` keep
their separate purposes.

## Dataset definitions and catalogue

```python
from dtcc_core import datasets

print(datasets.smoke)                # Compact descriptor
datasets.smoke.info()                # Description and parameter table
help_text = datasets.smoke.info(print=False)
datasets.info()                      # Grouped dataset catalogue
catalogue = datasets.info(print=False)
```

Parameter help uses `show_options()`, including the discovered schema of remote
datasets, without contacting their service. `describe()` remains the structured
metadata API. Code that previously used `str(dataset)` for parameter help or
`str(sensors)` / `str(vehicles)` for detailed reports should use
`.info(print=False)` instead. Display text is intended for people, not parsing.

## Extending native models

Model dataclasses use `@dataclass(repr=False)` so Python does not replace the
shared representation with a recursive dump. `_summary_items()` supplies cheap
`(name, value)` pairs. `_info_sections()` supplies detailed
`(heading, columns, rows)` sections; a section with `columns=None` supplies prose.
The shared `Model.info()` owns rendering, context attachment and print/return
behavior. These hooks are private implementation conventions, not a public
report schema. New subclasses should retain the default summary markers unless
they deliberately implement and verify complete constructor expressions.

Run the complete offline examples from the repository:

```sh
.venv/bin/python sandbox/model_display.py
```
