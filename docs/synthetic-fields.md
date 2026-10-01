# Synthetic volume fields

Generate geometry and fields using the bounds and coordinate reference of the
saved flagship model, without network access:

```sh
uv run python scripts/generate_fields.py
```

Outputs go to `data/fields/` relative to the checkout, even when run from
`scripts/`. Use `--output-dir PATH` to choose another directory. `-nx`, `-ny`,
`-nz` (also `--nx`, `--ny`, `--nz`) count cube subdivisions along each axis.
Defaults are **128 × 128 × 36**. Edit `DEFAULT_NX`, `DEFAULT_NY`, and
`DEFAULT_NZ` at the top of the script to change them. Time and path defaults
are also declared there as uppercase constants, and CLI flags override them.
Over the Delft domain, equal cell spacing would use
`ceil(128 * 277.152 / 2000) = 18` vertical subdivisions. The default `nz=36`
doubles that count, giving approximately **15.625 × 15.625 × 7.6987 m** cells.
These are fixed defaults; each axis can be overridden independently.

The default domain is read from `data/flagship/flagship.dtcc`. For the Delft
flagship it is E **83,570–85,570**, N **445,800–447,800**, and Z
approximately **−13.650499–263.501501 m**, in **EPSG:7415** (RD New + NAP).
The flagship generator doubles its padded vertical extent upward; this script
inherits that saved extent without applying another height multiplier.
Both `.dtcc` and VTU files store these world coordinates, so they overlay the
city directly. Use `--flagship-model PATH` for another saved flagship City.
Its root transform must be identity and its three domain extents positive.
Generate the [flagship model](flagship-model.md) first if that file is missing
or predates the expanded domain;
the field script reads its bounds and CRS once and does not modify it.
There are `(nx+1)(ny+1)(nz+1)` vertices, `nx*ny*nz` grid cells and
`6*nx*ny*nz` tetrahedra. Each cube has six positively oriented tetrahedra
around a common body diagonal; adjacent cubes share face triangulations.

Each snapshot contains exactly two **vertex** fields: `velocity` is a
float32 vector with projected `(easting, northing, up)` components in m/s; `pressure` is
float32 signed gauge pressure in Pa. Both representations use identical
coordinates and sampled values. The native grid is a `VolumeGrid`; its VTU
reference uses hexahedra, preserving the grid cells and point values.

## Synthetic construction

World positions are mapped once to `[-π, π]` along each axis before sampling
the analytic formulas. This keeps the waves and vortices spread across the
city volume without feeding large projected coordinate values into the
polynomial terms. Velocity components retain their synthetic m/s values in
the projected basis; they are not multiplied by domain lengths.

The velocity starts with the proposed XY sine/polynomial terms
`u = -sin(y) + 0.1*(x² - 2xy)` and
`v = sin(x) + 0.1*(y² - 2xy)`, shifting their arguments with time and height.
Cross-coupled helical waves add vertical motion, changing flow direction and
recirculation. Three orbiting Gaussian vortex cores with different axes and
opposite rotation directions add concentrated swirls and axial jets. Pressure
combines traveling waves, positive and negative moving cores, and a velocity
magnitude contribution. The scalar and every vector component vary in space
and time, providing useful slices, isosurfaces, glyphs and streamlines.

These are analytic visualization fixtures, not a fluid solver or measurements.
The fields do not enforce incompressibility, boundary conditions or a physical
pressure/velocity relation. They are finite everywhere, without masks or NaNs.
Very coarse resolutions cannot resolve all the small vortices and waves;
use `-nx 32 -ny 32 -nz 32` for a quick preview at lower resolution.

## Time and output files

DTCC Model's `Field` stores `(N,)` scalar or `(N, dim)` vector arrays. Neither
the Python model nor the canonical protobuf schema has a native time axis or
time-series reader. The flagship encodes time using ordinary snapshot Objects
and metadata. This script writes separate geometry files per time, avoiding a
model/schema change and keeping each file directly loadable as a mesh or grid.

`-nt` counts snapshots, including `t=0` and `t=T`; `-T` (also `--duration`)
sets both the sampled duration and the field's period in seconds. Every phase
uses integer harmonics, so the last frame matches the first for a smooth loop.
For example, the defaults produce 21 evenly spaced times from 0 to 10 seconds,
with **42 `.dtcc` files** and 42 matching VTU files:

```text
fields_tet_0000.dtcc   fields_grid_0000.dtcc
fields_tet_0000.vtu    fields_grid_0000.vtu
...
fields_tet_0020.dtcc   fields_grid_0020.dtcc
fields_tet_0020.vtu    fields_grid_0020.vtu
fields_tet.pvd        fields_grid.pvd
fields-series.json
```

At the default spatial resolution, each time contains **615,717 vertices**,
**589,824 grid cells**, and **3,538,944 tetrahedra**. The native tetrahedral
arrays alone occupy approximately **77.5 MiB per snapshot**, plus **9.4 MiB**
for the native grid fields. The 21-time native series totals roughly **1.8 GiB**;
VTU references require additional storage. Use `-nt 1` to inspect one full
resolution snapshot before generating the whole series.

Use `-nt 1` for **exactly two `.dtcc` files** at `t=0`. All subdivision and
snapshot counts must be positive integers, and `T` must be finite and positive.
Resolutions whose mesh arrays exceed the canonical protobuf limit are rejected
before allocation. Geometry is built once, and only one time's field arrays
are held during generation; file size and export time still scale with the
number of snapshots. Geometry is repeated in each file. Running again replaces
matching generated filenames; older unreferenced frames may remain in the
directory, so use the manifest or PVD collections to select the current series.

`fields-series.json` records the times in seconds and relative filenames. Keep
it alongside the `.dtcc` files: bare geometry files have no timestamp property.
For example:

```python
import json
from pathlib import Path
from dtcc_core import io

directory = Path("data/fields")
series = json.loads((directory / "fields-series.json").read_text())
frame = series["frames"][5]
mesh = io.load_model(directory / frame["tet"]["dtcc"])
grid = io.load_model(directory / frame["grid"]["dtcc"])
print(frame["time_seconds"], mesh.fields, grid.fields)
```

## ParaView reference

Open `fields_tet.pvd` or `fields_grid.pvd`, then click **Apply**. The
[ParaView PVD reader](https://www.paraview.org/paraview-docs/v5.11.0/python/paraview.simple.PVDReader.html)
uses the collection's explicit times for animation. Select `pressure` or
`velocity` (magnitude or a component) for coloring, and use the animation
controls to play through the loop. A **Slice** reveals interior structure;
**Contour** on pressure, **Glyph** on velocity and **Stream Tracer** with
interior seeds offer complementary reference views. Streamlines show the flow
at the current snapshot, not particle trajectories over time. Keep color
ranges fixed across times when comparing animation frames.
