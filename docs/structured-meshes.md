# Structured meshes

The `RectangleMesh` and `BoxMesh` factory functions create native meshes
without FEniCSx or an external mesher:

```python
from dtcc_core import BoxMesh, RectangleMesh

surface = RectangleMesh(8, 12)  # Mesh on the unit square, z=0
volume = BoxMesh(8, 12, 6)      # VolumeMesh on the unit cube
```

Arguments count subdivisions, not vertices, and must be positive integers.
Coordinates are local metres. Each rectangle has two counter-clockwise
triangles; each box has six positively oriented tetrahedra. Neighbouring boxes
share matching facets. These are ordinary DTCC model objects, supporting model
inspection, transformations and canonical serialization.

`BoxMesh` includes outward-oriented boundary triangles with labels xmin=-3,
xmax=-4, ymin=-5, ymax=-6, zmin=-1 and zmax=-2. `RectangleMesh` is a triangular
surface mesh; it does not assign volume boundary labels.

The [tagged facets guide](volume-mesh-boundaries.md) describes how `.dtcc`,
`.dtccpkg` and XDMF preserve boundary labels. VTK/VTU export and `to_cpp()`
omit them.
