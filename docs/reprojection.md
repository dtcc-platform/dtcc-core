# Reprojecting geometry with fields

Reprojection changes horizontal coordinates while keeping Z, connectivity,
markers and attached field values unchanged. Field names, descriptions, units,
dimensions and associations are preserved. Values remain attached to the same
vertices, faces, cells or samples; no interpolation or resampling is performed.

For example, a simulation temperature stays 20 degrees Celsius when its node's
location is expressed in longitude and latitude instead of projected metres:

```python
import numpy as np
from dtcc_core.model import Field, Mesh
from dtcc_core.reproject.reproject import reproject_mesh

mesh = Mesh(
    vertices=np.array([
        [500000., 6500000., 100.],
        [500010., 6500000., 100.],
        [500000., 6500010., 100.],
    ]),
    faces=np.array([[0, 1, 2]]),
    fields=[Field(
        name="temperature", unit="degC", association="vertex",
        values=np.array([20., 21., 22.]),
    )],
)
mapped = reproject_mesh(mesh, "EPSG:3006", "EPSG:4326")
np.testing.assert_array_equal(mapped.fields[0].values, mesh.fields[0].values)
```

`reproject_mesh` supports both `Mesh` and `VolumeMesh`. The same preservation
rule applies to `reproject_pointcloud`, `reproject_surface`,
`reproject_multisurface` (including constituent surfaces), and supported geometry
representations passed through `reproject_object`.

A CRS change returns an independent copy, updates its declared CRS and clears
derived bounds. The original geometry and fields are unchanged. A same-CRS call
returns the original object and emits no vector warning.

## Vector components and limitations

Fields with `dim > 1` emit a `UserWarning` on a CRS change:

> Vector field components are preserved unchanged. Their orientation and units
> are not transformed to the target coordinate system.

For example, velocity components in m/s stay in m/s with their original component
convention. They do not become angular rates when the target CRS uses degrees.
Keeping components can be suitable for local visualization where axis orientation
and scale closely agree; callers requiring a different vector basis must transform
the components separately. Scalar coordinate components and coordinate-dependent
statistics also remain unchanged, without an automatic warning. Core does not
infer physical meaning from field names or units.

This operation does not rotate vectors, rescale densities, transform vertical
datums or guarantee a mesh suitable for further simulation in the target CRS.
It supports horizontal two-axis CRSs only. Coordinates must already be in the
declared source frame. Existing restrictions on local affine transforms, stored
normals, semantic regions, DatasetContext and nested object frames still apply;
these are rejected rather than silently altered or discarded.
