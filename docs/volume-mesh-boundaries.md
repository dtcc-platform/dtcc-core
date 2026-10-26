# Tagged volume-mesh facets

`VolumeMesh.boundary_faces` is an integer (N, 3) array of vertex indices;
`boundary_markers` is an aligned integer (N,) array. A mapping from three vertex
indices to an integer marker is also accepted. Native loading normalizes either
input to aligned arrays. Cell markers remain separate.

Canonical `.dtcc` and `.dtccpkg` artifacts preserve these arrays, including
negative marker values. XDMF/HDF5 stores the same facet-label pairs. Every
facet must contain three distinct, in-range vertices; duplicate facets,
noninteger labels and mismatched lengths are rejected before a native file is
replaced. This validates the representation, not whether a tagged face is an
exterior boundary or whether a particular PDE boundary condition is appropriate.

The wire schema adds optional fields 4 and 5 to VolumeMesh. Untagged artifacts
omit both fields and remain readable. Tagged artifacts require this reader;
older strict Core readers reject their unknown fields. Regenerate Python
bindings from `dtcc_core/schemas/dtcc.proto` with a Protobuf 5.29-compatible
compiler. Do not edit generated bindings manually.

```python
mesh.boundary_markers = {(0, 1, 2): -2, (0, 1, 3): 17}
mesh.save("mesh.dtcc")
restored = dtcc_core.io.load_volume_mesh("mesh.dtcc")
labels = dict(zip(map(tuple, restored.boundary_faces), restored.boundary_markers))
```

A dataset package additionally retains its DatasetContext. Native mesh files
remain standalone artifacts, rather than numerical restart checkpoints.
