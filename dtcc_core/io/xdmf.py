"""
XDMF template for writing a tetrahedral volume mesh with boundary faces and markers.

The template expects three HDF5 datasets:
- Mesh geometry at `{h5file}:/Mesh/mesh/geometry` (n_pts x 3)
- Mesh topology at `{h5file}:/Mesh/mesh/topology` (n_tets x 4)
- Boundary face topology and markers at `{h5file}:/MeshTags/boundary_markers/topology` (n_facets x 3)
  and `{h5file}:/MeshTags/boundary_markers/values` (n_facets)

Placeholders:
- `h5file`: Path to the HDF5 file (as referenced from the XDMF)
- `n_pts`: Number of points
- `n_tets`: Number of tetrahedra
- `n_facets`: Number of boundary facets
"""


XDMF_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<Xdmf Version="3.0" xmlns:xi="https://www.w3.org/2001/XInclude">
  <Domain>

    <!-- 1) Volume mesh, no attributes here -->
    <Grid Name="mesh" GridType="Uniform">
      <Topology TopologyType="Tetrahedron"
                NumberOfElements="{n_tets}" NodesPerElement="4">
        <DataItem Format="HDF" NumberType="Int" Dimensions="{n_tets} 4">
          {h5file}:/Mesh/mesh/topology
        </DataItem>
      </Topology>
      <Geometry GeometryType="XYZ">
        <DataItem Format="HDF" NumberType="Float" Dimensions="{n_pts} 3">
          {h5file}:/Mesh/mesh/geometry
        </DataItem>
      </Geometry>
    </Grid>

    <!-- 2) Facet markers on a separate grid -->
    <Grid Name="boundary_markers" GridType="Uniform">
      <!-- re-use the same points -->
      <xi:include xpointer="xpointer(/Xdmf/Domain/Grid/Geometry)"/>

      <Topology TopologyType="Triangle"
                NumberOfElements="{n_facets}" NodesPerElement="3">
        <DataItem Format="HDF" NumberType="Int" Dimensions="{n_facets} 3">
          {h5file}:/MeshTags/boundary_markers/topology
        </DataItem>
      </Topology>

      <Attribute Name="boundary_markers"
                 AttributeType="Scalar"
                 Center="Cell">
        <DataItem Format="HDF"
                  NumberType="Int"
                  Dimensions="{n_facets}">
          {h5file}:/MeshTags/boundary_markers/values
        </DataItem>
      </Attribute>
    </Grid>

  </Domain>
</Xdmf>
"""
