#ifndef DTCC_COLUMN_MESH_H
#define DTCC_COLUMN_MESH_H

#include "VertexSmoother.h"
#include "model/City.h"
#include "model/GridField.h"
#include "model/Mesh.h"
#include "model/Surface.h"
#include "model/Vector.h"

namespace DTCC_BUILDER
{

typedef struct ColumnIndex
{
  std::size_t column;
  std::size_t index;

  ColumnIndex() = default;

  ColumnIndex(std::size_t column, std::size_t index) : column(column), index(index) {}
} ColumnIndex;

class ColumnSimplex
{
public:
  // Vertex indices
  ColumnIndex v0{};
  ColumnIndex v1{};
  ColumnIndex v2{};
  ColumnIndex v3{};

  // Create default simplex
  ColumnSimplex() = default;

  // Create simplex and optionally sort vertices
  ColumnSimplex(ColumnIndex v0, ColumnIndex v1, ColumnIndex v2, ColumnIndex v3)
  {
    this->v0 = v0;
    this->v1 = v1;
    this->v2 = v2;
    this->v3 = v3;
  }
};

/// ColumnMesh represents a tetrahedral mesh in 3D created by
// extruding a 2D mesh in columns.
class ColumnMesh : public Printable
{
  // private:
  //   const Mesh &_ground_mesh;

public:
  /// Vector of vectors of Vertices.
  std::vector<std::vector<Vector3D>> vertices{};

  std::vector<size_t> vertices_offset;

  /// Vector of vectors of cells (tetrahedra)
  std::vector<std::vector<ColumnSimplex>> cells{};

  /// Vector of vectors of cell markers
  std::vector<std::vector<int>> markers{};

  const size_t num_cell_columns{};

  const size_t num_vertex_columns{};

  const size_t num_marker_columns{};

  ColumnMesh() = default;

  ColumnMesh(const Mesh &ground_mesh)
      : num_cell_columns(ground_mesh.faces.size()), num_vertex_columns(ground_mesh.vertices.size()),
        num_marker_columns(ground_mesh.markers.size()),
        vertices_offset(ground_mesh.vertices.size() + 1)
  {
    assert((num_vertex_columns > 0) && "Empty ground mesh. It has no faces connecting vertices");
    assert((num_cell_columns > 0) && "Empty ground mesh. It has no faces connecting faces");

    vertices.resize(num_vertex_columns);
    cells.resize(num_cell_columns);
    markers.resize(num_marker_columns);

    vertices_offset.resize(num_vertex_columns + 1, 0);
  }
  virtual ~ColumnMesh() {} // make the destructor virtual

  // Pretty-print
  std::string __str__() const override
  {
    return "ColumnMesh mesh with " + str(vertices.size()) + " vertex columns and " +
           str(cells.size()) + " cell columns";
  }

  VolumeMesh to_volume_mesh()
  {
    VolumeMesh volume_mesh;

    const size_t volume_mesh_num_vertices = vertices_offset.back() + vertices.back().size();
    // Add Vertices
    volume_mesh.vertices.reserve(volume_mesh_num_vertices);
    for (size_t j = 0; j < num_vertex_columns; j++)
    {
      for (size_t k = 0; k < vertices[j].size(); k++)
      {
        volume_mesh.vertices.push_back(vertices[j][k]);
      }
    }

    // Add Cells
    for (size_t i = 0; i < num_cell_columns; i++)
    {
      for (size_t j = 0; j < cells[i].size(); j++)
      {
        size_t vc0 = vertices_offset[cells[i][j].v0.column] + cells[i][j].v0.index;
        size_t vc1 = vertices_offset[cells[i][j].v1.column] + cells[i][j].v1.index;
        size_t vc2 = vertices_offset[cells[i][j].v2.column] + cells[i][j].v2.index;
        size_t vc3 = vertices_offset[cells[i][j].v3.column] + cells[i][j].v3.index;

        volume_mesh.cells.push_back(Simplex3D(vc0, vc1, vc2, vc3));
      }
    }

    volume_mesh.markers.reserve(volume_mesh_num_vertices);
    for (size_t j = 0; j < markers.size(); j++)
    {
      for (size_t k = 0; k < markers[j].size(); k++)
      {
        volume_mesh.markers.push_back(markers[j][k]);
      }
    }

    return volume_mesh;
  }

  // Update vertex coordinates from a volume mesh
  void _update_vertices(VolumeMesh &volume_mesh)
  {
    for (size_t i = 0; i < num_vertex_columns; i++)
    {
      const size_t start_index = this->vertices_offset[i];
      const size_t end_index = this->vertices_offset[i + 1];
      const size_t num_vertices = end_index - start_index;

      this->vertices[i].clear();
      this->vertices[i].reserve(num_vertices); // Reserve space to avoid reallocations

      for (size_t j = start_index; j < end_index; j++)
      {
        // Use reference to avoid unnecessary copies
        const Vector3D &v = volume_mesh.vertices[j];
        this->vertices[i].push_back(v);
      }
    }
  }
};

} // namespace DTCC_BUILDER

#endif