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
public:
  /// Vector of vectors of vertices
  std::vector<std::vector<Vector3D>> vertices{};

  /// Vector of vectors of cells (tetrahedra)
  std::vector<std::vector<ColumnSimplex>> cells{};

  /// Vector of vectors of vertex markers
  std::vector<std::vector<int>> markers{};

  // Vector of vertex offsets
  std::vector<size_t> vertices_offset;

  // Default constructor
  ColumnMesh() = default;

  // Constructor
  ColumnMesh(const Mesh &ground_mesh)
  {
    vertices.resize(ground_mesh.vertices.size());
    cells.resize(ground_mesh.faces.size());
    markers.resize(ground_mesh.markers.size());
    vertices_offset.resize(ground_mesh.vertices.size() + 1, 0);
  }

  // Destructor
  virtual ~ColumnMesh() {}

  // Get vertex by column index
  const Vector3D &vertex(const ColumnIndex &index) const
  {
    return vertices[index.column][index.index];
  }

  // Get cell centroid
  Vector3D cell_centroid(const ColumnSimplex &cell) const
  {
    const Vector3D &v0 = vertex(cell.v0);
    const Vector3D &v1 = vertex(cell.v1);
    const Vector3D &v2 = vertex(cell.v2);
    const Vector3D &v3 = vertex(cell.v3);

    return (v0 + v1 + v2 + v3) / 4.0;
  }

  // Convert to volume mesh
  VolumeMesh to_volume_mesh()
  {
    VolumeMesh volume_mesh;

    // Add vertices
    const size_t volume_mesh_num_vertices = vertices_offset.back() + vertices.back().size();
    volume_mesh.vertices.reserve(volume_mesh_num_vertices);
    for (size_t j = 0; j < vertices.size(); j++)
    {
      for (size_t k = 0; k < vertices[j].size(); k++)
      {
        volume_mesh.vertices.push_back(vertices[j][k]);
      }
    }

    // Add cells
    for (size_t i = 0; i < cells.size(); i++)
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

    // Add markers
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

  // Convert to volume mesh (including trimming)
  VolumeMesh to_volume_mesh(const std::vector<std::vector<bool>> &keep_cells)
  {
    VolumeMesh volume_mesh;

    // Add vertices
    const size_t volume_mesh_num_vertices = vertices_offset.back() + vertices.back().size();
    volume_mesh.vertices.reserve(volume_mesh_num_vertices);
    for (size_t j = 0; j < vertices.size(); j++)
    {
      for (size_t k = 0; k < vertices[j].size(); k++)
      {
        volume_mesh.vertices.push_back(vertices[j][k]);
      }
    }

    // Add cells
    for (size_t i = 0; i < cells.size(); i++)
    {
      for (size_t j = 0; j < cells[i].size(); j++)
      {
        // Skip if cell should be trimmed
        if (!keep_cells[i][j])
          continue;

        size_t vc0 = vertices_offset[cells[i][j].v0.column] + cells[i][j].v0.index;
        size_t vc1 = vertices_offset[cells[i][j].v1.column] + cells[i][j].v1.index;
        size_t vc2 = vertices_offset[cells[i][j].v2.column] + cells[i][j].v2.index;
        size_t vc3 = vertices_offset[cells[i][j].v3.column] + cells[i][j].v3.index;

        volume_mesh.cells.push_back(Simplex3D(vc0, vc1, vc2, vc3));
      }
    }

    // Add markers
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
    for (size_t i = 0; i < vertices.size(); i++)
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

  // Pretty-print
  std::string __str__() const override
  {
    return "ColumnMesh mesh with " + str(vertices.size()) + " vertex columns and " +
           str(cells.size()) + " cell columns";
  }
};

} // namespace DTCC_BUILDER

#endif
