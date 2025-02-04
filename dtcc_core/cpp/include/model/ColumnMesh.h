#ifndef DTCC_COLUMN_MESH_H
#define DTCC_COLUMN_MESH_H

#include <map>

#include "VertexSmoother.h"
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

  // Layer index (integer height in terms of smallest layer height)
  size_t layer_index{};

  // Create default simplex
  ColumnSimplex() = default;

  // Create simplex and optionally sort vertices
  ColumnSimplex(ColumnIndex v0, ColumnIndex v1, ColumnIndex v2, ColumnIndex v3, size_t layer_index)
  {
    this->v0 = v0;
    this->v1 = v1;
    this->v2 = v2;
    this->v3 = v3;
    this->layer_index = layer_index;
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

  // Vector of number of prisms in each cell column
  std::vector<size_t> num_prisms;

  // Vector of vertex offsets
  std::vector<size_t> vertices_offset;

  // Number of layers of minimum height
  size_t num_min_layers{};

  // Number of layers of maximum height
  size_t num_max_layers{};

  // Default constructor
  ColumnMesh() = default;

  // Constructor
  ColumnMesh(const Mesh &ground_mesh)
  {
    vertices.resize(ground_mesh.vertices.size());
    cells.resize(ground_mesh.faces.size());
    markers.resize(ground_mesh.markers.size());
    num_prisms.resize(ground_mesh.faces.size(), 0);
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

  // Get cell height
  double cell_height(const ColumnSimplex &cell) const
  {
    const Vector3D &v0 = vertex(cell.v0);
    const Vector3D &v1 = vertex(cell.v1);
    const Vector3D &v2 = vertex(cell.v2);
    const Vector3D &v3 = vertex(cell.v3);

    return std::max({v0.z, v1.z, v2.z, v3.z});
  }

  // Get layer index of vertex
  size_t layer_index(const ColumnIndex &index) const
  {
    const size_t layer_step = num_min_layers / (vertices[index.column].size() - 1);
    return index.index * layer_step;
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

    // Add cells and renumber vertices
    std::unordered_map<size_t, size_t> vertex_indices;
    std::unordered_map<size_t, ColumnIndex> column_indices;
    for (size_t i = 0; i < cells.size(); i++)
    {
      for (size_t j = 0; j < cells[i].size(); j++)
      {
        // Skip if cell should be trimmed
        if (!keep_cells[i][j])
          continue;

        // Get vertex indices
        size_t i0 = vertices_offset[cells[i][j].v0.column] + cells[i][j].v0.index;
        size_t i1 = vertices_offset[cells[i][j].v1.column] + cells[i][j].v1.index;
        size_t i2 = vertices_offset[cells[i][j].v2.column] + cells[i][j].v2.index;
        size_t i3 = vertices_offset[cells[i][j].v3.column] + cells[i][j].v3.index;

        // Renumber vertices
        i0 = renumber_vertex(i0, cells[i][j].v0, vertex_indices, column_indices);
        i1 = renumber_vertex(i1, cells[i][j].v1, vertex_indices, column_indices);
        i2 = renumber_vertex(i2, cells[i][j].v2, vertex_indices, column_indices);
        i3 = renumber_vertex(i3, cells[i][j].v3, vertex_indices, column_indices);

        // Add cell
        volume_mesh.cells.push_back(Simplex3D(i0, i1, i2, i3));
      }
    }

    // Add vertices
    volume_mesh.vertices.resize(vertex_indices.size());
    for (const auto &it : vertex_indices)
    {
      const auto column_index = column_indices[it.first];
      volume_mesh.vertices[it.second] = vertices[column_index.column][column_index.index];
    }

    // Add markers
    volume_mesh.markers.resize(vertex_indices.size());
    for (const auto &it : vertex_indices)
    {
      const auto column_index = column_indices[it.first];
      volume_mesh.markers[it.second] = markers[column_index.column][column_index.index];
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

private:
  // Renumber vertex (if not already renumbered)
  size_t renumber_vertex(size_t vertex_index, ColumnIndex column_index,
                         std::unordered_map<size_t, size_t> &vertex_indices,
                         std::unordered_map<size_t, ColumnIndex> &column_indices)
  {
    // Check if vertex has already been renumbered
    const auto it = vertex_indices.find(vertex_index);
    if (it != vertex_indices.end())
      return it->second;

    // Renumber vertex
    const size_t new_index = vertex_indices.size();
    vertex_indices[vertex_index] = new_index;
    column_indices[vertex_index] = column_index;
    return new_index;
  }
};

} // namespace DTCC_BUILDER

#endif
