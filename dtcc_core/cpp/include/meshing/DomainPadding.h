#ifndef DTCC_VOLUME_MESH_DOMAIN_PADDING_H
#define DTCC_VOLUME_MESH_DOMAIN_PADDING_H

#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <utility>
#include <vector>

#include "../Logging.h"
#include "../model/Mesh.h"
#include "../model/VolumeMesh.h"

namespace DTCC_BUILDER
{
namespace VolumeMeshing
{
namespace DomainPaddingUtilities
{

static inline Mesh extract_top_mesh(const VolumeMesh &volume_mesh, std::vector<size_t> &new_to_old_index)
{
  info("Extracting surface from the top boundary of the domain");
  size_t num_vertices = 0;
  // std::vector<size_t> new_to_old_index;
  std::unordered_map<size_t, size_t> old_to_new_index;
  std::vector<Vector3D> top_vertices;
  std::vector<Simplex2D> top_faces;

  for (size_t i = 0; i < volume_mesh.vertices.size(); i++)
  {
    if (volume_mesh.markers[i] == -3)
    {
      top_vertices.push_back(volume_mesh.vertices[i]);
      old_to_new_index[i] = num_vertices;
      new_to_old_index.push_back(i);
      ++num_vertices;
    }
  }
  info("Top mesh vertices: " + str(top_vertices.size()));
  std::vector<int> top_markers(top_vertices.size(), -3);

  for (size_t i = 0; i < volume_mesh.cells.size(); i++)
  {
    const Simplex3D cell = volume_mesh.cells[i];

    const std::array<size_t, 4> cell_vertices = {cell.v0, cell.v1, cell.v2, cell.v3};

    std::vector<size_t> face_vertices_top;
    for (const auto v : cell_vertices)
    {
      if (volume_mesh.markers[v] == -3)
        face_vertices_top.push_back(old_to_new_index[v]);
    }

    if (face_vertices_top.size() == 3)
    {
      top_faces.push_back(
          Simplex2D(face_vertices_top[0], face_vertices_top[1], face_vertices_top[2]));
    }
  }
  info("Top mesh faces: " + str(top_faces.size()));
  Mesh top_mesh;
  top_mesh.vertices = top_vertices;
  top_mesh.faces = top_faces;

  info(top_mesh);

  return top_mesh;
}

static inline VolumeMesh weld_meshes(const VolumeMesh &volume_mesh, const VolumeMesh &padding_mesh,
                       std::vector<int> &vertex_matches)
{
  info("Merging Volume Mesh with Padding Mesh");
  // Merge the two meshes
  VolumeMesh merged_mesh;
  const size_t volume_mesh_num_vertices = volume_mesh.vertices.size();
  merged_mesh.vertices.insert(merged_mesh.vertices.end(), volume_mesh.vertices.begin(),
                              volume_mesh.vertices.end());
  merged_mesh.cells.insert(merged_mesh.cells.end(), volume_mesh.cells.begin(),
                           volume_mesh.cells.end());
  merged_mesh.markers.insert(merged_mesh.markers.end(), volume_mesh.markers.begin(),
                             volume_mesh.markers.end());

  // Replace old top vertex markers
  for (auto &m : merged_mesh.markers)
  {
    if (m == -3)
      m = -5;
  }

  size_t vertices_added = 0;
  for (size_t i = 0; i < padding_mesh.vertices.size(); i++)
  {
    if (vertex_matches[i] < 0)
    {
      merged_mesh.vertices.push_back(padding_mesh.vertices[i]);
      merged_mesh.markers.push_back(padding_mesh.markers[i]);
      vertex_matches[i] = volume_mesh_num_vertices + vertices_added;
      ++vertices_added;
    }
  }

  for (size_t i = 0; i < padding_mesh.cells.size(); i++)
  {
    const std::array<size_t, 4> cell = {padding_mesh.cells[i].v0, padding_mesh.cells[i].v1,
                                        padding_mesh.cells[i].v2, padding_mesh.cells[i].v3};
    std::array<size_t, 4> new_cell = {0, 0, 0, 0};
    for (size_t j = 0; j < 4; j++)
    {
      new_cell[j] = vertex_matches[cell[j]];
    }
    merged_mesh.cells.push_back(Simplex3D(new_cell[0], new_cell[1], new_cell[2], new_cell[3]));
  }

  return merged_mesh;
}

} // namespace DomainPaddingUtilities

} // namespace VolumeMeshing

} // namespace DTCC_BUILDER

#endif