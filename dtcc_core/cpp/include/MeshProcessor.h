// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_MESH_PROCESSOR_H
#define DTCC_MESH_PROCESSOR_H

#include <limits>
#include <set>
#include <unordered_map>
#include <utility>

#include "Hashing.h"
#include "Logging.h"

#include "KDTreeVectorOfVectorsAdaptor.h"
#include "nanoflann.hpp"

#include "DisjointSet.h"
#include "BoundingBox.h"
#include "BoundingBoxTree.h"

#include "model/Mesh.h"
#include "model/Simplices.h"
#include "model/Polygon.h"
#include "model/Vector.h"

namespace DTCC_BUILDER
{

class MeshProcessor
{
public:


static inline void compute_mesh_domain_markers(Mesh &mesh, const std::vector<Polygon> &subdomains)
{
  info("Computing domain markers...");
  Timer timer("compute_domain_markers");

  BoundingBoxTree2D search_tree;
  std::vector<BoundingBox2D> bounding_boxes;
  bounding_boxes.reserve(subdomains.size());
  for (const auto &subdomain : subdomains)
  {
    bounding_boxes.emplace_back(subdomain);
  }
  search_tree.build(bounding_boxes);

  mesh.markers.resize(mesh.faces.size());
  std::fill(mesh.markers.begin(), mesh.markers.end(), -2);

  std::vector<bool> is_building_vertex(mesh.vertices.size());
  std::fill(is_building_vertex.begin(), is_building_vertex.end(), false);

  if (!subdomains.empty())
  {
    for (size_t i = 0; i < mesh.faces.size(); i++)
    {
      const Vector3D c_3d = mesh.mid_point(i);
      const Vector2D c_2d(c_3d.x, c_3d.y);
      std::vector<size_t> indices = search_tree.find(Vector2D(c_2d));

      if (!indices.empty())
      {
        for (const auto &index : indices)
        {
          if (Geometry::polygon_contains_2d(subdomains[index], c_2d))
          {
            mesh.markers[i] = index;
            const Simplex2D &T = mesh.faces[i];
            is_building_vertex[T.v0] = true;
            is_building_vertex[T.v1] = true;
            is_building_vertex[T.v2] = true;
            break;
          }
        }
      }
    }

    for (size_t i = 0; i < mesh.faces.size(); i++)
    {
      const Simplex2D &T = mesh.faces[i];
      const bool touches_building =
          (is_building_vertex[T.v0] || is_building_vertex[T.v1] || is_building_vertex[T.v2]);

      if (touches_building && mesh.markers[i] == -2)
        mesh.markers[i] = -1;
    }
  }
}
  /// Merge meshes into a single mesh
  static Mesh compact_mesh(const Mesh &mesh)
  {
    return compact_mesh_vertices(mesh);
  }

  /// Merge meshes into a single mesh
  static Mesh merge_meshes(const std::vector<Mesh> &meshes, bool weld = false, double snap = 0)
  {
    // info("Merging " + str(meshes.size()) + " meshes into a single mesh...");
    Timer timer("merge_meshes");

    // Create empty mesh
    Mesh mesh;

    // Count the number of vertices and cells
    size_t num_vertices = 0;
    size_t num_cells = 0;
    for (const auto &m : meshes)
    {
      num_vertices += m.vertices.size();
      num_cells += m.faces.size();
    }

    // Allocate arrays
    mesh.vertices.resize(num_vertices);
    mesh.faces.resize(num_cells);
    mesh.markers.resize(num_cells, default_marker());
    // Merge data
    size_t vertex_offset = 0;
    size_t face_offset = 0;
    for (const auto &src_mesh : meshes)
    {
      const size_t current_vertex_offset = vertex_offset;

      for (size_t j = 0; j < src_mesh.faces.size(); ++j)
      {
        Simplex2D c = src_mesh.faces[j];
        c.v0 += current_vertex_offset;
        c.v1 += current_vertex_offset;
        c.v2 += current_vertex_offset;
        mesh.faces[face_offset] = c;

        int marker = default_marker();
        if (j < src_mesh.markers.size())
          marker = src_mesh.markers[j];
        mesh.markers[face_offset] = marker;

        ++face_offset;
      }

      for (size_t j = 0; j < src_mesh.vertices.size(); ++j)
        mesh.vertices[current_vertex_offset + j] = src_mesh.vertices[j];

      vertex_offset += src_mesh.vertices.size();
    }
    if (snap > 0)
      mesh = snap_vertices(mesh, snap);
    if (weld)
      mesh = weld_mesh(mesh);
    return mesh;
  }

  // return a list of naked edges and the faces that contain them
  static std::vector<std::pair<Simplex1D, Simplex2D>>
  find_naked_edges(const std::vector<Simplex2D> &faces)
  {
    std::unordered_map<Simplex1D, int, Simplex1DHash> edge_counts;
    std::unordered_map<Simplex1D, size_t, Simplex1DHash> edge_face;
    // Count the number of times each edge appears in the faces
    for (const auto &face : faces)
    {
      for (int i = 0; i < 3; ++i)
      {
        // info("edge: " + str(face[i]) + " " + str(face[(i + 1) % 3]));
        Simplex1D edge(face[i], face[(i + 1) % 3], true);
        if (edge_counts.find(edge) == edge_counts.end())
        {
          edge_counts.insert({edge, 1});
          edge_face.insert({edge, face[(i + 2) % 3]});
        }
        else
        {
          edge_counts[edge] += 1;
        }
      }
    }

    // Find the edges that appear only once
    std::vector<std::pair<Simplex1D, Simplex2D>> naked_edges;
    for (auto it = edge_counts.begin(); it != edge_counts.end(); ++it)
    {
      if (it->second == 1)
      {
        auto e = it->first;

        auto edge_pair = std::make_pair(e, Simplex2D(e.v0, e.v1, edge_face[e]));
        naked_edges.push_back(edge_pair);
      }
    }

    return naked_edges;
  }

  static Mesh weld_mesh(const Mesh &mesh)
  {
    // size_t num_vertices = mesh.vertices.size();
    Timer timer("weld_mesh");

    // Create empty mesh
    Mesh welded_mesh;
    std::unordered_map<Vector3D, size_t, Vector3DHash> vertex_map;
    std::unordered_map<size_t, size_t> vertex_idx_map;
    for (size_t i = 0; i < mesh.vertices.size(); ++i)
    {
      const Vector3D &vertex = mesh.vertices[i];
      if (vertex_map.find(vertex) == vertex_map.end())
      {
        vertex_map.insert({vertex, welded_mesh.vertices.size()});
        vertex_idx_map.insert({i, vertex_map.size() - 1});
        welded_mesh.vertices.push_back(vertex);
      }
      else
      {
        vertex_idx_map.insert({i, vertex_map[vertex]});
      }
    }
    welded_mesh.markers.reserve(mesh.faces.size());
    for (size_t i = 0; i < mesh.faces.size(); ++i)
    {
      const Simplex2D &face = mesh.faces[i];
      welded_mesh.faces.push_back(Simplex2D(vertex_idx_map[face.v0], vertex_idx_map[face.v1],
                                            vertex_idx_map[face.v2]));

      int marker = default_marker();
      if (i < mesh.markers.size())
        marker = mesh.markers[i];
      welded_mesh.markers.push_back(marker);
    }
    for (size_t i = 0; i < mesh.normals.size(); ++i)
    {
      welded_mesh.normals.push_back(mesh.normals[i]);
    }
    // info("welded " + str(num_vertices) + " vertices to " +
    //     str(welded_mesh.vertices.size()) + " vertices");
    return compact_mesh(welded_mesh);
  }

  static Mesh snap_vertices(const Mesh &mesh, double snap_distance)
  {
    Mesh snapped_mesh = mesh;
    // snap_distance *= snap_distance;

    auto merge_candidates = DisjointSet(mesh.vertices.size());

    typedef KDTreeVectorOfVectorsAdaptor<std::vector<Vector3D>, double, 3 /* dims */> my_kd_tree_t;
    my_kd_tree_t vert_index(3 /*dim*/, mesh.vertices, 10 /* max leaf */);
    vert_index.index->buildIndex();
    for (size_t i = 0; i < mesh.vertices.size(); ++i)
    {
      auto &pt = mesh.vertices[i];
      std::vector<double> query_pt{pt.x, pt.y, pt.z};
      auto neighbours = vert_index.radius_query(&query_pt[0], snap_distance);
      if (neighbours.size() > 1)
      {
        auto first_idx = neighbours[0].first;
        for (size_t j = 1; j < neighbours.size(); ++j)
        {
          size_t idx = neighbours[j].first;
          merge_candidates.unionSets(first_idx, idx);
        }
      }
    }

    auto merge_sets = merge_candidates.getSets();
//    size_t num_merged = 0;
    for (const auto &ms : merge_sets)
    {
      auto merge_group = ms.second;
      if (merge_group.size() > 1)
      {
        size_t target = ms.first;
        for (auto &face : snapped_mesh.faces)
        {

          auto face_normal = Geometry::face_normal(face, snapped_mesh);
          bool snapped = false;

          if ( (face.v0 != target) && (std::find(merge_group.begin(), merge_group.end(), face.v0) != merge_group.end()))
          {
            face.v0 = target;
            snapped = true;
          }
          if ((face.v1 != target) && (std::find(merge_group.begin(), merge_group.end(), face.v1) != merge_group.end()))
          {
            face.v1 = target;
            snapped = true;
          }
          if ((face.v2 != target) && (std::find(merge_group.begin(), merge_group.end(), face.v2) != merge_group.end()) )
          {
            face.v2 = target;
            snapped = true;
          }
          if (snapped)
          {
            // check and fix normals
            auto snapped_face_normal = Geometry::face_normal(face, snapped_mesh);
            if (face_normal.dot(snapped_face_normal) < 0)
            {
              std::swap(face.v1, face.v2);
            }

          }
        }
      }
    }
//    info("Merging " + str(num_merged) + " vertices");
    if (snapped_mesh.markers.size() < snapped_mesh.faces.size())
      snapped_mesh.markers.resize(snapped_mesh.faces.size(), default_marker());
    else if (snapped_mesh.markers.size() > snapped_mesh.faces.size())
      snapped_mesh.markers.resize(snapped_mesh.faces.size());

    return compact_mesh(snapped_mesh);
  }

private:
  static Mesh compact_mesh_vertices(const Mesh &mesh)
  {
    if (mesh.faces.empty())
      return mesh;

    Mesh compacted_mesh;
    compacted_mesh.vertices.reserve(mesh.faces.size() * 3);
    compacted_mesh.faces.reserve(mesh.faces.size());
    compacted_mesh.markers.reserve(mesh.faces.size());
    compacted_mesh.normals = mesh.normals;

    std::unordered_map<size_t, size_t> vertex_idx_map;
    for (size_t i = 0; i < mesh.faces.size(); ++i)
    {
      const Simplex2D &face = mesh.faces[i];

      auto map_vertex = [&](size_t vertex_index) -> size_t {
        auto it = vertex_idx_map.find(vertex_index);
        if (it != vertex_idx_map.end())
          return it->second;

        const size_t new_index = compacted_mesh.vertices.size();
        vertex_idx_map[vertex_index] = new_index;
        compacted_mesh.vertices.push_back(mesh.vertices[vertex_index]);
        return new_index;
      };

      // Call map_vertex in a fixed order: argument evaluation order is
      // unspecified in C++ (gcc and clang differ), and map_vertex assigns
      // new vertex indices by first use.
      const size_t new_v0 = map_vertex(face.v0);
      const size_t new_v1 = map_vertex(face.v1);
      const size_t new_v2 = map_vertex(face.v2);
      compacted_mesh.faces.push_back(Simplex2D(new_v0, new_v1, new_v2));

      int marker = default_marker();
      if (i < mesh.markers.size())
        marker = mesh.markers[i];
      compacted_mesh.markers.push_back(marker);
    }

    return compacted_mesh;
  }

  static int default_marker()
  {
    return std::numeric_limits<int>::lowest();
  }
};

} // namespace DTCC_BUILDER

#endif
