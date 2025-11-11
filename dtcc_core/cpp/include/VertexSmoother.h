// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_VERTEX_SMOOTHER_H
#define DTCC_VERTEX_SMOOTHER_H

#include <unordered_set>

#include "Logging.h"
#include "Timer.h"
#include "model/Mesh.h"

namespace DTCC_BUILDER
{

class VertexSmoother
{
public:
  // Smooth mesh
  static void
  _smooth_mesh(Mesh &mesh, size_t num_smoothings, bool z_only = false)
  {
    info("Smoothing mesh...");
    Timer timer("smooth_mesh");

    // build vertex connectivity
    info("Building vertex connectivity");
    const size_t num_vertices = mesh.vertices.size();
    std::vector<std::unordered_set<size_t>> vertex_neighbors(num_vertices);
    auto vert_con_t = Timer("vertex connectivity");
    for (const auto &T : mesh.faces)
    {
      vertex_neighbors[T.v0].insert(T.v1);
      vertex_neighbors[T.v0].insert(T.v2);
      vertex_neighbors[T.v1].insert(T.v2);
      vertex_neighbors[T.v1].insert(T.v0);
      vertex_neighbors[T.v2].insert(T.v0);
      vertex_neighbors[T.v2].insert(T.v1);
    }
    vert_con_t.stop();

    // Smooth by setting each vertex coordinate to average of neighbors
    for (size_t n = 0; n < num_smoothings; n++)
    {
      info("Smoothing iteration " + str(n));
      for (size_t i = 0; i < num_vertices; i++)
      {
        Vector3D p{};
        for (const auto &j : vertex_neighbors[i])
          p += Vector3D(mesh.vertices[j]);
        p /= static_cast<float>(vertex_neighbors[i].size());
        if (z_only)
          mesh.vertices[i].z = p.z;
        else
          mesh.vertices[i] = p;
      }
    }
  }

  // Smooth 3D mesh
  static void smooth_mesh(VolumeMesh &volume_mesh, size_t num_smoothings)
  {
    info("Smoothing volume mesh...");
    Timer timer("smooth_mesh");

    // build vertex connectivity
    info("Building vertex connectivity");
    const size_t num_vertices = volume_mesh.vertices.size();
    std::vector<std::unordered_set<size_t>> vertex_neighbors(num_vertices);
    for (const auto &T : volume_mesh.cells)
    {
      vertex_neighbors[T.v0].insert(T.v1);
      vertex_neighbors[T.v0].insert(T.v2);
      vertex_neighbors[T.v0].insert(T.v3);
      vertex_neighbors[T.v1].insert(T.v2);
      vertex_neighbors[T.v1].insert(T.v3);
      vertex_neighbors[T.v1].insert(T.v0);
      vertex_neighbors[T.v2].insert(T.v3);
      vertex_neighbors[T.v2].insert(T.v0);
      vertex_neighbors[T.v2].insert(T.v1);
      vertex_neighbors[T.v3].insert(T.v0);
      vertex_neighbors[T.v3].insert(T.v1);
      vertex_neighbors[T.v3].insert(T.v2);
    }

    // Smooth by setting each vertex coordinate to average of neighbors
    for (size_t n = 0; n < num_smoothings; n++)
    {
      info("Smoothing iteration " + str(n));
      for (size_t i = 0; i < num_vertices; i++)
      {
        Vector3D p;
        for (const auto &j : vertex_neighbors[i])
          p += Vector3D(volume_mesh.vertices[j]);
        p /= static_cast<float>(vertex_neighbors[i].size());
        volume_mesh.vertices[i] = p;
      }
    }
  }

  static void smooth_mesh(Mesh &mesh, size_t num_smoothings, bool fix_building_vertices = false,
                            bool z_only = false)
  {
    info("Smoothing mesh...");
    Timer timer("smooth_mesh");

    const size_t num_vertices = mesh.vertices.size();
    const size_t num_faces = mesh.faces.size();
    if (num_faces == 0 || num_vertices == 0)
      return;

    auto vert_con_t = Timer("vertex connectivity");

    std::vector<std::unordered_set<size_t>> vertex_neighbors(num_vertices);
    // neighbors.shrink_to_fit();

    for (const auto &T : mesh.faces)
    {
      vertex_neighbors[T.v0].insert(T.v1);
      vertex_neighbors[T.v0].insert(T.v2);
      vertex_neighbors[T.v1].insert(T.v0);
      vertex_neighbors[T.v1].insert(T.v2);
      vertex_neighbors[T.v2].insert(T.v0);
      vertex_neighbors[T.v2].insert(T.v1);
    }
    vert_con_t.stop();

    std::vector<char> freeze_vertex;
    if (fix_building_vertices)
    {
      if (mesh.markers.size() != num_faces)
      {
        error("smooth_mesh: mesh.markers size does not match faces; cannot freeze "
              "building/platform vertices safely.");
      }
      else
      {
        freeze_vertex.assign(num_vertices, 0);
        for (size_t fi = 0; fi < num_faces; ++fi)
        {
          const int marker = mesh.markers[fi];
          if (marker >= 0)
          {
            const auto &f = mesh.faces[fi];
            freeze_vertex[f.v0] = 1;
            freeze_vertex[f.v1] = 1;
            freeze_vertex[f.v2] = 1;
          }
        }
      }
    }

    // std::vector<Vector3D> next(mesh.vertices.begin(), mesh.vertices.end());

    for (size_t it = 0; it < num_smoothings; ++it)
    {
      info("Smoothing iteration " + str(it));
      for (size_t i = 0; i < num_vertices; ++i)
      {
        if (fix_building_vertices && !freeze_vertex.empty() && freeze_vertex[i])
        {
          // next[i] = mesh.vertices[i]; // keep exactly
          continue;
        }

        const auto &N = vertex_neighbors[i];
        if (N.empty())
        {
          // next[i] = mesh.vertices[i]; // degree-0: keep
          continue;
        }

        Vector3D avg{0, 0, 0};
        for (size_t j : N)
          avg += Vector3D(mesh.vertices[j]);

        const double invN = 1.0 / static_cast<double>(N.size());
        avg *= invN;

        if (z_only)
        {
          // next[i] = mesh.vertices[i];
          mesh.vertices[i].z = avg.z;
        }
        else
        {
          mesh.vertices[i] = avg;
        }
      }
      // mesh.vertices.swap(next); // commit iteration
    }
  }

  // Smooth grid field
  static GridField smooth_field(const GridField &field, size_t num_smoothings)
  {
    info("Smoothing grid field...");
    Timer timer("smooth_field");

    // Create copy of field
    GridField _field{field};

    // Neighbor indices
    std::vector<size_t> indices{};
    indices.reserve(4);

    // Smooth by setting each value to average of neighbors
    for (size_t n = 0; n < num_smoothings; n++)
    {
      info("Smoothing iteration " + str(n));
      for (size_t i = 0; i < _field.values.size(); i++)
      {
        // Get neighbors
        indices.clear();
        _field.grid.index_to_boundary(i, indices);

        // Compute average
        double value = 0.0;
        for (const size_t &j : indices)
          value += _field.values[j];
        value /= static_cast<double>(indices.size());
        _field.values[i] = value;
      }
    }

    return _field;
  }
};

} // namespace DTCC_BUILDER

#endif
