#pragma once
#ifndef DTCC_VOLUME_MESH_COLUMN_MESH_PROCESSING_H
#define DTCC_VOLUME_MESH_COLUMN_MESH_PROCESSING_H

#include "meshing/BuilderMesh.h"
#include "model/VolumeMesh.h"
#include "model/ColumnMesh.h"

namespace DTCC_BUILDER
{
namespace VolumeMeshingUtilities
{ 

inline void connect_column_mesh_cells(BuilderMesh &mesh, ColumnMesh &column_mesh)
  {

    for (size_t i = 0; i < mesh.faces.size(); i++)
    {
      // Get face vertices and colors
      const Simplex2D face_simplex = mesh.faces[i];
      const std::array<size_t, 3> face = {face_simplex.v0, face_simplex.v1, face_simplex.v2};
      const std::array<int, 3> _vertex_colors = {mesh.vertex_colors[face_simplex.v0],
                                                 mesh.vertex_colors[face_simplex.v1],
                                                 mesh.vertex_colors[face_simplex.v2]};

      // Calculate offsets and sizes
      const std::array<size_t, 3> col_offsets = {0, 0, 0};
      const std::array<size_t, 3> col_sizes = {
          column_mesh.vertices_offset[face[0] + 1] - column_mesh.vertices_offset[face[0]],
          column_mesh.vertices_offset[face[1] + 1] - column_mesh.vertices_offset[face[1]],
          column_mesh.vertices_offset[face[2] + 1] - column_mesh.vertices_offset[face[2]]};

      // Calculate number of prisms
      const size_t num_prisms =
          (col_sizes.back() - 1) / (1 << (mesh.face_colors[i] - _vertex_colors[2]));
      column_mesh.num_prisms[i] = num_prisms;

      // Create prism iterator
      std::vector<std::array<size_t, 4>> prism_iterator(num_prisms);
      for (size_t j = 0; j < num_prisms; j++)
      {
        const size_t layer_index = (j + 1) * column_mesh.num_min_layers / num_prisms;
        prism_iterator[j] = {col_offsets[0] + j * (1 << (mesh.face_colors[i] - _vertex_colors[0])),
                             col_offsets[1] + j * (1 << (mesh.face_colors[i] - _vertex_colors[1])),
                             col_offsets[2] + j * (1 << (mesh.face_colors[i] - _vertex_colors[2])),
                             layer_index};
      }

      // Add tetrahedrons based on face partition type
      switch (mesh.face_partitions[i])
      {
      case 0: // 6 vertices, 3 tetrahedrons
      {
        for (const auto &ar : prism_iterator)
        {
          const size_t k = ar[0];
          const size_t l = ar[1];
          const size_t m = ar[2];
          const size_t n = ar[3];

          ColumnIndex bot_triangle_0(face[0], k);
          ColumnIndex bot_triangle_1(face[1], l);
          ColumnIndex bot_triangle_2(face[2], m);
          ColumnIndex top_triangle_0(face[0], k + 1);
          ColumnIndex top_triangle_1(face[1], l + 1);
          ColumnIndex top_triangle_2(face[2], m + 1);

          ColumnSimplex K0(bot_triangle_0, bot_triangle_1, bot_triangle_2, top_triangle_2, n);
          ColumnSimplex K1(bot_triangle_0, top_triangle_1, bot_triangle_1, top_triangle_2, n);
          ColumnSimplex K2(bot_triangle_0, top_triangle_0, top_triangle_1, top_triangle_2, n);

          column_mesh.cells[i].emplace_back(K0);
          column_mesh.cells[i].emplace_back(K1);
          column_mesh.cells[i].emplace_back(K2);
        }
      }
      break;

      case 1: // 7 vertices, 4 tetrahedrons
      {
        for (const auto &ar : prism_iterator)
        {
          const size_t k = ar[0];
          const size_t l = ar[1];
          const size_t m = ar[2];
          const size_t n = ar[3];

          ColumnIndex bot_triangle_0(face[0], k);
          ColumnIndex bot_triangle_1(face[1], l);
          ColumnIndex bot_triangle_2(face[2], m);
          ColumnIndex mid_triangle_0(face[0], k + 1);
          ColumnIndex top_triangle_0(face[0], k + 2);
          ColumnIndex top_triangle_1(face[1], l + 1);
          ColumnIndex top_triangle_2(face[2], m + 1);

          ColumnSimplex K0(bot_triangle_0, bot_triangle_1, bot_triangle_2, mid_triangle_0, n);
          ColumnSimplex K1(bot_triangle_1, top_triangle_2, bot_triangle_2, mid_triangle_0, n);
          ColumnSimplex K2(bot_triangle_1, top_triangle_2, mid_triangle_0, top_triangle_1, n);
          ColumnSimplex K3(mid_triangle_0, top_triangle_0, top_triangle_1, top_triangle_2, n);

          column_mesh.cells[i].emplace_back(K0);
          column_mesh.cells[i].emplace_back(K1);
          column_mesh.cells[i].emplace_back(K2);
          column_mesh.cells[i].emplace_back(K3);
        }
      }
      break;

      case 2: // 8 vertices, 5 tetrahedrons
      {
        for (const auto &ar : prism_iterator)
        {
          const size_t k = ar[0];
          const size_t l = ar[1];
          const size_t m = ar[2];
          const size_t n = ar[3];

          ColumnIndex bot_triangle_0(face[0], k);
          ColumnIndex bot_triangle_1(face[1], l);
          ColumnIndex bot_triangle_2(face[2], m);
          ColumnIndex mid_triangle_0(face[0], k + 1);
          ColumnIndex mid_triangle_1(face[1], l + 1);
          ColumnIndex top_triangle_0(face[0], k + 2);
          ColumnIndex top_triangle_1(face[1], l + 2);
          ColumnIndex top_triangle_2(face[2], m + 1);

          ColumnSimplex K0(bot_triangle_0, bot_triangle_1, bot_triangle_2, mid_triangle_1, n);
          ColumnSimplex K1(bot_triangle_0, bot_triangle_2, mid_triangle_0, mid_triangle_1, n);
          ColumnSimplex K2(bot_triangle_2, top_triangle_2, mid_triangle_0, mid_triangle_1, n);
          ColumnSimplex K3(top_triangle_0, top_triangle_2, top_triangle_1, mid_triangle_0, n);
          ColumnSimplex K4(top_triangle_1, top_triangle_2, mid_triangle_1, mid_triangle_0, n);

          column_mesh.cells[i].emplace_back(K0);
          column_mesh.cells[i].emplace_back(K1);
          column_mesh.cells[i].emplace_back(K2);
          column_mesh.cells[i].emplace_back(K3);
          column_mesh.cells[i].emplace_back(K4);
        }
      }
      break;

      default:
        error("Unhandled partition type: " + str(mesh.face_partitions[i]));
        break;
      }
    }
  }

} // namespace VolumeMeshingUtilities

} // namespace DTCC_BUILDER

#endif