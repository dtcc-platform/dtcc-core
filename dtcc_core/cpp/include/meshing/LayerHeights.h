// include/meshing/LayerHeights.h
#pragma once
#ifndef DTCC_VOLUME_MESH_LAYER_HEIGHTS_H
#define DTCC_VOLUME_MESH_LAYER_HEIGHTS_H

#include <algorithm>
#include <cmath>
// #include <stack>
#include <vector>

#include "Geometry.h"
#include "meshing/BuilderMesh.h"

namespace DTCC_BUILDER
{
namespace VolumeMeshingUtilities
{ 
  // Compute layer heights for all faces in the ground mesh
  [[nodiscard]]
  inline std::vector<double> compute_layer_heights(BuilderMesh &mesh);

  namespace Utilities
  {
    // Compute ideal layer height for a regular tetrahedron
    inline double ideal_layer_height(double area)
    {
      return std::pow(2.0, 1.5) * std::pow(3.0, -0.75) * std::sqrt(area);
    }

    // Compute closest layer height index to a given height
    size_t closest_layer_height(double h, const std::vector<double> &layer_heights)
    {
      size_t min_index = 0;
      double min_diff = std::abs(std::log(h / layer_heights[0]));
      for (size_t i = 1; i < layer_heights.size(); i++)
      {
        const double diff = std::abs(std::log(h / layer_heights[i]));
        if (diff < min_diff)
        {
          min_index = i;
          min_diff = diff;
        }
      }
      return min_index;
    }

    // Assign face colors (closest layer height index)
    void assign_face_colors(const std::vector<double> &areas,
                            const std::vector<double> &layer_heights, BuilderMesh &mesh)

    {
      // Iterate over ground faces
      // mesh.face_colors.resize(areas.size());
      for (size_t i = 0; i < areas.size(); i++)
      {
        // Assign closest layer height index
        const double h = ideal_layer_height(areas[i]);
        const size_t face_color = closest_layer_height(h, layer_heights);
        mesh.face_colors[i] = face_color;
      }
    }

    // Check layer heights (deviation from ideal)
    void check_layer_heights(const std::vector<double> &areas,
                            const std::vector<double> &layer_heights, const BuilderMesh &mesh)
    {
      double max_error = 0.0;
      for (size_t i = 0; i < areas.size(); i++)
      {
        const double h = ideal_layer_height(areas[i]);
        const double H = layer_heights[mesh.face_colors[i]];
        const double e = std::abs(h - H) / h;
        max_error = std::max(max_error, e);
      }
      info("Max layer height error: " + str(100 * max_error, 2L) + "%");
    }
  } // namespace Utilities

  inline std::vector<double> compute_layer_heights(BuilderMesh &mesh){
    
 
    // Compute face areas
    const size_t num_faces = mesh.faces.size();
    std::vector<double> areas(num_faces);
    for (std::size_t i = 0; i < num_faces; i++)
    {
      const auto &v0 = mesh.vertices[mesh.faces[i].v0];
      const auto &v1 = mesh.vertices[mesh.faces[i].v1];
      const auto &v2 = mesh.vertices[mesh.faces[i].v2];
      areas[i] = Geometry::triangle_area(v0, v1, v2);
    }

    // Compute ideal min and max layer heights based on areas
    const double min_area = *std::min_element(areas.begin(), areas.end());
    const double max_area = *std::max_element(areas.begin(), areas.end());
    const double _min_height = Utilities::ideal_layer_height(min_area);
    const double _max_height = Utilities::ideal_layer_height(max_area);
    info("Ideal layer heights: [" + str(_min_height) + ", " + str(_max_height) + "]");

    // Compute optimal dyadic layer heights based on ideal min and max
    const double rho = _max_height / _min_height;
    const double mid = std::sqrt(_min_height * _max_height);
    const size_t steps = static_cast<int>(std::log2(rho) + 0.5);
    double min_height = mid / std::pow(2, steps / 2.0);

    std::vector<double> layer_heights;
    for (size_t i = 0; i < steps + 1; i++)
      layer_heights.push_back(min_height * std::pow(2.0, i));
    const double max_height = layer_heights.back();
    info("Adjusted layer heights: [" + str(min_height) + ", " + str(max_height) + "]");
    info("Number of layers: " + str(layer_heights.size()));

    // Assign face colors (closest layer height index)
    info("Assigning face colors...");
    Utilities::assign_face_colors(areas, layer_heights, mesh);
    Utilities::check_layer_heights(areas, layer_heights, mesh);

    // Build vertex and face mappings
    info("Building mapping from vertices to faces...");
    mesh.build_vertex_to_face_mapping();
    info("Building mapping from faces to faces...");
    mesh.build_face_to_face_mapping();

    // Iteratively reassign colors to avoid big jumps
    info("Reassigning colors to avoid big jumps...");
    size_t iteration = 0;
    const size_t max_color_iterations = 10;
    while (mesh.check_face_colors() > 0)
    {
      mesh.reassign_face_colors();
      if (++iteration == max_color_iterations)
      {
        error("Reached max color iterations");
      }
    }
    Utilities::check_layer_heights(areas, layer_heights, mesh);

    // Assign vertex colors and sort by color
    info("Assigning vertex colors...");
    mesh.assign_vertex_colors();
    mesh.sort_faces_by_vertex_color_and_index();

    // Assign face partitions
    info("Assigning face partitions...");
    mesh.assign_face_partitions();

    // Eliminate type 3 partitions
    info("Eliminating type 3 partitions...");
    mesh.eliminate_type_3_partitions();

    // Double-check face colors
    if (mesh.check_face_colors() > 0)
      error("Found big jumps after partition elimination");

    return layer_heights;
  }

} // namespace VolumeMeshingUtilities

} // namespace DTCC_BUILDER
#endif