#pragma once
#ifndef DTCC_VOLUME_MESH_BUILDER_H
#define DTCC_VOLUME_MESH_BUILDER_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <stack>
#include <tuple>
#include <vector>

#include "Logging.h"
#include "Smoother.h"
#include "Timer.h"

#include "model/ColumnMesh.h"
#include "model/GridField.h"
#include "model/Mesh.h"
#include "model/Surface.h"
#include "model/Vector.h"

#include "meshing/DomainPadding.h"
#include "meshing/BuilderMesh.h"
#include "meshing/LayerHeights.h"
#include "meshing/MeshImprovement.h"
#include "meshing/ColumnMeshProcessing.h"
#include "meshing/MeshLogging.h"

#include "MeshQualityMetrics.h"

namespace DTCC_BUILDER
{

class VolumeMeshBuilder
{
public:
  double domain_height;
  double top_height;

private:
  // FIXME: Use static functions instead, avoid all the many private
  // variables (cleaner)

  ColumnMesh _column_mesh;

  const std::vector<Surface> _buildings;

  std::vector<double> _building_ground_height;

  std::vector<double> adj_building_height;

  const GridField &_dem;

  BuilderMesh _ground_mesh;

  BoundingBox2D mesh_bounds;

  Vector3D mesh_center;

  Vector3D mesh_origin;

public:
  // Constructor
  VolumeMeshBuilder(const std::vector<Surface> &buildings, const GridField &dem, Mesh &ground_mesh,
                    double domain_height)
      : _buildings(buildings), _dem(dem), _ground_mesh(ground_mesh), domain_height(domain_height),
        _column_mesh(ground_mesh)
  {

    assert((!_ground_mesh.vertices.empty()) && "Ground mesh has no vertices");
    assert((!_ground_mesh.faces.empty()) && "Ground mesh has no faces");
    assert((!_ground_mesh.markers.empty()) && "Ground mesh has no markers");

    mesh_bounds = BoundingBox2D(ground_mesh.vertices);
    mesh_center = Vector3D((mesh_bounds.P.x + mesh_bounds.Q.x) / 2.0,
                           (mesh_bounds.P.y + mesh_bounds.Q.y) / 2.0,0.0);
    mesh_origin = Vector3D(mesh_bounds.P.x, mesh_bounds.P.y, 0.0); 
    

    compute_building_ground_heights();
  }

  // Destructor
  ~VolumeMeshBuilder() {}

  VolumeMesh build(const size_t smoother_iterations = 1000,
                   const double smoother_relative_tolerance = 0.001,
                   const double domain_padding_height = 0.0,
                   const double aspect_ratio_threshold = 10.0, const size_t debug_step = 6)
  {
    info("Building volume mesh...");

    MeshLogger logger;
    double min_ar, max_ar, median_ar = 0.0;
    // Layer ground mesh
    logger.step_begin("1: Ground mesh layering");
    VolumeMesh volume_mesh = layer_ground_mesh(_ground_mesh);
    logger.step_stop();
    std::tie(min_ar, max_ar, median_ar) = Geometry::aspect_ratio(volume_mesh);
    logger.log_step(min_ar,median_ar,max_ar);

    // Debugging
    if (debug_step == 2)
      return volume_mesh;

    // Volume mesh smoothing
    logger.step_begin("2: Volume mesh smoothing (ground only)");
    const bool fix_top = false;
    // volume_mesh = Smoother::smooth_volume_mesh_poisson(volume_mesh, _buildings, _dem, 0.0, false,
    //                                                    fix_top, smoother_iterations,
    //                                                    smoother_relative_tolerance);
    volume_mesh = Smoother::smooth_volume_mesh_elastic(volume_mesh, _buildings, _dem, 0.0, false,
                                                       fix_top, smoother_iterations,
                                                       smoother_relative_tolerance, mesh_bounds);
    _column_mesh._update_vertices(volume_mesh);
    logger.step_stop();
    std::tie(min_ar, max_ar, median_ar) = Geometry::aspect_ratio(volume_mesh);
    logger.log_step(min_ar,median_ar,max_ar);
    // check_mesh_quality(volume_mesh, 3);

    // Debugging
    if (debug_step == 3)
      return volume_mesh;

    // Trim volume mesh
    logger.step_begin("3: Volume mesh trimming");
    info("Trimming volume mesh...");
    volume_mesh = trim_volume_mesh();
    logger.step_stop();
    std::tie(min_ar, max_ar, median_ar) = Geometry::aspect_ratio(volume_mesh);
    logger.log_step(min_ar,median_ar,max_ar);
    

    // Debugging
    if (debug_step == 4)
      return volume_mesh;

    // FIXME: Smooth mesh in-place instead of returning a new mesh

    // Smooth volume mesh (again)
    logger.step_begin("4: Volume mesh smoothing (ground and buildings)");
    volume_mesh = Smoother::smooth_volume_mesh_elastic(volume_mesh, _buildings, _dem, 0.0, true,
                                                       false, smoother_iterations,
                                                       smoother_relative_tolerance, mesh_bounds);
    logger.step_stop();
    std::tie(min_ar, max_ar, median_ar) = Geometry::aspect_ratio(volume_mesh);
    logger.log_step(min_ar,median_ar,max_ar);

    if (debug_step == 5)
      return volume_mesh;

    logger.step_begin("5: Volume mesh improvement (Removing slivers)");
    volume_mesh = VolumeMeshing::MeshImprovement::remove_tetrahedra(volume_mesh, aspect_ratio_threshold);
    logger.step_stop();
    std::tie(min_ar, max_ar, median_ar) = Geometry::aspect_ratio(volume_mesh);
    logger.log_step(min_ar,median_ar,max_ar);

    if (debug_step == 6)
      return volume_mesh;

    logger.step_begin("6: Add Volume mesh Domain padding");
    top_height = domain_height + _dem.max();
    volume_mesh = layer_padding_mesh(volume_mesh, domain_height, 2.0);
    logger.step_stop();
    std::tie(min_ar, max_ar, median_ar) = Geometry::aspect_ratio(volume_mesh);
    logger.log_step(min_ar,median_ar,max_ar);

    logger.summary();

    TetrahedronMeshQuality::check_volume_mesh(volume_mesh);
    return volume_mesh;
  }

private:

  // Check mesh quality
  std::tuple<double,double,double> check_mesh_quality(const VolumeMesh &volume_mesh, int step, bool write_to_file = false)
  {
    // Compute aspect ratios
    const auto aspect_ratios = Geometry::aspect_ratio(volume_mesh);
    const double min = std::get<0>(aspect_ratios);
    const double max = std::get<1>(aspect_ratios);
    const double median = std::get<2>(aspect_ratios);

    // Write aspect ratios to file for debugging
    const auto _aspect_ratios = Geometry::aspect_ratios(volume_mesh);
    if (write_to_file)
    {
      std::ofstream file("aspect_ratios_" + str(step) + ".txt");
      for (const auto &ar : _aspect_ratios)
        file << ar << std::endl;
      file.close();
    }
    // Print aspect ratios
    info("Mesh quality (aspect ratio): min = " + str(min, 3L) + ", max = " + str(max, 3L) +
         ", median = " + str(median, 3L));

    return std::tie(min,max,median);
  }

void report(const VolumeMesh &volume_mesh, double elapsed_time, int step)
{
    // compute your mesh quality metrics as before
    double min_ar, max_ar, median_ar;
    std::tie(min_ar, max_ar, median_ar) = check_mesh_quality(volume_mesh, step);

    // single static stream object, reused across calls
    static std::ofstream ofs;

    if (step == 2) {
        // every time we hit step 2, start fresh
        if (ofs.is_open()) {
            ofs.close();
        }
        ofs.clear();  // clear any flags
        ofs.open("mesh_report.txt", std::ios::out | std::ios::trunc);
        if (!ofs) {
            std::cerr << "[MESH REPORT] ERROR: cannot open mesh_report.txt for truncate\n";
        }
    }
    else {
        // for steps > 2, make sure it's open in append mode
        if (!ofs.is_open()) {
            ofs.open("mesh_report.txt", std::ios::out | std::ios::app);
            if (!ofs) {
                std::cerr << "[MESH REPORT] ERROR: cannot open mesh_report.txt for append\n";
            }
        }
    }

    // build a single formatted line
    std::ostringstream line;
    line << "step="      << step
         << " vertices=" << volume_mesh.vertices.size()
         << " cells="    << volume_mesh.cells.size()
         << " time_s="   << std::fixed << std::setprecision(3) << elapsed_time
         << " min_ar="   << std::fixed << std::setprecision(3) << min_ar
         << " median_ar="<< std::fixed << std::setprecision(3) << median_ar
         << " max_ar="   << std::fixed << std::setprecision(3) << max_ar
         << "\n";

    // always echo to console
    std::cout << "[MESH REPORT] " << line.str();

    // write (and flush) to file
    if (ofs) {
        ofs << line.str();
        ofs.flush();
    }
}

  // Layer ground mesh
  VolumeMesh layer_ground_mesh(BuilderMesh &mesh)
  { 
    
    const std::vector<double> layer_heights = VolumeMeshingUtilities::compute_layer_heights(_ground_mesh);
    // Check if layer heights are empty
    if (layer_heights.empty())
      error("Empty layer heights vector");
    // Compute max building height

    double max_building_height = layer_heights.back();
    for (size_t i = 0; i < _buildings.size(); i++)
    {
      double building_height = _buildings[i].max_height() - _building_ground_height[i];

      // Validate that each building's height is positive and log an error with specific details if
      // it's not. This ensures data integrity and prevents unexpected behavior caused by invalid
      // building heights.
      if (building_height <= 0.0)
        error("Building " + str(i) + " height less or equal to 0.0m (" + str(building_height) +
              "m).");
      max_building_height = std::max(max_building_height, building_height);
    }

    // Compute number of layers of max height and min height
    _column_mesh.num_max_layers = std::ceil(max_building_height / layer_heights.back()) + 2;
    _column_mesh.num_min_layers = _column_mesh.num_max_layers << (layer_heights.size() - 1);

    info("Layers used (Min layer heigt): "+ str(_column_mesh.num_min_layers));
    info("Layers used (Max layer heigt): "+ str(_column_mesh.num_max_layers));
    // Layer vertices in columns
    for (size_t j = 0; j < mesh.vertices.size(); j++)
    {
      const Vector3D &vg = mesh.vertices[j];
      const double h = layer_heights[mesh.vertex_colors[j]];
      const size_t col_size = (_column_mesh.num_min_layers >> mesh.vertex_colors[j]) + 1;
      for (size_t i = 0; i < col_size; i++)
      {
        Vector3D v(vg.x, vg.y, i * h);
        _column_mesh.vertices[j].push_back(v);
      }
      _column_mesh.vertices_offset[j + 1] = _column_mesh.vertices_offset[j] + col_size;
    }

    // Set Column Mesh Cells
    VolumeMeshingUtilities::connect_column_mesh_cells(mesh, _column_mesh);

    // Set markers
    auto mesh_vertex_markers = mesh.face_to_vertex_markers();
    for (size_t i = 0; i < mesh.vertices.size(); i++)
    {
      auto &markers = _column_mesh.markers[i];
      // Choose default marker: if mesh_vertex_markers[i] is >= 0 then -4 (building wall), otherwise
      // -5 (Other)
      if (mesh_vertex_markers[i] >= 0)
      {
        // markers.resize(_column_mesh.vertices[i].size(), -4);
        markers.resize(_column_mesh.vertices[i].size(), mesh_vertex_markers[i]);
      }
      else
      {
        // Initialize all markers to Other (-5)
        markers.resize(_column_mesh.vertices[i].size(), -5);
        // Set first (bottom) marker to Ground (-2) if ground otherwise Halo (-1)
        markers.front() = mesh_vertex_markers[i] == -2 ? -2 : -1;
      }

      // Set last (top) marker to Top (-3)
      markers.back() = -3;
    }

    return _column_mesh.to_volume_mesh();
  }

  // Trim volume mesh
  VolumeMesh trim_volume_mesh()
  {
    // The volume mesh is trimmed by removing cells that are below the top of
    // each building. To find out where we can trim the cells in each column
    // we first find which layer indices are common to all cells in each
    // building. This tells us at which layer index height we can trim the
    // building (without cutting any prisms). The cut is made at the height
    // closest to the building height.

    // Initialize allowed building trimming layer indices and distances
    std::vector<std::vector<std::pair<bool, double>>> building_layer_indices(_buildings.size());
    for (size_t i = 0; i < _buildings.size(); i++)
    {
      building_layer_indices[i].resize(_column_mesh.num_min_layers,
                                       {true, std::numeric_limits<double>::max()});
    }

    // Iterate over cell columns
    for (size_t i = 0; i < _column_mesh.cells.size(); i++)
    {
      // Skip if not building
      const int marker = _ground_mesh.markers[i];
      if (marker < 0)
        continue;

      // Mark which layer indices are *not* allowed for this column
      const size_t cell_layer_step = _column_mesh.num_min_layers / _column_mesh.num_prisms[i];
      for (size_t j = 0; j < _column_mesh.num_min_layers; j++)
      {
        if ((j + 1) % cell_layer_step != 0)
          building_layer_indices[marker][j] = {false, 0.0};
      }

      // Get building height
      const double building_height = _buildings[marker].max_height();

      // Iterate over vertex columns for current cell column
      const auto &face = _ground_mesh.faces[i];
      for (size_t j : {face.v0, face.v1, face.v2})
      {
        // Compute vertex layer step
        const size_t vertex_layer_step =
            _column_mesh.num_min_layers / (_column_mesh.vertices[j].size() - 1);

        // Iterate over vertices in column
        for (size_t k = 1; k < _column_mesh.vertices[j].size(); k++)
        {
          // Skip if not allowed layer index
          const size_t layer_index = k * vertex_layer_step;
          if (!building_layer_indices[marker][layer_index - 1].first)
            continue;

          // Check for closest layer to building height
          const double z = _column_mesh.vertices[j][k].z;
          const double d = std::abs(z - building_height);

          if (d < building_layer_indices[marker][layer_index - 1].second)
            building_layer_indices[marker][layer_index - 1] = {true, d};
        }
      }
    }

    // Compute building trimming layer indices based on allowed layer indices
    std::vector<size_t> trimming_layer_indices(_buildings.size(), 0);
    for (size_t i = 0; i < building_layer_indices.size(); i++)
    {
      // Find minimal distance among allowed layer indices
      size_t min_index = 0;
      double min_distance = std::numeric_limits<double>::max();
      for (size_t j = 0; j < building_layer_indices[i].size(); j++)
      {
        if (building_layer_indices[i][j].first &&
            building_layer_indices[i][j].second < min_distance)
        {
          min_index = j;
          min_distance = building_layer_indices[i][j].second;
        }
      }

      // Sanity check
      if (min_distance == std::numeric_limits<double>::max())
        error("No allowed layer indices for building " + str(i));

      // Set trimming index
      trimming_layer_indices[i] = min_index + 1;
    }

    // Mark which cells should be trimmed (or rather kept)
    std::vector<std::vector<bool>> keep_cells(_column_mesh.cells.size());
    for (size_t i = 0; i < _column_mesh.cells.size(); i++)
    {
      // Keep all cells by default
      keep_cells[i].resize(_column_mesh.cells[i].size(), true);

      // Skip if not building
      const int marker = _ground_mesh.markers[i];
      if (marker < 0)
        continue;

      // Iterate over cells in column
      for (size_t j = 0; j < _column_mesh.cells[i].size(); j++)
      {
        // Trim cell if inside trimming layer
        const auto &cell = _column_mesh.cells[i][j];
        if (cell.layer_index <= trimming_layer_indices[marker])
          keep_cells[i][j] = false;

        // // Mark top vertices as building
        // Assuming cell.v0, cell.v1, cell.v2, cell.v3 have the same type
        for (const auto &vertex : {cell.v0, cell.v1, cell.v2, cell.v3})
        {
          if (_column_mesh.layer_index(vertex) == trimming_layer_indices[marker])
          {
            // Mark this vertex with the given marker.
            //_column_mesh.markers[vertex.column][vertex.index] = marker;
            _column_mesh.markers[vertex.column][vertex.index] = _buildings.size() + marker;

            // Comment in the following to change marker types above buildings.
            const size_t col_size = _column_mesh.vertices[vertex.column].size();
            for (size_t v = vertex.index + 1; v < col_size; v++)
            { 
              _column_mesh.markers[vertex.column][v] = -4;
              //_column_mesh.markers[vertex.column][v] = -5;
            }
            _column_mesh.markers[vertex.column].back() = -3;
          }
        }
      }
    }

    // auto mesh_vertex_markers = face_to_vertex_markers();
    // for (size_t i = 0; i < _ground_mesh.vertices.size(); i++)
    // {
    //   const int marker = mesh_vertex_markers[i];
    //   if (marker < -1)
    //     continue;

    //   // std::cout<< "Vertex Column: " << i << " Ground Marker: "<< marker << "\t| ";
    //   for (size_t j = 0; j < _column_mesh.vertices[i].size(); j++)
    //   {
    //     std::cout << _column_mesh.vertices[i][j].x <<", "<< _column_mesh.vertices[i][j].y << ","
    //     <<_column_mesh.vertices[i][j].z << "," << _column_mesh.markers[i][j] << std::endl;
    //   }
    // }

    return _column_mesh.to_volume_mesh(keep_cells);
  }

  /// Computes ground heights at the building centroids
 void compute_building_ground_heights()
  {

    const size_t num_buildings = _buildings.size();
    _building_ground_height.resize(num_buildings);

    for (size_t i = 0; i < _buildings.size(); i++)
    {

      const Vector3D centroid = Geometry::surface_centroid(_buildings[i]);

      _building_ground_height[i] = _dem(centroid);
    }
  }

  VolumeMesh layer_padding_mesh(const VolumeMesh &volume_mesh, double top_height,
                                double max_scale = 2.0)
  {
    // Get max layer height
    std::vector<size_t> new_to_old_index;
    Mesh _top_mesh = VolumeMeshing::DomainPaddingUtilities::extract_top_mesh(volume_mesh, new_to_old_index);
    BuilderMesh mesh(_top_mesh);
    double top_surface_max_z = 0.0;
    double top_surface_min_z = std::numeric_limits<double>::max();
    for (auto v : mesh.vertices)
    {
      top_surface_max_z = std::max(top_surface_max_z, v.z);
      top_surface_min_z = std::min(top_surface_min_z, v.z);
    }

    
    const auto layer_heights = VolumeMeshingUtilities::compute_layer_heights(mesh);
    const double max_layer_height = layer_heights.back();

    if (top_height <= top_surface_max_z){
      info("Domain is taller than input top height!");
      top_height = top_surface_max_z + max_layer_height;
    }
    double padding_height = top_height - top_surface_max_z;
    
    std::vector<int> vertex_matches;
    if (new_to_old_index.size() != mesh.vertices.size())
    {
      std::cout << "number of vertices between old and new mesh " << new_to_old_index.size()
                << " | " << mesh.vertices.size() << std::endl;
      error("Mismatch in number of vertices between old and new mesh");
    }

    ColumnMesh padding_mesh(_top_mesh);
    
    
    // Number of Max height layers needed to cover padding height.
    // const size_t n =  static_cast<size_t>(2 * padding_height / (max_layer_height * (max_scale + 1)));
    const size_t n = std::max(
        static_cast<size_t>(2u),
        static_cast<size_t>(2 * padding_height / (max_layer_height * (max_scale + 1)))
      );

    info("Top surface max z:" + str(top_surface_max_z) + " Top domain height " + str(top_height) + " Padding Height: " + str(padding_height));
    info("Padding Domain with " + str(n) + " max layers scaled from 1 to " + str(max_scale));
    // std::vector<size_t> offset_before_padding = _column_mesh.vertices_offset;

    for (size_t i = 0; i < mesh.vertices.size(); i++)
    {
      const double top_vertex_z = mesh.vertices[i].z;
      const size_t color = mesh.vertex_colors[i];

      // Number of sub-layers for this column
      const size_t n_k = static_cast<size_t>(max_layer_height / layer_heights[color]);
      const double s_1 = n > 1 ? (max_scale - 1) / static_cast<double>(n - 1) : 0.0;

      // Precompute the expected unscaled total increment (H) for this column.
      double H = n_k * layer_heights[color] * ((n - 1) + s_1 * ((n - 1) * n / 2.0));

      // Compute the normalization scale factor.
      double scale_norm_factor = (H > 0.0) ? ((top_height - top_vertex_z) / H) : 1.0;

      // Now, generate vertices using the scaled increments.
      double cumulative_height = 0.0;
      padding_mesh.vertices[i].push_back(mesh.vertices[i]);
      vertex_matches.push_back(new_to_old_index[i]);
      for (size_t j = 1; j < n; j++)
      {
        // const double s = 1.0 + j*(max_scale - 1)/(n-1);
        const double s = (1.0 + j * s_1) * scale_norm_factor;

        for (size_t k = 0; k < n_k; k++)
        {
          cumulative_height += s * layer_heights[color];
          Vector3D v(mesh.vertices[i].x, mesh.vertices[i].y, top_vertex_z + cumulative_height);
          padding_mesh.vertices[i].push_back(v);
          vertex_matches.push_back(-1);
        }
      }
      padding_mesh.vertices_offset[i + 1] =
          padding_mesh.vertices_offset[i] + padding_mesh.vertices[i].size();
    }

    VolumeMeshingUtilities::connect_column_mesh_cells(mesh, padding_mesh);

    padding_mesh.markers.resize(mesh.vertices.size());
    for (size_t i = 0; i < mesh.vertices.size(); i++)
    {
      auto &markers = padding_mesh.markers[i];

      // Initialize all markers to Other (-5)
      markers.resize(padding_mesh.vertices[i].size(), -5);
      // Set last (top) marker to Top (-3)
      markers.back() = -3;
    }
    auto _padding_mesh = padding_mesh.to_volume_mesh();
    info("Volume Mesh " + str(volume_mesh));
    info("Padding Mesh " + str(_padding_mesh));
    
    return VolumeMeshing::DomainPaddingUtilities::weld_meshes(volume_mesh, _padding_mesh, vertex_matches);
  }

};

} // namespace DTCC_BUILDER

#endif
