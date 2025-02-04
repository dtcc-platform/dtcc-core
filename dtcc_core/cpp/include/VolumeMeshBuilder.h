#ifndef DTCC_VOLUME_MESH_BUILDER_H
#define DTCC_VOLUME_MESH_BUILDER_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <stack>
#include <tuple>
#include <vector>

#include "Geometry.h"
#include "Logging.h"
#include "MeshProcessor.h"
#include "Smoother.h"
#include "Timer.h"
#include "VertexSmoother.h"
#include "model/GridField.h"
#include "model/Mesh.h"
#include "model/Surface.h"
#include "model/Vector.h"

#include "model/ColumnMesh.h"

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

  Mesh &_ground_mesh;

  // Stores vertex to face mapping for ground mesh.
  std::vector<std::unordered_set<int>> vf;

  // Stores vertex to face mapping for ground mesh.
  std::vector<std::unordered_set<int>> ff;

  std::vector<int> vertex_markers;

  std::vector<int> face_colors;
  std::vector<int> vertex_colors;

  std::vector<size_t> vert_skip;

  std::vector<int> face_partitions;

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

    compute_building_ground_heights();
  }

  // Destructor
  ~VolumeMeshBuilder() {}

  VolumeMesh build(const size_t smoother_iterations = 1000,
                   const double smoother_relative_tolerance = 0.001,
                   const double domain_padding_height = 0.0, const size_t debug_step = 5)
  {
    info("Building volume mesh...");

    Timer t3_1("Step 3.1: Layer height computation");
    info("Computing layer heights...");
    const auto layer_heights = compute_layer_heights(_ground_mesh);
    t3_1.stop();
    t3_1.print();

    // Check if layer heights are empty
    if (layer_heights.empty())
      error("Empty layer heights vector");

    // Layer ground mesh
    Timer t3_2("Step 3.2: Ground mesh layering");
    info("Layering ground mesh...");
    VolumeMesh volume_mesh = layer_ground_mesh(layer_heights);
    t3_2.stop();
    t3_2.print();
    check_mesh_quality(volume_mesh, 2);

    // Debugging
    if (debug_step == 2)
      return volume_mesh;

    // Volume mesh smoothing
    Timer t3_3("Step 3.3: Volume mesh smoothing (ground only)");
    info("Smoothing volume mesh...");
    volume_mesh = Smoother::smooth_volume_mesh(volume_mesh, _buildings, _dem, 0.0, false, false,
                                               smoother_iterations, smoother_relative_tolerance);
    _column_mesh._update_vertices(volume_mesh);
    t3_3.stop();
    t3_3.print();
    check_mesh_quality(volume_mesh, 3);

    // Debugging
    if (debug_step == 3)
      return volume_mesh;

    // Trim volume mesh
    Timer t3_4("Step 3.4: Volume mesh trimming");
    info("Trimming volume mesh...");
    volume_mesh = trim_volume_mesh();
    t3_4.stop();
    t3_4.print();
    check_mesh_quality(volume_mesh, 4);

    // Debugging
    if (debug_step == 4)
      return volume_mesh;

    // FIXME: Smooth mesh in-place instead of returning a new mesh

    // Smooth volume mesh (again)
    Timer t3_5("Step 3.5: Volume mesh smoothing (ground and buildings)");
    info("Smoothing volume mesh...");
    volume_mesh = Smoother::smooth_volume_mesh(volume_mesh, _buildings, _dem, 0.0, true, false,
                                               smoother_iterations, smoother_relative_tolerance);
    t3_5.stop();
    t3_5.print();
    check_mesh_quality(volume_mesh, 5);

    // FIXME: Not used yet (only for padding)
    top_height = domain_height + _dem.max();

    return volume_mesh;
  }

private:
  // Compute layer heights for all faces in the ground mesh
  std::vector<double> compute_layer_heights(Mesh &ground_mesh)
  {
    // Compute face areas
    const size_t num_faces = ground_mesh.faces.size();
    std::vector<double> areas(num_faces);
    for (std::size_t i = 0; i < num_faces; i++)
    {
      const auto &v0 = ground_mesh.vertices[ground_mesh.faces[i].v0];
      const auto &v1 = ground_mesh.vertices[ground_mesh.faces[i].v1];
      const auto &v2 = ground_mesh.vertices[ground_mesh.faces[i].v2];
      areas[i] = Geometry::triangle_area(v0, v1, v2);
    }

    // Compute ideal min and max layer heights based on areas
    const double min_area = *std::min_element(areas.begin(), areas.end());
    const double max_area = *std::max_element(areas.begin(), areas.end());
    const double _min_height = ideal_layer_height(min_area);
    const double _max_height = ideal_layer_height(max_area);
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
    assign_face_colors(areas, layer_heights);
    check_layer_heights(areas, layer_heights);

    // Build vertex and face mappings
    info("Building mapping from vertices to faces...");
    build_vertex_to_face_mapping();
    info("Building mapping from faces to faces...");
    build_face_to_face_mapping();

    // Iteratively reassign colors to avoid big jumps
    info("Reassigning colors to avoid big jumps...");
    size_t iteration = 0;
    const size_t max_color_iterations = 10;
    while (check_face_colors() > 0)
    {
      reassign_face_colors();
      if (++iteration == max_color_iterations)
      {
        error("Reached max color iterations");
      }
    }
    check_layer_heights(areas, layer_heights);

    // Assign vertex colors and sort by color
    info("Assigning vertex colors...");
    assign_vertex_colors(layer_heights);
    sort_faces_by_vertex_color_and_index();

    // Assign face partitions
    info("Assigning face partitions...");
    assign_face_partitions();

    // Eliminate type 3 partitions
    info("Eliminating type 3 partitions...");
    eliminate_type_3_partitions();

    // Double-check face colors
    if (check_face_colors() > 0)
      error("Found big jumps after partition elimination");

    return layer_heights;
  }

  // Compute ideal layer height for a regular tetrahedron
  double ideal_layer_height(double area)
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

  // Check layer heights (deviation from ideal)
  void check_layer_heights(const std::vector<double> &areas,
                           const std::vector<double> &layer_heights)
  {
    double max_error = 0.0;
    for (size_t i = 0; i < areas.size(); i++)
    {
      const double h = ideal_layer_height(areas[i]);
      const double H = layer_heights[face_colors[i]];
      const double e = std::abs(h - H) / h;
      max_error = std::max(max_error, e);
    }
    info("Max layer height error: " + str(100 * max_error, 2L) + "%");
  }

  // Check mesh quality
  void check_mesh_quality(const VolumeMesh &volume_mesh, int step)
  {
    // Compute aspect ratios
    const auto aspect_ratios = Geometry::aspect_ratio(volume_mesh);
    const double min = std::get<0>(aspect_ratios);
    const double max = std::get<1>(aspect_ratios);
    const double median = std::get<2>(aspect_ratios);

    // Write aspect ratios to file for debugging
    const auto _aspect_rations = Geometry::aspect_ratios(volume_mesh);
    std::ofstream file("aspect_ratios_" + str(step) + ".txt");
    for (const auto &ar : _aspect_rations)
      file << ar << std::endl;
    file.close();

    // Print aspect ratios
    info("Mesh quality (aspect ratio): min = " + str(min, 3L) + ", max = " + str(max, 3L) +
         ", median = " + str(median, 3L));
  }

  // Assign face colors (closest layer height index)
  void assign_face_colors(const std::vector<double> &areas,
                          const std::vector<double> &layer_heights)
  {
    // Iterate over ground faces
    face_colors.resize(areas.size());
    for (size_t i = 0; i < areas.size(); i++)
    {
      // Assign closest layer height index
      const double h = ideal_layer_height(areas[i]);
      const size_t face_color = closest_layer_height(h, layer_heights);
      face_colors[i] = face_color;
    }
  }

  // Reassign face colors to avoid big jumps
  void reassign_face_colors()
  {
    for (size_t i = 0; i < ff.size(); i++)
    {
      for (const auto &j : ff[i])
      {
        const auto diff = face_colors[i] - face_colors[j];
        if (diff > 1)
        {
          face_colors[i] -= diff - 1;
        }
      }
    }
  }

  // Check face colors (big jumps)
  size_t check_face_colors()
  {
    size_t num_big_jumps = 0;
    for (size_t i = 0; i < ff.size(); i++)
    {
      size_t _num_big_jumps = 0;
      for (const auto &j : ff[i])
      {
        const auto jump = face_colors[i] - face_colors[j];
        if (jump > 1)
          _num_big_jumps++;
      }
      if (_num_big_jumps > 0)
        num_big_jumps++;
    }

    double percentage = 100.0 * num_big_jumps / ff.size();
    info("Big jumps: " + str(num_big_jumps) + " / " + str(ff.size()) + " (" + str(percentage, 2L) +
         "%)");

    return num_big_jumps;
  }

  /// Assign vertex colors based on minimum neighbor face color
  void assign_vertex_colors(const std::vector<double> &layer_heights)
  {
    const size_t num_vertices = _ground_mesh.vertices.size();
    const size_t num_faces = _ground_mesh.faces.size();
    vertex_colors.resize(num_vertices, layer_heights.size());
    for (size_t i = 0; i < num_faces; i++)
    {
      const auto &face = _ground_mesh.faces[i];
      for (const auto &j : {face.v0, face.v1, face.v2})
        vertex_colors[j] = std::min(vertex_colors[j], face_colors[i]);
    }
  }

  // Sort faces by vertex color and index
  void sort_faces_by_vertex_color_and_index()
  {
    for (auto &face : _ground_mesh.faces)
    {
      size_t c0 = vertex_colors[face.v0];
      size_t c1 = vertex_colors[face.v1];
      size_t c2 = vertex_colors[face.v2];
      if (c0 > c1)
      {
        std::swap(c0, c1);
        std::swap(face.v0, face.v1);
      }
      if (c0 > c2)
      {
        std::swap(c0, c2);
        std::swap(face.v0, face.v2);
      }
      if (c1 > c2)
      {
        std::swap(c1, c2);
        std::swap(face.v1, face.v2);
      }
      // Compare and swap first and second elements
      if (c0 > c1 || (c0 == c1 && face.v0 > face.v1))
      {
        std::swap(c0, c1);
        std::swap(face.v0, face.v1);
      }

      // Compare and swap first and third elements
      if (c0 > c2 || (c0 == c2 && face.v0 > face.v2))
      {
        std::swap(c0, c2);
        std::swap(face.v0, face.v2);
      }

      // Compare and swap second and third elements
      if (c1 > c2 || (c1 == c2 && face.v1 > face.v2))
      {
        std::swap(c1, c2);
        std::swap(face.v1, face.v2);
      }
    }
  }

  // Assign face (prism) partitions based on face and vertex colors.
  //
  // Partition 0: 6 vertices, 3 tetrahedrons
  // Partition 1: 7 vertices, 4 tetrahedrons
  // Partition 2: 8 vertices, 5 tetrahedrons
  // Partition 3: 9 vertices, 6 tetrahedrons (reduntant)
  void assign_face_partitions()
  {
    const size_t num_faces = _ground_mesh.faces.size();
    face_partitions.resize(num_faces);
    for (size_t i = 0; i < num_faces; i++)
    {
      const auto &face = _ground_mesh.faces[i];
      const size_t c0 = vertex_colors[face.v0];
      const size_t c1 = vertex_colors[face.v1];
      const size_t c2 = vertex_colors[face.v2];
      face_partitions[i] = 3 * face_colors[i] - (c0 + c1 + c2);
    }
  }

  // Eliminate type 3 partitions which are just two stacked prisms of type 0
  void eliminate_type_3_partitions()
  {
    const size_t num_faces = _ground_mesh.faces.size();
    for (size_t i = 0; i < num_faces; i++)
    {
      if (face_partitions[i] == 3)
      {
        face_colors[i]--;
        face_partitions[i] = 0;
      }
    }
  }

  // Build mapping from vertices to faces
  void build_vertex_to_face_mapping()
  {
    const size_t num_faces = _ground_mesh.faces.size();
    vf.resize(num_faces);
    for (size_t i = 0; i < num_faces; i++)
    {
      vf[_ground_mesh.faces[i].v0].insert(i);
      vf[_ground_mesh.faces[i].v1].insert(i);
      vf[_ground_mesh.faces[i].v2].insert(i);
    }
  }

  /// Build mapping from faces to faces
  void build_face_to_face_mapping()
  {
    const size_t num_faces = _ground_mesh.faces.size();
    ff.resize(num_faces);
    for (size_t i = 0; i < num_faces; i++)
    {
      for (size_t j = 0; j < 3; j++)
      {
        for (const int &v : vf[_ground_mesh.faces[i][j]])
        {
          ff[i].insert(v);
        }
      }
    }
  }

  // Layer ground mesh
  VolumeMesh layer_ground_mesh(const std::vector<double> &layer_heights)
  {
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
    _column_mesh.num_max_layers = std::ceil(max_building_height / layer_heights.back());
    _column_mesh.num_min_layers = _column_mesh.num_max_layers << (layer_heights.size() - 1);

    // Layer vertices in columns
    for (size_t j = 0; j < _ground_mesh.vertices.size(); j++)
    {
      const Vector3D &vg = _ground_mesh.vertices[j];
      const double h = layer_heights[vertex_colors[j]];
      const size_t col_size = (_column_mesh.num_min_layers >> vertex_colors[j]) + 1;
      for (size_t i = 0; i < col_size; i++)
      {
        Vector3D v(vg.x, vg.y, i * h);
        _column_mesh.vertices[j].push_back(v);
      }
      _column_mesh.vertices_offset[j + 1] = _column_mesh.vertices_offset[j] + col_size;
    }

    // Layer cells in columns
    for (size_t i = 0; i < _ground_mesh.faces.size(); i++)
    {
      // Get face vertices and colors
      const Simplex2D face_simplex = _ground_mesh.faces[i];
      const std::array<size_t, 3> face = {face_simplex.v0, face_simplex.v1, face_simplex.v2};
      const std::array<int, 3> _vertex_colors = {vertex_colors[face_simplex.v0],
                                                 vertex_colors[face_simplex.v1],
                                                 vertex_colors[face_simplex.v2]};

      // Calculate offsets and sizes
      const std::array<size_t, 3> col_offsets = {0, 0, 0};
      const std::array<size_t, 3> col_sizes = {
          _column_mesh.vertices_offset[face[0] + 1] - _column_mesh.vertices_offset[face[0]],
          _column_mesh.vertices_offset[face[1] + 1] - _column_mesh.vertices_offset[face[1]],
          _column_mesh.vertices_offset[face[2] + 1] - _column_mesh.vertices_offset[face[2]]};

      // Calculate number of prisms
      const size_t num_prisms =
          (col_sizes.back() - 1) / (1 << (face_colors[i] - _vertex_colors[2]));
      _column_mesh.num_prisms[i] = num_prisms;

      // Create prism iterator
      std::vector<std::array<size_t, 4>> prism_iterator(num_prisms);
      for (size_t j = 0; j < num_prisms; j++)
      {
        const size_t layer_index = (j + 1) * _column_mesh.num_min_layers / num_prisms;
        prism_iterator[j] = {col_offsets[0] + j * (1 << (face_colors[i] - _vertex_colors[0])),
                             col_offsets[1] + j * (1 << (face_colors[i] - _vertex_colors[1])),
                             col_offsets[2] + j * (1 << (face_colors[i] - _vertex_colors[2])),
                             layer_index};
      }

      // Add tetrahedrons based on face partition type
      switch (face_partitions[i])
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

          _column_mesh.cells[i].emplace_back(K0);
          _column_mesh.cells[i].emplace_back(K1);
          _column_mesh.cells[i].emplace_back(K2);
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

          _column_mesh.cells[i].emplace_back(K0);
          _column_mesh.cells[i].emplace_back(K1);
          _column_mesh.cells[i].emplace_back(K2);
          _column_mesh.cells[i].emplace_back(K3);
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

          _column_mesh.cells[i].emplace_back(K0);
          _column_mesh.cells[i].emplace_back(K1);
          _column_mesh.cells[i].emplace_back(K2);
          _column_mesh.cells[i].emplace_back(K3);
          _column_mesh.cells[i].emplace_back(K4);
        }
      }
      break;

      default:
        error("Unhandled partition type: " + str(face_partitions[i]));
        break;
      }
    }

    // Set markers
    auto mesh_vertex_markers = face_to_vertex_markers();
    for (size_t i = 0; i < _ground_mesh.vertices.size(); i++)
    {
      // Initialize all markers to Other (-4)
      auto &markers = _column_mesh.markers[i];
      markers.resize(_column_mesh.vertices[i].size(), -4);

      // Set last (top) marker to Top (-3)
      markers.back() = -3;

      // Set first (bottom) marker to Ground (-2) if ground otherwise Halo (-1)
      markers.front() = mesh_vertex_markers[i] == -2 ? -2 : -1;
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

        // Mark top vertices as building
        if (_column_mesh.layer_index(cell.v0) == trimming_layer_indices[marker])
          _column_mesh.markers[cell.v0.column][cell.v0.index] = marker;
        if (_column_mesh.layer_index(cell.v1) == trimming_layer_indices[marker])
          _column_mesh.markers[cell.v1.column][cell.v1.index] = marker;
        if (_column_mesh.layer_index(cell.v2) == trimming_layer_indices[marker])
          _column_mesh.markers[cell.v2.column][cell.v2.index] = marker;
        if (_column_mesh.layer_index(cell.v3) == trimming_layer_indices[marker])
          _column_mesh.markers[cell.v3.column][cell.v3.index] = marker;
      }
    }

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

  /// Mapping face markers to Vertex markers for ground mesh.
  ///
  /// Each vertex adopts the max marker value from the faces it belongs to.
  ///
  /// Ground Mesh Markers:
  /// -2: Ground
  /// -1: Building halos
  ///  0: Building 0
  ///  1: Building 1
  ///  etc (non-negative integers mark faces inside buildings)
  std::vector<int> face_to_vertex_markers()
  {
    const size_t num_vertices = _ground_mesh.vertices.size();
    const size_t num_faces = _ground_mesh.faces.size();

    vertex_markers.resize(num_vertices, -2);

    if (!_ground_mesh.markers.size())
    {
      error("Ground mesh has no face Markers. Treating all "
            "faces as "
            "ground");
      return vertex_markers;
    }

    for (size_t f = 0; f < num_faces; f++)
    {
      if (_ground_mesh.markers[f] < -2)
      {
        info("Problem problem with marker:" + str(_ground_mesh.markers[f]));
      }

      const std::array<size_t, 3> I = {_ground_mesh.faces[f].v0, _ground_mesh.faces[f].v1,
                                       _ground_mesh.faces[f].v2};

      vertex_markers[I[0]] = std::max(vertex_markers[I[0]], _ground_mesh.markers[f]);
      vertex_markers[I[1]] = std::max(vertex_markers[I[1]], _ground_mesh.markers[f]);
      vertex_markers[I[2]] = std::max(vertex_markers[I[2]], _ground_mesh.markers[f]);
    }

    return vertex_markers;
  }

  VolumeMesh add_domain_padding(const std::vector<double> &layer_heights,
                                const double padding_height, double max_scale = 3.0)
  {
    // Get max layer height
    const double max_layer_height = layer_heights.back();

    // Number of Max height layers needed to cover padding height.
    const size_t n = static_cast<size_t>(2 * padding_height / (max_layer_height * (max_scale + 1)));

    info("Padding Domain with " + str(n) + " max layers scaled from 1 to " + str(max_scale));
    std::vector<size_t> offset_before_padding = _column_mesh.vertices_offset;

    for (size_t i = 0; i < _ground_mesh.vertices.size(); i++)
    {
      // Z coordinate of the last vertex of each column.
      const double top_vertex_z = _column_mesh.vertices[i].back().z;

      const size_t n_k = static_cast<size_t>(max_layer_height / layer_heights[vertex_colors[i]]);
      double layer_h = 0.0;
      const double s_1 = n > 1 ? (max_scale - 1) / (n - 1) : 0;
      for (size_t j = 1; j < n; j++)
      {
        // const double s = 1.0 + j*(max_scale - 1)/(n-1);
        const double s = 1.0 + j * s_1;

        for (size_t k = 0; k < n_k; k++)
        {
          layer_h += s * layer_heights[vertex_colors[i]];
          Vector3D v(_ground_mesh.vertices[i].x, _ground_mesh.vertices[i].y,
                     top_vertex_z + layer_h);
          _column_mesh.vertices[i].push_back(v);
        }
      }
      _column_mesh.vertices_offset[i + 1] =
          _column_mesh.vertices_offset[i] + _column_mesh.vertices[i].size();
    }

    size_t volume_mesh_num_cells = 0;
    for (size_t i = 0; i < _ground_mesh.faces.size(); i++)
    {
      const Simplex2D face_simplex = _ground_mesh.faces[i];
      std::array<size_t, 3> face = {face_simplex.v0, face_simplex.v1, face_simplex.v2};
      std::array<int, 3> v_colors = {vertex_colors[face_simplex.v0], vertex_colors[face_simplex.v1],
                                     vertex_colors[face_simplex.v2]};

      const std::array<size_t, 3> column_offsets = {
          offset_before_padding[face[0] + 1] - offset_before_padding[face[0]] - 1,
          offset_before_padding[face[1] + 1] - offset_before_padding[face[1]] - 1,
          offset_before_padding[face[2] + 1] - offset_before_padding[face[2]] - 1,
      };
      const std::array<size_t, 3> column_len = {
          _column_mesh.vertices_offset[face[0] + 1] - _column_mesh.vertices_offset[face[0]],
          _column_mesh.vertices_offset[face[1] + 1] - _column_mesh.vertices_offset[face[1]],
          _column_mesh.vertices_offset[face[2] + 1] - _column_mesh.vertices_offset[face[2]]};
      const size_t num_prisms =
          (column_len.back() - column_offsets.back() - 1) / (1 << (face_colors[i] - v_colors[2]));

      std::vector<std::array<size_t, 3>> prism_iterator(num_prisms);

      for (size_t j = 0; j < num_prisms; j++)
      {
        prism_iterator[j] = {column_offsets[0] + j * (1 << (face_colors[i] - v_colors[0])),
                             column_offsets[1] + j * (1 << (face_colors[i] - v_colors[1])),
                             column_offsets[2] + j * (1 << (face_colors[i] - v_colors[2]))};
      }
      switch (face_partitions[i])
      {
      case 0:
      {
        for (const auto &ar : prism_iterator)
        {
          const size_t k = ar[0];
          const size_t l = ar[1];
          const size_t m = ar[2];

          ColumnIndex bot_triangle_0(face[0], k), bot_triangle_1(face[1], l),
              bot_triangle_2(face[2], m);
          ColumnIndex top_triangle_0(face[0], k + 1), top_triangle_1(face[1], l + 1),
              top_triangle_2(face[2], m + 1);

          ColumnSimplex K0(bot_triangle_0, bot_triangle_1, bot_triangle_2, top_triangle_2, 0);
          ColumnSimplex K1(bot_triangle_0, top_triangle_1, bot_triangle_1, top_triangle_2, 0);
          ColumnSimplex K2(bot_triangle_0, top_triangle_0, top_triangle_1, top_triangle_2, 0);

          _column_mesh.cells[i].emplace_back(K0);
          _column_mesh.cells[i].emplace_back(K1);
          _column_mesh.cells[i].emplace_back(K2);
          volume_mesh_num_cells += 3;
        }
      }
      break;
      case 1:
      {
        for (const auto &ar : prism_iterator)
        {
          const size_t k = ar[0];
          const size_t l = ar[1];
          const size_t m = ar[2];

          ColumnIndex bot_triangle_0(face[0], k), bot_triangle_1(face[1], l),
              bot_triangle_2(face[2], m);
          ColumnIndex mid_triangle_0(face[0], k + 1);
          ColumnIndex top_triangle_0(face[0], k + 2), top_triangle_1(face[1], l + 1),
              top_triangle_2(face[2], m + 1);

          ColumnSimplex K0(bot_triangle_0, bot_triangle_1, bot_triangle_2, mid_triangle_0, 0);
          ColumnSimplex K1(bot_triangle_1, top_triangle_2, bot_triangle_2, mid_triangle_0, 0);
          ColumnSimplex K2(bot_triangle_1, top_triangle_2, mid_triangle_0, top_triangle_1, 0);
          ColumnSimplex K3(mid_triangle_0, top_triangle_0, top_triangle_1, top_triangle_2, 0);

          _column_mesh.cells[i].emplace_back(K0);
          _column_mesh.cells[i].emplace_back(K1);
          _column_mesh.cells[i].emplace_back(K2);
          _column_mesh.cells[i].emplace_back(K3);
          volume_mesh_num_cells += 4;
        }
      }
      break;
      case 2:
      {
        for (const auto &ar : prism_iterator)
        {
          const size_t k = ar[0];
          const size_t l = ar[1];
          const size_t m = ar[2];

          ColumnIndex bot_triangle_0(face[0], k), bot_triangle_1(face[1], l),
              bot_triangle_2(face[2], m);
          ColumnIndex mid_triangle_0(face[0], k + 1), mid_triangle_1(face[1], l + 1);
          ColumnIndex top_triangle_0(face[0], k + 2), top_triangle_1(face[1], l + 2),
              top_triangle_2(face[2], m + 1);

          ColumnSimplex K0(bot_triangle_0, bot_triangle_1, bot_triangle_2, mid_triangle_1, 0);
          ColumnSimplex K1(bot_triangle_0, bot_triangle_2, mid_triangle_0, mid_triangle_1, 0);
          ColumnSimplex K2(bot_triangle_2, top_triangle_2, mid_triangle_0, mid_triangle_1, 0);
          ColumnSimplex K3(top_triangle_0, top_triangle_2, top_triangle_1, mid_triangle_0, 0);
          ColumnSimplex K4(top_triangle_1, top_triangle_2, mid_triangle_1, mid_triangle_0, 0);

          _column_mesh.cells[i].emplace_back(K0);
          _column_mesh.cells[i].emplace_back(K1);
          _column_mesh.cells[i].emplace_back(K2);
          _column_mesh.cells[i].emplace_back(K3);
          _column_mesh.cells[i].emplace_back(K4);
          volume_mesh_num_cells += 5;
        }
      }
      break;

      case 3:
      {
        for (const auto &ar : prism_iterator)
        {
          const size_t k = ar[0];
          const size_t l = ar[1];
          const size_t m = ar[2];

          ColumnIndex bot_triangle_0(face[0], k), bot_triangle_1(face[1], l),
              bot_triangle_2(face[2], m);
          ColumnIndex mid_triangle_0(face[0], k + 1), mid_triangle_1(face[1], l + 1),
              mid_triangle_2(face[2], m + 1);
          ColumnIndex top_triangle_0(face[0], k + 2), top_triangle_1(face[1], l + 2),
              top_triangle_2(face[2], m + 2);

          ColumnSimplex K0(bot_triangle_0, bot_triangle_1, bot_triangle_2, mid_triangle_2, 0);
          ColumnSimplex K1(bot_triangle_0, mid_triangle_1, bot_triangle_1, mid_triangle_2, 0);
          ColumnSimplex K2(bot_triangle_0, mid_triangle_0, mid_triangle_1, mid_triangle_2, 0);
          ColumnSimplex K3(mid_triangle_0, mid_triangle_1, mid_triangle_2, top_triangle_2, 0);
          ColumnSimplex K4(mid_triangle_0, top_triangle_1, mid_triangle_1, top_triangle_2, 0);
          ColumnSimplex K5(mid_triangle_0, top_triangle_0, top_triangle_1, top_triangle_2, 0);

          _column_mesh.cells[i].emplace_back(K0);
          _column_mesh.cells[i].emplace_back(K1);
          _column_mesh.cells[i].emplace_back(K2);
          _column_mesh.cells[i].emplace_back(K3);
          _column_mesh.cells[i].emplace_back(K4);
          _column_mesh.cells[i].emplace_back(K5);
          volume_mesh_num_cells += 6;
        }
        break;
      }
      default:
        error("Face Coloring Error: Large layer height "
              "difference");
        break;
      }
    }
    info("Cells added to the Volume Mesh during domain padding: " + str(volume_mesh_num_cells));
    return _column_mesh.to_volume_mesh();
  }
};

} // namespace DTCC_BUILDER

#endif
