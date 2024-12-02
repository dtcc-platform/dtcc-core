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
#include "model/City.h"
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

  std::vector<int> max_building_colors;

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
    VolumeMesh volume_mesh = layer_ground_mesh(layer_heights);
    t3_2.stop();
    t3_2.print();

    // Debugging
    if (debug_step <= 2)
      return volume_mesh;

    // Volume mesh smoothing
    Timer t3_3("Step 3.3: Volume mesh smoothing (ground only)");
    volume_mesh = Smoother::smooth_volume_mesh(volume_mesh, _buildings, _dem, 0.0, false, false,
                                               smoother_iterations, smoother_relative_tolerance);
    t3_3.stop();
    t3_3.print();

    // Debugging
    if (debug_step == 3)
      return volume_mesh;

    _column_mesh._update_vertices(volume_mesh);

    // Trim volume mesh
    Timer t3_4("Step 3.4: Volume mesh trimming");
    volume_mesh = trim_volume_mesh(layer_heights);
    t3_4.stop();
    t3_4.print();

    // Debugging
    if (debug_step == 4)
      return volume_mesh;

    // Compute top height
    const double top_height = domain_height + _dem.max();

    // Smooth volume mesh (again)
    Timer t3_5("Step 3.5: Volume mesh smoothing (ground and buildings)");
    volume_mesh =
        Smoother::smooth_volume_mesh(volume_mesh, _buildings, _dem, top_height, true, true,
                                     smoother_iterations, smoother_relative_tolerance);
    t3_5.stop();
    t3_5.print();

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
    for (int i = 0; i < steps + 1; i++)
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
    sort_faces_by_vertex_color();

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

  // Assign face colors (closest layer height index)
  void assign_face_colors(const std::vector<double> &areas,
                          const std::vector<double> &layer_heights)
  {
    // Array of common building colors
    std::vector<size_t> building_colors(_buildings.size(), layer_heights.size());

    // Iterate over ground faces
    face_colors.resize(areas.size());
    for (size_t i = 0; i < areas.size(); i++)
    {
      // Assign closest layer height index
      const double h = ideal_layer_height(areas[i]);
      const size_t face_color = closest_layer_height(h, layer_heights);
      face_colors[i] = face_color;

      // Save minimum color for each building
      const int marker = _ground_mesh.markers[i];
      if (marker >= 0)
        building_colors[marker] = std::min(building_colors[marker], face_color);
    }

    // FIXME: Is this necessary? Much better quality if we do this, but why?

    // Assign common colors to all faces of the same building
    // for (size_t i = 0; i < areas.size(); i++)
    // {
    //   const int marker = _ground_mesh.markers[i];
    //   if (marker >= 0)
    //     face_colors[i] = building_colors[marker];
    // }
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

  // Sort faces by vertex color
  void sort_faces_by_vertex_color()
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
    double max_height = 0.0;
    for (size_t i = 0; i < _buildings.size(); i++)
    {
      const double height = _buildings[i].max_height() - _building_ground_height[i];
      max_height = std::max(max_height, height);
    }

    // Compute layering height (cover tallest building and respect layer heights)
    const double H = layer_heights.back();
    const double layering_height = std::ceil(max_height / H) * H;

    // Layer vertices in columns
    for (size_t j = 0; j < _ground_mesh.vertices.size(); j++)
    {
      const Vector3D &vg = _ground_mesh.vertices[j];
      const double h = layer_heights[vertex_colors[j]];
      const size_t col_size = static_cast<size_t>((layering_height / h)) + 1;
      std::cout << "Column size: " << col_size << std::endl;
      std::cout << "Column height: " << h * (col_size - 1) << std::endl;

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

      // Create prism iterator
      const size_t num_prisms =
          (col_sizes.back() - 1) / (1 << (face_colors[i] - _vertex_colors[2]));
      std::vector<std::array<size_t, 3>> prism_iterator(num_prisms);
      for (size_t j = 0; j < num_prisms; j++)
      {
        prism_iterator[j] = {col_offsets[0] + j * (1 << (face_colors[i] - _vertex_colors[0])),
                             col_offsets[1] + j * (1 << (face_colors[i] - _vertex_colors[1])),
                             col_offsets[2] + j * (1 << (face_colors[i] - _vertex_colors[2]))};
      }

      // Add tetrahedrons based on face partition type
      switch (face_partitions[i])
      {
      case 0:
      {
        for (const auto &ar : prism_iterator)
        {
          const size_t k = ar[0];
          const size_t l = ar[1];
          const size_t m = ar[2];

          ColumnIndex bot_triangle_0(face[0], k);
          ColumnIndex bot_triangle_1(face[1], l);
          ColumnIndex bot_triangle_2(face[2], m);
          ColumnIndex top_triangle_0(face[0], k + 1);
          ColumnIndex top_triangle_1(face[1], l + 1);
          ColumnIndex top_triangle_2(face[2], m + 1);

          ColumnSimplex K0(bot_triangle_0, bot_triangle_1, bot_triangle_2, top_triangle_2);
          ColumnSimplex K1(bot_triangle_0, top_triangle_1, bot_triangle_1, top_triangle_2);
          ColumnSimplex K2(bot_triangle_0, top_triangle_0, top_triangle_1, top_triangle_2);

          _column_mesh.cells[i].emplace_back(K0);
          _column_mesh.cells[i].emplace_back(K1);
          _column_mesh.cells[i].emplace_back(K2);
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

          ColumnIndex bot_triangle_0(face[0], k);
          ColumnIndex bot_triangle_1(face[1], l);
          ColumnIndex bot_triangle_2(face[2], m);
          ColumnIndex mid_triangle_0(face[0], k + 1);
          ColumnIndex top_triangle_0(face[0], k + 2);
          ColumnIndex top_triangle_1(face[1], l + 1);
          ColumnIndex top_triangle_2(face[2], m + 1);

          ColumnSimplex K0(bot_triangle_0, bot_triangle_1, bot_triangle_2, mid_triangle_0);
          ColumnSimplex K1(bot_triangle_1, top_triangle_2, bot_triangle_2, mid_triangle_0);
          ColumnSimplex K2(bot_triangle_1, top_triangle_2, mid_triangle_0, top_triangle_1);
          ColumnSimplex K3(mid_triangle_0, top_triangle_0, top_triangle_1, top_triangle_2);

          _column_mesh.cells[i].emplace_back(K0);
          _column_mesh.cells[i].emplace_back(K1);
          _column_mesh.cells[i].emplace_back(K2);
          _column_mesh.cells[i].emplace_back(K3);
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

          ColumnIndex bot_triangle_0(face[0], k);
          ColumnIndex bot_triangle_1(face[1], l);
          ColumnIndex bot_triangle_2(face[2], m);
          ColumnIndex mid_triangle_0(face[0], k + 1);
          ColumnIndex mid_triangle_1(face[1], l + 1);
          ColumnIndex top_triangle_0(face[0], k + 2);
          ColumnIndex top_triangle_1(face[1], l + 2);
          ColumnIndex top_triangle_2(face[2], m + 1);

          ColumnSimplex K0(bot_triangle_0, bot_triangle_1, bot_triangle_2, mid_triangle_1);
          ColumnSimplex K1(bot_triangle_0, bot_triangle_2, mid_triangle_0, mid_triangle_1);
          ColumnSimplex K2(bot_triangle_2, top_triangle_2, mid_triangle_0, mid_triangle_1);
          ColumnSimplex K3(top_triangle_0, top_triangle_2, top_triangle_1, mid_triangle_0);
          ColumnSimplex K4(top_triangle_1, top_triangle_2, mid_triangle_1, mid_triangle_0);

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

  // Mapping face markers to Vertex markers for ground mesh.
  //
  // Each vertex adopts the min marker value from the faces it belongs to.
  //
  // Ground Mesh Markers:
  // -2: Ground
  // -1: Building halos
  //  0: Building 0
  //  1: Building 1
  //  etc (non-negative integers mark faces inside buildings)
  std::vector<int> face_to_vertex_markers_min()
  {
    vertex_markers.resize(_ground_mesh.vertices.size(), std::numeric_limits<int>::max());

    if (!_ground_mesh.markers.size())
    {
      error("Ground mesh has no face Markers. Treating all "
            "faces as "
            "ground");
      return vertex_markers;
    }

    for (size_t f = 0; f < _ground_mesh.markers.size(); f++)
    {
      if (_ground_mesh.markers[f] < -2)
      {
        info("Problem problem with marker:" + str(_ground_mesh.markers[f]));
      }

      const std::array<size_t, 3> I = {_ground_mesh.faces[f].v0, _ground_mesh.faces[f].v1,
                                       _ground_mesh.faces[f].v2};

      vertex_markers[I[0]] = std::min(vertex_markers[I[0]], _ground_mesh.markers[f]);
      vertex_markers[I[1]] = std::min(vertex_markers[I[1]], _ground_mesh.markers[f]);
      vertex_markers[I[2]] = std::min(vertex_markers[I[2]], _ground_mesh.markers[f]);
    }

    return vertex_markers;
  }

  /// Adjust the input building heights to the maximum layer height used by
  /// the vertices that belong to the building.
  std::vector<double> adjust_building_heights(const std::vector<double> &layer_heights)
  {
    const size_t num_buildings = _buildings.size();
    max_building_colors.resize(num_buildings, 0);
    // build map from buildings to cells in 2D mesh

    info("Building map from buildings to cells in 2D mesh");
    std::vector<std::vector<size_t>> building_faces(num_buildings);
    for (size_t i = 0; i < _ground_mesh.markers.size(); i++)
    {
      const int building_index = _ground_mesh.markers[i];
      if (building_index >= 0)
      {
        building_faces[building_index].push_back(i);
      }
    }

    for (size_t building_index = 0; building_index < num_buildings; building_index++)
    {
      // Note: We could either work with face or vertex
      // colors/layer-heights
      for (const size_t &f : building_faces[building_index])
      {
        // Working with face colors... This choice also
        // changes the calculation for the number of prisms
        // created.
        max_building_colors[building_index] =
            std::max({face_colors[f], max_building_colors[building_index]});
      }
    }

    // Loop over buildings
    std::vector<double> adj_building_heights(num_buildings, 0);
    for (size_t building_index = 0; building_index < num_buildings; building_index++)
    {
      const double fitting_layer_height = layer_heights[max_building_colors[building_index]];

      double n =
          (_buildings[building_index].max_height() - _building_ground_height[building_index]) /
          fitting_layer_height; // CHANGE ASAP WROOONG

      adj_building_heights[building_index] = fitting_layer_height;
      if (n > 1)
      {
        adj_building_heights[building_index] *= std::floor(n);
      }
    }

    std::cout << "Adjusted Building Heights Size: " << adj_building_heights.size() << std::endl;
    return adj_building_heights;
  }

  std::vector<size_t> get_building_trimming_index(double buffer = 0.0)
  {
    const size_t num_buildings = _buildings.size();
    std::vector<size_t> trimming_index(num_buildings, 0);

    // Looping over all Vertex columns of the columnar volume mesh.
    for (size_t i = 0; i < _ground_mesh.vertices.size(); i++)
    {
      int marker = vertex_markers[i];
      if (marker >= 0 && vertex_colors[i] == max_building_colors[marker])
      {
        size_t j = 0;
        while (_column_mesh.vertices[i][j + 1].z <= (_buildings[marker].max_height() + buffer))
        {
          j++;
        }
        trimming_index[marker] =
            std::max(trimming_index[marker], static_cast<size_t>(j * (1 << vertex_colors[i])));
      }
    }

    for (size_t i = 0; i < trimming_index.size(); i++)
    {
      if (trimming_index[i] == std::numeric_limits<size_t>::max())
      {
        trimming_index[i] = 1 << max_building_colors[i];
      }
    }
    return trimming_index;
  }

  std::vector<double> get_adjusted_building_heights(const std::vector<double> &layer_heights)
  {
    const size_t num_buildings = _buildings.size();
    auto trimming_index = get_building_trimming_index();

    std::vector<double> adj_b_heights(num_buildings, 0.0);
    for (size_t i = 0; i < num_buildings; i++)
    {
      const double building_height = _buildings[i].max_height() - _building_ground_height[i];
      adj_b_heights[i] = trimming_index[i] * layer_heights[max_building_colors[i]];

      if ((adj_b_heights[i] - building_height) > layer_heights[max_building_colors[i]])
        adj_b_heights[i] -= layer_heights[max_building_colors[i]];
      else if ((building_height - adj_b_heights[i]) > layer_heights[max_building_colors[i]])
        adj_b_heights[i] += layer_heights[max_building_colors[i]];
    }

    return adj_b_heights;
  }

  VolumeMesh trim_volume_mesh(const std::vector<double> &layer_heights)
  {
    info("Trimming volume mesh...");

    // Adjust Building heights to the largest layer height assigned to
    // buildings faces.
    std::vector<double> adj_building_heights = adjust_building_heights(layer_heights);
    std::vector<int> min_vertex_markers = face_to_vertex_markers_min();
    std::cout << "Short test to find error. Ajusted buisdfsdfflding "
                 "Heights: "
              << std::endl;

    auto trimming_index = get_building_trimming_index();

    // In this vector we store how many minimum layer vertices in each
    // column are skipped to trim mesh and create a building
    vert_skip.resize(_ground_mesh.vertices.size(), 0);
    for (size_t i = 0; i < _ground_mesh.vertices.size(); i++)
    {
      if (min_vertex_markers[i] >= 0)
      {

        const double building_adj_height = adj_building_heights[min_vertex_markers[i]];
        const double layer_h = layer_heights[vertex_colors[i]];

        size_t col_begin = static_cast<size_t>((building_adj_height / layer_h));
        vert_skip[i] = col_begin;

        std::vector<Vector3D> trimmed_column;
        const size_t col_end = _column_mesh.vertices[i].size();
        for (size_t j = col_begin; j < col_end; j++)
        {
          Vector3D v = _column_mesh.vertices[i][j];
          trimmed_column.push_back(v);
        }
        _column_mesh.vertices[i] = trimmed_column;
      }
      _column_mesh.vertices_offset[i + 1] =
          _column_mesh.vertices_offset[i] + _column_mesh.vertices[i].size();
    }

    std::cout << "Short test to find error. Ajusted building Heights: "
              << adj_building_heights.size() << std::endl;

    for (size_t i = 0; i < _ground_mesh.markers.size(); i++)
    {
      // No trimming needed for cell columns that are not marked
      // as buildings.
      if (_ground_mesh.markers[i] < 0)
      {
        continue;
      }

      _column_mesh.cells[i].clear();
      const Simplex2D face_simplex = _ground_mesh.faces[i];
      const int face_marker = _ground_mesh.markers[i];
      std::array<size_t, 3> face = {face_simplex.v0, face_simplex.v1, face_simplex.v2};

      std::array<int, 3> v_colors = {vertex_colors[face_simplex.v0], vertex_colors[face_simplex.v1],
                                     vertex_colors[face_simplex.v2]};

      const std::array<size_t, 3> column_offsets = {0, 0, 0};

      const std::array<size_t, 3> col_skip = {vert_skip[face[0]] * (1 << v_colors[0]),
                                              vert_skip[face[1]] * (1 << v_colors[1]),
                                              vert_skip[face[2]] * (1 << v_colors[2])};

      const std::array<size_t, 3> column_len = {
          _column_mesh.vertices_offset[face[0] + 1] - _column_mesh.vertices_offset[face[0]],
          _column_mesh.vertices_offset[face[1] + 1] - _column_mesh.vertices_offset[face[1]],
          _column_mesh.vertices_offset[face[2] + 1] - _column_mesh.vertices_offset[face[2]]};

      const std::array<size_t, 3> column_len_nrml = {(column_len[0] - 1) * (1 << v_colors[0]),
                                                     (column_len[1] - 1) * (1 << v_colors[1]),
                                                     (column_len[2] - 1) * (1 << v_colors[2])};

      size_t min_nrml_col_len =
          std::min({column_len_nrml[0], column_len_nrml[1], column_len_nrml[2]});

      size_t max_col_skip = std::max({col_skip[0], col_skip[1], col_skip[2]});

      size_t face_col_skip =
          static_cast<size_t>((adj_building_heights[face_marker] / layer_heights[0]));

      // We do that for the corner of the buildings. Although no
      // vertices are trimmed we should not create cells because
      // we eliminate the corners of the buildings.
      min_nrml_col_len -= face_col_skip - max_col_skip;
      max_col_skip = face_col_skip; //(trimming_index[face_marker])* 1
                                    //<<(max_building_colors[face_marker]);

      const std::array<size_t, 3> prism_start_indexes{
          static_cast<size_t>((max_col_skip - col_skip[0]) / (1 << v_colors[0])),
          static_cast<size_t>((max_col_skip - col_skip[1]) / (1 << v_colors[1])),
          static_cast<size_t>((max_col_skip - col_skip[2]) / (1 << v_colors[2]))};

      const size_t num_prisms = static_cast<size_t>(min_nrml_col_len / (1 << face_colors[i]));

      std::vector<std::array<size_t, 3>> prism_iterator(num_prisms);

      for (size_t j = 0; j < num_prisms; j++)
      {
        prism_iterator[j] = {
            column_offsets[0] + prism_start_indexes[0] + (j << (face_colors[i] - v_colors[0])),
            column_offsets[1] + prism_start_indexes[1] + (j << (face_colors[i] - v_colors[1])),
            column_offsets[2] + prism_start_indexes[2] + (j << (face_colors[i] - v_colors[2]))};
      }

      // FIXME: Reuse code above in layer_ground_mesh

      switch (face_partitions[i])
      {
      case 0:
      {
        for (const auto &ar : prism_iterator)
        {
          const size_t k = ar[0];
          const size_t l = ar[1];
          const size_t m = ar[2];

          ColumnIndex bot_triangle_0(face[0], k);
          ColumnIndex bot_triangle_1(face[1], l);
          ColumnIndex bot_triangle_2(face[2], m);
          ColumnIndex top_triangle_0(face[0], k + 1);
          ColumnIndex top_triangle_1(face[1], l + 1);
          ColumnIndex top_triangle_2(face[2], m + 1);

          ColumnSimplex K0(bot_triangle_0, bot_triangle_1, bot_triangle_2, top_triangle_2);
          ColumnSimplex K1(bot_triangle_0, top_triangle_1, bot_triangle_1, top_triangle_2);
          ColumnSimplex K2(bot_triangle_0, top_triangle_0, top_triangle_1, top_triangle_2);

          _column_mesh.cells[i].emplace_back(K0);
          _column_mesh.cells[i].emplace_back(K1);
          _column_mesh.cells[i].emplace_back(K2);
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

          ColumnIndex bot_triangle_0(face[0], k);
          ColumnIndex bot_triangle_1(face[1], l);
          ColumnIndex bot_triangle_2(face[2], m);
          ColumnIndex mid_triangle_0(face[0], k + 1);
          ColumnIndex top_triangle_0(face[0], k + 2);
          ColumnIndex top_triangle_1(face[1], l + 1);
          ColumnIndex top_triangle_2(face[2], m + 1);

          ColumnSimplex K0(bot_triangle_0, bot_triangle_1, bot_triangle_2, mid_triangle_0);
          ColumnSimplex K1(bot_triangle_1, top_triangle_2, bot_triangle_2, mid_triangle_0);
          ColumnSimplex K2(bot_triangle_1, top_triangle_2, mid_triangle_0, top_triangle_1);
          ColumnSimplex K3(mid_triangle_0, top_triangle_0, top_triangle_1, top_triangle_2);

          _column_mesh.cells[i].emplace_back(K0);
          _column_mesh.cells[i].emplace_back(K1);
          _column_mesh.cells[i].emplace_back(K2);
          _column_mesh.cells[i].emplace_back(K3);
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

          ColumnIndex bot_triangle_0(face[0], k);
          ColumnIndex bot_triangle_1(face[1], l);
          ColumnIndex bot_triangle_2(face[2], m);
          ColumnIndex mid_triangle_0(face[0], k + 1);
          ColumnIndex mid_triangle_1(face[1], l + 1);
          ColumnIndex top_triangle_0(face[0], k + 2);
          ColumnIndex top_triangle_1(face[1], l + 2);
          ColumnIndex top_triangle_2(face[2], m + 1);

          ColumnSimplex K0(bot_triangle_0, bot_triangle_1, bot_triangle_2, mid_triangle_1);
          ColumnSimplex K1(bot_triangle_0, bot_triangle_2, mid_triangle_0, mid_triangle_1);
          ColumnSimplex K2(bot_triangle_2, top_triangle_2, mid_triangle_0, mid_triangle_1);
          ColumnSimplex K3(top_triangle_0, top_triangle_2, top_triangle_1, mid_triangle_0);
          ColumnSimplex K4(top_triangle_1, top_triangle_2, mid_triangle_1, mid_triangle_0);

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

    const auto mesh_vertex_markers = face_to_vertex_markers();

    for (size_t i = 0; i < _ground_mesh.vertices.size(); i++)
    {
      const int marker = mesh_vertex_markers[i];
      std::vector<int> tmp(_column_mesh.vertices[i].size(), -4);
      if (marker < 0)
      {
        tmp.front() = marker;
      }
      else
      {
        const double layer_h = layer_heights[vertex_colors[i]];
        const size_t idx =
            static_cast<size_t>((adj_building_heights[marker] / layer_h)) - vert_skip[i];
        if (idx == 0)
        {
          tmp.front() = marker;
        }
        else
        {
          tmp.front() = -1;
          tmp[idx] = marker;
        }
        tmp.back() = -3;
      }
      _column_mesh.markers[i] = tmp;
    }

    return _column_mesh.to_volume_mesh();
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

          ColumnSimplex K0(bot_triangle_0, bot_triangle_1, bot_triangle_2, top_triangle_2);
          ColumnSimplex K1(bot_triangle_0, top_triangle_1, bot_triangle_1, top_triangle_2);
          ColumnSimplex K2(bot_triangle_0, top_triangle_0, top_triangle_1, top_triangle_2);

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

          ColumnSimplex K0(bot_triangle_0, bot_triangle_1, bot_triangle_2, mid_triangle_0);
          ColumnSimplex K1(bot_triangle_1, top_triangle_2, bot_triangle_2, mid_triangle_0);
          ColumnSimplex K2(bot_triangle_1, top_triangle_2, mid_triangle_0, top_triangle_1);
          ColumnSimplex K3(mid_triangle_0, top_triangle_0, top_triangle_1, top_triangle_2);

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

          ColumnSimplex K0(bot_triangle_0, bot_triangle_1, bot_triangle_2, mid_triangle_1);
          ColumnSimplex K1(bot_triangle_0, bot_triangle_2, mid_triangle_0, mid_triangle_1);
          ColumnSimplex K2(bot_triangle_2, top_triangle_2, mid_triangle_0, mid_triangle_1);
          ColumnSimplex K3(top_triangle_0, top_triangle_2, top_triangle_1, mid_triangle_0);
          ColumnSimplex K4(top_triangle_1, top_triangle_2, mid_triangle_1, mid_triangle_0);

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

          ColumnSimplex K0(bot_triangle_0, bot_triangle_1, bot_triangle_2, mid_triangle_2);
          ColumnSimplex K1(bot_triangle_0, mid_triangle_1, bot_triangle_1, mid_triangle_2);
          ColumnSimplex K2(bot_triangle_0, mid_triangle_0, mid_triangle_1, mid_triangle_2);
          ColumnSimplex K3(mid_triangle_0, mid_triangle_1, mid_triangle_2, top_triangle_2);
          ColumnSimplex K4(mid_triangle_0, top_triangle_1, mid_triangle_1, top_triangle_2);
          ColumnSimplex K5(mid_triangle_0, top_triangle_0, top_triangle_1, top_triangle_2);

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
