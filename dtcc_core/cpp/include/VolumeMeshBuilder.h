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

#include "meshing/DomainPadding.h"
#include "meshing/MeshImprovement.h"

namespace DTCC_BUILDER
{



class BuilderMesh
{
  public:
    // References to the original mesh's attributes.

    /// Array of vertices
    std::vector<Vector3D> &vertices;

    /// Array of faces (triangles)
    std::vector<Simplex2D> &faces;

    /// Array of cell markers
    std::vector<int> &markers;

    // Stores vertex to face mapping for ground mesh.
    std::vector<std::unordered_set<int>> vf;

    // Stores vertex to face mapping for ground mesh.
    std::vector<std::unordered_set<int>> ff;

    std::vector<int> face_colors;

    std::vector<int> vertex_colors;

    std::vector<int> face_partitions;



    // BuilderMesh(){}

    BuilderMesh(Mesh &mesh): 
          faces(mesh.faces), 
          vertices(mesh.vertices), 
          markers(mesh.markers),
          face_colors(mesh.faces.size()),
          vertex_colors(mesh.vertices.size(),mesh.faces.size()),
          face_partitions(mesh.faces.size())
    {}

    ~BuilderMesh(){}

    // Build mapping from vertices to faces
  void build_vertex_to_face_mapping()
  {
    const size_t num_faces = faces.size();
    vf.resize(num_faces);
    for (size_t i = 0; i < num_faces; i++)
    {
      vf[faces[i].v0].insert(i);
      vf[faces[i].v1].insert(i);
      vf[faces[i].v2].insert(i);
    }
  }

  /// Build mapping from faces to faces
  void build_face_to_face_mapping()
  {
    const size_t num_faces = faces.size();
    ff.resize(num_faces);
    for (size_t i = 0; i < num_faces; i++)
    {
      for (size_t j = 0; j < 3; j++)
      {
        for (const int &v : vf[faces[i][j]])
        {
          ff[i].insert(v);
        }
      }
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
    const size_t num_vertices = vertices.size();
    const size_t num_faces = faces.size();
    
    std::vector<int> vertex_markers(num_vertices, -2);

    if (!markers.size())
    {
      error("Ground mesh has no face Markers. Treating all "
            "faces as "
            "ground");
      return vertex_markers;
    }

    for (size_t f = 0; f < num_faces; f++)
    {
      if (markers[f] < -2)
      {
        error(" Not allowd value for marker:" + str(markers[f]));
      }

      const std::array<size_t, 3> I = {faces[f].v0,faces[f].v1,
                                       faces[f].v2};

      vertex_markers[I[0]] = std::max(vertex_markers[I[0]], markers[f]);
      vertex_markers[I[1]] = std::max(vertex_markers[I[1]], markers[f]);
      vertex_markers[I[2]] = std::max(vertex_markers[I[2]], markers[f]);
    }

    return vertex_markers;
  }

};





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

  // // Stores vertex to face mapping for ground mesh.
  // std::vector<std::unordered_set<int>> vf;

  // // Stores vertex to face mapping for ground mesh.
  // std::vector<std::unordered_set<int>> ff;

  
  Vector2D mesh_center;
  // std::vector<int> face_colors;
  // std::vector<int> vertex_colors;

  // std::vector<int> face_partitions;

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
    mesh_center = Vector2D((mesh_bounds.P.x + mesh_bounds.Q.x) / 2.0,
                              (mesh_bounds.P.y + mesh_bounds.Q.y) / 2.0);


    
    compute_building_ground_heights();
  }

  // Destructor
  ~VolumeMeshBuilder() {}

  VolumeMesh build(const size_t smoother_iterations = 1000,
                   const double smoother_relative_tolerance = 0.001,
                   const double domain_padding_height = 0.0, 
                   const double aspect_ratio_threshold = 10.0,
                   const size_t debug_step = 6)
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
    VolumeMesh volume_mesh = layer_ground_mesh(_ground_mesh,layer_heights);
    t3_2.stop();
    t3_2.print();
    check_mesh_quality(volume_mesh, 2);

    // Debugging
    if (debug_step == 2)
      return volume_mesh;

    // Volume mesh smoothing
    Timer t3_3("Step 3.3: Volume mesh smoothing (ground only)");
    info("Smoothing volume mesh...");
    const bool fix_top = false;
    
    std::cout << "Ground Mesh bounding box: " << mesh_bounds.__str__() << std::endl;
    volume_mesh = Smoother::smooth_volume_mesh_elastic(volume_mesh, _buildings, _dem, 0.0, false, fix_top,
                                               smoother_iterations, smoother_relative_tolerance, mesh_bounds);
    // volume_mesh = Smoother::smooth_volume_mesh_amgcl(volume_mesh, _buildings, _dem, 0.0, false, false,
    //                                            smoother_iterations, smoother_relative_tolerance);
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
    volume_mesh = Smoother::smooth_volume_mesh_elastic(volume_mesh, _buildings, _dem, 0.0, true, false,
                                               smoother_iterations, smoother_relative_tolerance, mesh_bounds);
    t3_5.stop();
    t3_5.print();
    check_mesh_quality(volume_mesh, 5);
    
    if (debug_step == 5)
      return volume_mesh;
    
    Timer t3_6("Step 3.6: Volume mesh improvement (Removing slivers)");
    info("Refining volume mesh...");
    volume_mesh = VolumeMeshImprovement::remove_tetrahedra(volume_mesh,aspect_ratio_threshold);
    t3_6.stop();
    t3_6.print();
    check_mesh_quality(volume_mesh, 6);
    
    if (debug_step == 6)
      return volume_mesh;
    
    Timer t3_7("Step 3.7: Add Volume mesh Domain padding");
    info("Adding Domain padding...");

    top_height = domain_height + _dem.max();
    volume_mesh = layer_padding_mesh(volume_mesh,domain_height);
    t3_7.stop();
    t3_7.print();
    check_mesh_quality(volume_mesh, 7);

    return volume_mesh;
  }

private:
  // Compute layer heights for all faces in the ground mesh
  std::vector<double> compute_layer_heights(BuilderMesh &mesh)
  {
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
    assign_face_colors(areas, layer_heights, mesh);
    check_layer_heights(areas,layer_heights, mesh);

    // Build vertex and face mappings
    info("Building mapping from vertices to faces...");
    mesh.build_vertex_to_face_mapping();
    info("Building mapping from faces to faces...");
    mesh.build_face_to_face_mapping();

    // Iteratively reassign colors to avoid big jumps
    info("Reassigning colors to avoid big jumps...");
    size_t iteration = 0;
    const size_t max_color_iterations = 10;
    while (check_face_colors(mesh) > 0)
    {
      reassign_face_colors(mesh);
      if (++iteration == max_color_iterations)
      {
        error("Reached max color iterations");
      }
    }
    check_layer_heights(areas,layer_heights,mesh);

    // Assign vertex colors and sort by color
    info("Assigning vertex colors...");
    assign_vertex_colors(mesh);
    sort_faces_by_vertex_color_and_index(mesh);

    // Assign face partitions
    info("Assigning face partitions...");
    assign_face_partitions(mesh);

    // Eliminate type 3 partitions
    info("Eliminating type 3 partitions...");
    eliminate_type_3_partitions(mesh );

    // Double-check face colors
    if (check_face_colors(mesh) > 0)
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
                           const std::vector<double> &layer_heights,
                           const BuilderMesh &mesh)
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

  // Check mesh quality
  void check_mesh_quality(const VolumeMesh &volume_mesh, int step, bool write_to_file = false)
  {
    // Compute aspect ratios
    const auto aspect_ratios = Geometry::aspect_ratio(volume_mesh);
    const double min = std::get<0>(aspect_ratios);
    const double max = std::get<1>(aspect_ratios);
    const double median = std::get<2>(aspect_ratios);

    // Write aspect ratios to file for debugging
    const auto _aspect_rations = Geometry::aspect_ratios(volume_mesh);
    if (write_to_file)
    {
      std::ofstream file("aspect_ratios_" + str(step) + ".txt");
      for (const auto &ar : _aspect_rations)
        file << ar << std::endl;
      file.close();
    }
    // Print aspect ratios
    info("Mesh quality (aspect ratio): min = " + str(min, 3L) + ", max = " + str(max, 3L) +
         ", median = " + str(median, 3L));
  }

  // Assign face colors (closest layer height index)
  void assign_face_colors(const std::vector<double> &areas,
                          const std::vector<double> &layer_heights,
                          BuilderMesh &mesh)

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

  // Reassign face colors to avoid big jumps
  void reassign_face_colors(BuilderMesh &mesh)
  {
    for (size_t i = 0; i < mesh.ff.size(); i++)
    {
      for (const auto &j : mesh.ff[i])
      {
        const auto diff = mesh.face_colors[i] - mesh.face_colors[j];
        if (diff > 1)
        {
          mesh.face_colors[i] -= diff - 1;
        }
      }
    }
  }

  // Check face colors (big jumps)
  size_t check_face_colors(const BuilderMesh &mesh)
  {
    size_t num_big_jumps = 0;
    for (size_t i = 0; i < mesh.ff.size(); i++)
    {
      size_t _num_big_jumps = 0;
      for (const auto &j : mesh.ff[i])
      {
        const auto jump = mesh.face_colors[i] - mesh.face_colors[j];
        if (jump > 1)
          _num_big_jumps++;
      }
      if (_num_big_jumps > 0)
        num_big_jumps++;
    }

    double percentage = 100.0 * num_big_jumps / mesh.ff.size();
    info("Big jumps: " + str(num_big_jumps) + " / " + str(mesh.ff.size()) + " (" + str(percentage, 2L) +
         "%)");

    return num_big_jumps;
  }

  /// Assign vertex colors based on minimum neighbor face color
  void assign_vertex_colors(BuilderMesh &mesh)
  {
    // const size_t num_vertices = mesh.vertices.size();
    const size_t num_faces = mesh.faces.size();
    // vertex_colors.resize(num_vertices, layer_heights.size());
    for (size_t i = 0; i < num_faces; i++)
    {
      const auto &face = mesh.faces[i];
      for (const auto &j : {face.v0, face.v1, face.v2})
        mesh.vertex_colors[j] = std::min(mesh.vertex_colors[j], mesh.face_colors[i]);
    }
  }

  // Sort faces by vertex color and index
  void sort_faces_by_vertex_color_and_index(BuilderMesh &mesh)
  {
    for (auto &face : mesh.faces)
    {
      size_t c0 = mesh.vertex_colors[face.v0];
      size_t c1 = mesh.vertex_colors[face.v1];
      size_t c2 = mesh.vertex_colors[face.v2];
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
  void assign_face_partitions( BuilderMesh &mesh)
  {
    const size_t num_faces = mesh.faces.size();
    // face_partitions.resize(num_faces);
    for (size_t i = 0; i < num_faces; i++)
    {
      const auto &face = mesh.faces[i];
      const size_t c0 = mesh.vertex_colors[face.v0];
      const size_t c1 = mesh.vertex_colors[face.v1];
      const size_t c2 = mesh.vertex_colors[face.v2];
      mesh.face_partitions[i] = 3 * mesh.face_colors[i] - (c0 + c1 + c2);
    }
  }

  // Eliminate type 3 partitions which are just two stacked prisms of type 0
  void eliminate_type_3_partitions(BuilderMesh &mesh)
  {
    const size_t num_faces = mesh.face_partitions.size();
    for (size_t i = 0; i < num_faces; i++)
    {
      if (mesh.face_partitions[i] == 3)
      {
        mesh.face_colors[i]--;
        mesh.face_partitions[i] = 0;
      }
    }
  }

  void connect_column_mesh_cells(BuilderMesh &mesh, ColumnMesh &column_mesh){

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
  
  // Layer ground mesh
  VolumeMesh layer_ground_mesh(BuilderMesh &mesh, const std::vector<double> layer_heights)
  {
    // Compute max building height
    double max_building_height =layer_heights.back();
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
    _column_mesh.num_max_layers = std::ceil(max_building_height /layer_heights.back()) +1;
    _column_mesh.num_min_layers = _column_mesh.num_max_layers << (layer_heights.size() - 1);

    std::cout << "Number of layers [MAX LAYER HEIGHT]: " << _column_mesh.num_max_layers << std::endl;
    std::cout << "Number of layers [MIN LAYER HEIGHT]: " << _column_mesh.num_min_layers << std::endl;
    // Layer vertices in columns
    for (size_t j = 0; j < mesh.vertices.size(); j++)
    {
      const Vector3D &vg = mesh.vertices[j];
      const double h =layer_heights[mesh.vertex_colors[j]];
      const size_t col_size = (_column_mesh.num_min_layers >> mesh.vertex_colors[j]) + 1;
      for (size_t i = 0; i < col_size; i++)
      {
        Vector3D v(vg.x, vg.y, i * h);
        _column_mesh.vertices[j].push_back(v);
      }
      _column_mesh.vertices_offset[j + 1] = _column_mesh.vertices_offset[j] + col_size;
    }

    // Set Column Mesh Cells
    connect_column_mesh_cells(mesh, _column_mesh);
    
    // Set markers
    auto mesh_vertex_markers = mesh.face_to_vertex_markers();
    for (size_t i = 0; i < mesh.vertices.size(); i++)
    { 
      auto &markers = _column_mesh.markers[i];
      // Choose default marker: if mesh_vertex_markers[i] is >= 0 then -4 (building wall), otherwise -5 (Other)
      if (mesh_vertex_markers[i] >= 0) {
        markers.resize(_column_mesh.vertices[i].size(), -4);
      } else {
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
            _column_mesh.markers[vertex.column][vertex.index] = marker;
            
            // // Comment in the following to change marker types above buildings.
            // const size_t col_size = _column_mesh.vertices[vertex.column].size();
            // for (size_t v = vertex.index + 1; v < col_size; v++)
            // {
            //   _column_mesh.markers[vertex.column][v] = -4;
            // }
            // _column_mesh.markers[vertex.column].back() = -3;
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
    //     std::cout << _column_mesh.vertices[i][j].x <<", "<< _column_mesh.vertices[i][j].y << "," <<_column_mesh.vertices[i][j].z << "," << _column_mesh.markers[i][j] << std::endl;
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

  

  VolumeMesh layer_padding_mesh(const VolumeMesh &volume_mesh,
                                const double top_height, 
                                double max_scale = 3.0)
  {
    // Get max layer height
    std::vector<size_t> new_to_old_index;
    Mesh _top_mesh = VolumeMeshDomainPadding::extract_top_mesh(volume_mesh,new_to_old_index);
    BuilderMesh mesh(_top_mesh);
    const auto layer_heights = compute_layer_heights(mesh);

    const double max_layer_height = layer_heights.back();
    std::vector<int> vertex_matches;

    if (new_to_old_index.size() != mesh.vertices.size()){
      std::cout << "number of vertices between old and new mesh " <<new_to_old_index.size() << " | "<<mesh.vertices.size() << std::endl;
      error("Mismatch in number of vertices between old and new mesh");
    }

    ColumnMesh padding_mesh(_top_mesh);
    double top_surface_mesh_z = 0.0;
    for(auto v: mesh.vertices){
      top_surface_mesh_z = std::max(top_surface_mesh_z, v.z);
    }
    const double padding_height = top_height - top_surface_mesh_z;

    // Number of Max height layers needed to cover padding height.
    const size_t n = static_cast<size_t>(2 * padding_height / (max_layer_height * (max_scale + 1)));

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
      double H = n_k * layer_heights[color] * ( (n - 1) + s_1 * ((n - 1) * n / 2.0) );
      
      // Compute the normalization scale factor.
      double scale_norm_factor = (H > 0.0) ? ( (top_height - top_vertex_z) / H ) : 1.0;
      
      // Now, generate vertices using the scaled increments.
      double cumulative_height = 0.0;
      padding_mesh.vertices[i].push_back(mesh.vertices[i]);
      vertex_matches.push_back(new_to_old_index[i]);
      for (size_t j = 1; j < n; j++)
      {
        // const double s = 1.0 + j*(max_scale - 1)/(n-1);
        const double s = (1.0 + j * s_1)* scale_norm_factor;

        for (size_t k = 0; k < n_k; k++)
        { 
          cumulative_height += s * layer_heights[color];
          Vector3D v(mesh.vertices[i].x, mesh.vertices[i].y,
                     top_vertex_z + cumulative_height);
          padding_mesh.vertices[i].push_back(v);
          vertex_matches.push_back(-1);
        }
      }
      padding_mesh.vertices_offset[i + 1] =
      padding_mesh.vertices_offset[i] + padding_mesh.vertices[i].size();
    }

    connect_column_mesh_cells(mesh,padding_mesh);

    
    auto _padding_mesh =  padding_mesh.to_volume_mesh();

    return weld_meshes(volume_mesh, _padding_mesh, vertex_matches);
    
  }

  VolumeMesh weld_meshes(const VolumeMesh &volume_mesh, const VolumeMesh &padding_mesh,  std::vector<int> &vertex_matches)
  {
    info("Merging Volume Mesh with Padding Mesh");
    // Merge the two meshes
    VolumeMesh merged_mesh;
    const size_t volume_mesh_num_vertices = volume_mesh.vertices.size();
    merged_mesh.vertices.insert(merged_mesh.vertices.end(), volume_mesh.vertices.begin(), volume_mesh.vertices.end());
    merged_mesh.cells.insert(merged_mesh.cells.end(), volume_mesh.cells.begin(), volume_mesh.cells.end());
    // merged_mesh.markers.insert(merged_mesh.markers.end(), volume_mesh.markers.begin(), volume_mesh.markers.end());

    size_t vertices_added = 0;
    for (size_t i = 0; i < padding_mesh.vertices.size(); i++)
    {
      if (vertex_matches[i] < 0)
      {
        merged_mesh.vertices.push_back(padding_mesh.vertices[i]);
        vertex_matches[i] = volume_mesh_num_vertices + vertices_added;
        ++vertices_added;
      }
    }

    for (size_t i = 0; i < padding_mesh.cells.size(); i++)
    {
      const std::array<size_t,4> cell = {padding_mesh.cells[i].v0, padding_mesh.cells[i].v1, padding_mesh.cells[i].v2, padding_mesh.cells[i].v3};
      std::array<size_t,4> new_cell = {0,0,0,0};
      for (size_t j = 0; j < 4; j++)
      {
        new_cell[j] = vertex_matches[cell[j]];
      }
      merged_mesh.cells.push_back(Simplex3D(new_cell[0], new_cell[1], new_cell[2], new_cell[3]));
    }
    

    return merged_mesh;
  }

};

} // namespace DTCC_BUILDER

#endif
