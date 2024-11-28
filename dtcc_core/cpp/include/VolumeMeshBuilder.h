#ifndef DTCC_VOLUME_MESH_BUILDER_H
#define DTCC_VOLUME_MESH_BUILDER_H

#include <cmath>
#include <iostream>
#include <map>
#include <stack>
#include <tuple>
#include <vector>
#include <algorithm>

#include "Eigen/Eigen"
#include "Eigen/Geometry"

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

#include "model/ColumnarVolumeMesh.h"

namespace DTCC_BUILDER
{

class VolumeMeshBuilder
{
public:
  double domain_height;

  double top_height;

private:

  ColumnarVolumeMesh _col_volume_mesh;

  VolumeMesh __volume_mesh__;

  const std::vector<Surface> _buildings;

  std::vector<double> _building_ground_height;

  std::vector<double> adj_building_height;



  //double adjusted_relaxation_height;

  const GridField &_dem;

  Mesh &_ground_mesh;

  /// Number of faces in ground mesh and number of cell columns for VolumeMesh
  /// creation
  const size_t num_ground_faces;

  /// Number of vertices in ground mesh and number of vertex columns for
  /// VolumeMesh creation
  const size_t num_ground_vertices;

  /// Number of marker columns for VolumeMesh creation
  const size_t num_ground_markers;

  // Stores vertex to face mapping for ground mesh.
  std::vector<std::unordered_set<int>> vf;

  // Stores vertex to face mapping for ground mesh.
  std::vector<std::unordered_set<int>> ff;

  std::vector<int> vertex_markers;

  // Number of different layer heights used for the mesh layering.
  int num_layers;

  /// Layer heights used for the mesh layering
  std::vector<double> layer_heights;

  /// Colors (ints) are assigned to each ground mesh face describing
  /// which layer height is ideal for the triangular face.
  std::vector<int> face_colors;

  std::vector<int> vertex_colors;

  std::vector<size_t> vert_skip;

  std::vector<int> face_partition;

  std::vector<int> min_building_colors;

  std::vector<int> max_building_colors;

public:
  // Constructor
  VolumeMeshBuilder(const std::vector<Surface> &buildings,
                    const GridField &dem,
                    Mesh &ground_mesh,
                    double domain_height)
      : _buildings(buildings),
        _dem(dem),
        _ground_mesh(ground_mesh),
        domain_height(domain_height),
        num_ground_vertices(ground_mesh.vertices.size()),
        num_ground_faces(ground_mesh.faces.size()),
        num_ground_markers(ground_mesh.markers.size()),
        _col_volume_mesh(ground_mesh)
  {
    assert((num_ground_vertices > 0) &&
           "Empty ground mesh. It has no faces connecting vertices");

    assert((num_ground_faces > 0) &&
           "Empty ground mesh. It has no faces connecting faces");

    assert((num_ground_markers > 0) &&
           "Empty ground mesh. It has no face markers");

    top_height = compute_top_height();

    compute_building_ground_heights();

    face_colors.reserve(num_ground_faces);


  }

  // Destructor
  ~VolumeMeshBuilder() {}

  VolumeMesh build(const size_t smoother_iterations = 1000,
                   const double smoother_relative_tolerance = 0.001,
                   const double domain_padding_height = 0.0,
                   const size_t debug_step = 5)
  {
    info("Building Volume mesh from input ground mesh.");
    info(_ground_mesh);

    Timer t3_1("Step 3.1: Layer height computation");
    compute_layer_heights();
    t3_1.stop();
    t3_1.print();

    // Checking if compute_layer_height() returned valid layer heights
    if (num_layers == 0)
      error("Error: Empty Layer heights Vector.");


    


    domain_height = compute_relaxation_height() ;
    info("Initial domain height adjusted to max building height.. "+ str(domain_height));
    top_height = compute_top_height();
    info("Top height adjusted to max building height.. "+ str(top_height));

    Timer t3_2("Step 3.2: Ground mesh layering");
    __volume_mesh__ = layer_ground_mesh();
    t3_2.stop();
    t3_2.print();

    // info(__volume_mesh__.__str__());
    if(debug_step <= 2){
      // if (padding_height > 0.0)
      //   __volume_mesh__ = add_domain_padding(padding_height);
      return __volume_mesh__;
    }

    Timer t3_3("Step 3.3: Smooth Volume Mesh (No fixed buildings)");
    __volume_mesh__ =
        Smoother::smooth_volume_mesh(__volume_mesh__,
                                     _buildings, _dem,
                                     top_height,
                                     false,
                                     smoother_iterations,
                                     smoother_relative_tolerance);
    t3_3.stop();
    t3_3.print();

    if(debug_step == 3){
      // if (padding_height > 0.0)
      //   __volume_mesh__ = add_domain_padding(padding_height);
      return __volume_mesh__;
    }
    _col_volume_mesh._update_vertices(__volume_mesh__);

    Timer t3_4("Step 3.4: Trim volume mesh.");
    __volume_mesh__ = trim_volume_mesh();
    t3_4.stop();
    t3_4.print();




    if(debug_step == 4){
      // if (padding_height > 0.0)
      //   __volume_mesh__ = add_domain_padding(padding_height);
      return __volume_mesh__;
    }

    Timer t3_5("Step 3.5: Smooth Volume Mesh with fixed buildings");
    __volume_mesh__ = Smoother::smooth_volume_mesh(__volume_mesh__,
                                                   _buildings,
                                                  _dem,
                                                  top_height,
                                                  true,
                                                  smoother_iterations,
                                                  smoother_relative_tolerance);
    t3_5.stop();
    t3_5.print();

    _col_volume_mesh._update_vertices(__volume_mesh__);

    // if (padding_height > 0.0) __volume_mesh__ = add_domain_padding(padding_height);
    info(__volume_mesh__.__str__());

    return __volume_mesh__;
  }

private:
  /// Returns the top height of our smoothed domain.
  double compute_top_height() { return domain_height + _dem.max(); }


  /// Computes the ground heights at the centroids of buildings.
  void compute_building_ground_heights(){

    const size_t num_buildings = _buildings.size();
    _building_ground_height.resize(num_buildings);

    for(size_t i=0; i< _buildings.size(); i++){

      const Vector3D centroid = Geometry::surface_centroid(_buildings[i]);

      _building_ground_height[i] = _dem(centroid);
    }
  }

  /// Computes the relaxation height for cells above buildings.
  ///
  /// This method calculates the relaxation height, which is used to determine the height of
  /// cells above buildings in the mesh where high resolution is not necessary (bigger cells).
  double compute_relaxation_height(double buffer = 0.0)
  {

    const double max_layer_height = layer_heights.back();
    const size_t num_buildings = _buildings.size();
    std::vector<double> building_heights(num_buildings, 0.0);

    for (size_t i = 0; i < num_buildings; i++)
    {
      building_heights[i] =
          _buildings[i].max_height() - _building_ground_height[i];
    }

    double _max_building_height = 0.0;
    auto max =
        std::max_element(building_heights.begin(), building_heights.end());
    if (max != building_heights.end())
    {
      _max_building_height = *max;
    }

    // Buffer should be a non-negative float number.
    if (buffer < 0)
      buffer = 0;

    double relaxation_height = _max_building_height + buffer;
    double adj_relaxation_height =
        (std::ceil(relaxation_height / max_layer_height) + 1) *
        max_layer_height;

    info("Max building height: " + str(_max_building_height) + " m.");
    info("Relaxation height for cells above buildings: " +
         str(relaxation_height) + " m.");
    info("Adjusted relaxation height for cells above buildings: " +
         str(adj_relaxation_height) + " m.");

    return adj_relaxation_height;
  }

  /// Builds the mapping from vertices to faces.
  ///
  /// This method populates the `vf` vector, which maps each vertex in the ground mesh
  /// to the set of faces that contain that vertex.
  void build_vertex_to_face_mapping()
  {
    vf.resize(num_ground_vertices);
    for (size_t i = 0; i < num_ground_faces; i++)
    {
      vf[_ground_mesh.faces[i].v0].insert(i);
      vf[_ground_mesh.faces[i].v1].insert(i);
      vf[_ground_mesh.faces[i].v2].insert(i);
    }
  }

  /// Builds the mapping from faces to faces.
  ///
  /// This method populates the `ff` vector, which maps each face in the ground mesh
  /// to the set of neighboring faces that share at least one vertex with it.
  void build_face_to_face_mapping()
  {
    if (vf.empty())
      build_vertex_to_face_mapping();

    ff.resize(num_ground_faces);
    for (size_t i = 0; i < num_ground_faces; i++)
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
    vertex_markers.resize(num_ground_vertices, -2);

    if (!_ground_mesh.markers.size())
    {
      error("Ground mesh has no face Markers. Treating all faces as ground");
      return vertex_markers;
    }

    for (size_t f = 0; f < num_ground_faces; f++)
    {
      if (_ground_mesh.markers[f] < -2)
      {
        info("Problem problem with marker:" + str(_ground_mesh.markers[f]));
      }

      const std::array<size_t, 3> I = {_ground_mesh.faces[f].v0,
                                       _ground_mesh.faces[f].v1,
                                       _ground_mesh.faces[f].v2};

      vertex_markers[I[0]] =
          std::max(vertex_markers[I[0]], _ground_mesh.markers[f]);
      vertex_markers[I[1]] =
          std::max(vertex_markers[I[1]], _ground_mesh.markers[f]);
      vertex_markers[I[2]] =
          std::max(vertex_markers[I[2]], _ground_mesh.markers[f]);
    }

    return vertex_markers;
  }

  /// Assign face colors based on layer heights
  void assign_face_colors(std::vector<double> &areas)
  {
    size_t num_buildings = _buildings.size();
    min_building_colors.resize(num_buildings,num_layers);
    // Assign layer heights to mesh (closest by quotient)
    for (size_t i = 0; i < num_ground_faces; i++)
    {
      double h = ideal_layer_height(areas[i]);
      std::vector<double> d(layer_heights.size());
      d.reserve(layer_heights.size());

      for (size_t ih = 0; ih < layer_heights.size(); ih++)
      {
        d[ih] = std::abs(std::log(h / layer_heights[ih]));
      }
      auto min_it = std::min_element(d.begin(), d.end());
      face_colors[i] = std::distance(d.begin(), min_it);

      if (_ground_mesh.markers[i] >= 0)
        min_building_colors[_ground_mesh.markers[i]] = std::min(min_building_colors[_ground_mesh.markers[i]],face_colors[i]);


    }

    for (size_t i = 0; i < num_ground_faces; i++)
    {
      if (_ground_mesh.markers[i] >= 0 && face_colors[i] - min_building_colors[_ground_mesh.markers[i]] > 1 ) {
        face_colors[i]--; // = min_building_colors[_ground_mesh.markers[i]];
      }
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

  // Check face colors to avoid big jumps
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
    info("Big jumps: " + str(num_big_jumps) + " / " + str(ff.size()) +
         " (" + str(percentage, 2L) + "%)");

    return num_big_jumps;
  }

  /// Assign vertex colors based on minimum neighbor face color
  void assign_vertex_colors()
  {
    vertex_colors.resize(num_ground_vertices, num_layers);
    for (size_t i = 0; i < num_ground_faces; i++)
    {
      const std::vector<size_t> face = {_ground_mesh.faces[i].v0,
                                        _ground_mesh.faces[i].v1,
                                        _ground_mesh.faces[i].v2};
      for (const auto &j : face)
      {
        if (vertex_colors[j] > face_colors[i])
          vertex_colors[j] = face_colors[i];
      }
    }
  }

  /// Reassigns vertex colors to ensure that we dont have layer height differences larger than 2
  ///  between each vertex column in a building and the minimum color used in the building
  ///
  /// @note Experimental.. not currently used
  void reassign_vertex_colors()
  {
    vertex_markers = face_to_vertex_markers();

    for (size_t i = 0; i < num_ground_vertices; i++)
    {

      if (vertex_markers[i]>= 0)
      {
        const int building_index = vertex_markers[i];
        if (vertex_colors[i] - min_building_colors[building_index] >= 2)
        {
          vertex_colors[i]--;
        }
      }

    }

  }


  /// Computes ideal layer height for a regular tetrahedron
  double ideal_layer_height(double area, const double scale = 2.0)
  {
    const double c = std::pow(2.0, 1.5) * std::pow(3.0, -0.75);
    double h = scale * c * std::sqrt(area);
    return h;
  }


  // Compute layer heights for all faces in the ground mesh
  void compute_layer_heights()
  {
    std::vector<double> areas(num_ground_faces);
    for (std::size_t i = 0; i < num_ground_faces; i++)
    {
      areas[i] = Geometry::triangle_area(
          _ground_mesh.vertices[_ground_mesh.faces[i].v0],
          _ground_mesh.vertices[_ground_mesh.faces[i].v1],
          _ground_mesh.vertices[_ground_mesh.faces[i].v2]);
    }

    // Compute ideal layer heights for smallest and largest mesh sizes
    double _min_area = 0.0;
    auto min = std::min_element(areas.begin(), areas.end());
    if (min != areas.end())
    {
      _min_area = *min;
    }

    //Remove the limiting of _min_area. Used for debugging.
    // Ground Mesh has creates extremely small triangles probably due to the handling
    // of building footprints by builder.
    const double min_area_threshold = 0.2;
    if (_min_area < min_area_threshold ){
      info("Minimum area bound by lower limit of "+str(min_area_threshold) + "| min_area: " + str(_min_area));
      _min_area = min_area_threshold;
    }

    double _max_area = 0.0;
    auto max = std::max_element(areas.begin(), areas.end());
    if (max != areas.end())
    {
      _max_area = *max;
    }

    info("Min Face Area: " + str(_min_area));
    info("Max Face Area: " + str(_max_area));

    double _min_height = ideal_layer_height(_min_area);
    double _max_height = ideal_layer_height(_max_area);

    info("Min Face Ideal Layer Height: " + str(_min_height));
    info("Max Face Ideal Layer Height: " + str(_max_height));

    // Compute dyadic mesh sizes to match min/max as close as possible
    double rho = _max_height / _min_height;
    double mid = std::sqrt(_min_height * _max_height);
    num_layers = static_cast<int>(std::log2(rho) + 0.5);
    double min_height = mid / std::pow(2, num_layers / 2);

    info("rho: " + str(rho) + " | mid: " + str(mid));
    info("Number of Layers: " + str(num_layers));

    // Create layer_heights array
    layer_heights.reserve(num_layers);
    for (int i = 0; i < num_layers; i++)
    {
      layer_heights.push_back(min_height * std::pow(2.0, i));
    }

    // Print layers and their height (May Be Removed)
    for (size_t i = 0; i < layer_heights.size(); i++)
    {
      std::cout << "Layer " << i << ": " << layer_heights[i] << "m"
                << std::endl;
    }

    //Test for max builing color.
    //ideal_building_colors();

    info("Assigning layer heights to mesh");
    // Assign layer heights to mesh (closest by quotient)
    assign_face_colors(areas);

    info("Building mapping from vertices to faces");
    build_vertex_to_face_mapping();

    info("Building mapping from faces to faces");
    build_face_to_face_mapping();

    // Iteratively reassign colors to avoid big jumps
    info("Reassigning colors to avoid big jumps");
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

    // Assign vertex colors based on minimum neighbor face color
    info("Assigning vertex colors");
    assign_vertex_colors();

    //Sort vertices in each face of the ground mesh in ascending index and color.
    mesh_faces_color_sort();

    assign_face_partitions();

    type_3_partition_elimination();
  }



  /// Utility function that sorts the vertices of a face according to index and
  /// color.
  void color_sort(std::array<size_t, 3> &vertices, std::array<int, 3> &colors)
  {
    // Combine colors and vertices into a vector of pairs
    std::vector<std::pair<size_t, int>> combinedList;
    for (size_t i = 0; i < colors.size(); ++i)
    {
      combinedList.emplace_back(vertices[i], colors[i]);
    }

    // Sort the combined list first based on the first element (vertices)
    std::sort(combinedList.begin(), combinedList.end(),
              [](const auto &lhs, const auto &rhs)
              { return lhs.first < rhs.first; });

    // Sort again based on the second element (colors)
    std::sort(combinedList.begin(), combinedList.end(),
              [](const auto &lhs, const auto &rhs)
              { return lhs.second < rhs.second; });

    // Unpack the sorted list back into colors and vertices
    for (size_t i = 0; i < combinedList.size(); ++i)
    {
      colors[i] = combinedList[i].second;
      vertices[i] = combinedList[i].first;
    }
  }

  /// Sorts the vertices of the `_ground_mesh` according to index and color.
  void mesh_faces_color_sort(){
    for (size_t i = 0; i < num_ground_faces; i++)
    {
      const Simplex2D face_simplex = _ground_mesh.faces[i];

      std::array<size_t, 3> face = {face_simplex.v0, face_simplex.v1,
                                    face_simplex.v2};
      std::array<int, 3> v_colors = {vertex_colors[face_simplex.v0],
                                     vertex_colors[face_simplex.v1],
                                     vertex_colors[face_simplex.v2]};

      color_sort(face, v_colors);

      _ground_mesh.faces[i].v0 = face[0];
      _ground_mesh.faces[i].v1 = face[1];
      _ground_mesh.faces[i].v2 = face[2];
    }
  }

  /// Assigns a partition to each face of the ground_mesh.
  ///
  /// Partition types dictate how the layering of cells above each face will be done.
  /// Layering is implemented by stacking repeated structures we call prisms. The number of cells
  /// and the way the vertices are connected to create those cells is dictated by the partition type
  ///
  /// We have the following types:
  /// `Face Partition 0` : It consists of 6 vertices and 3 tetrahedral cells.
  /// `Face Partition 1` : It consists of 7 vertices and 4 tetrahedral cells.
  /// `Face Partition 2` : It consists of 8 vertices and 5 tetrahedral cells.
  /// `Face Partition 3` : It consists of 9 vertices and 6 tetrahedral cells. (Redundant partition)
  void assign_face_partitions(){
    info("Assigning Partitions to ground mesh faces.");
    face_partition.resize(num_ground_faces);
    for (size_t f = 0; f < num_ground_faces; f++)
    {

      const Simplex2D face_simplex = _ground_mesh.faces[f];

      std::array<int, 3> v_colors = {vertex_colors[face_simplex.v0],
                                     vertex_colors[face_simplex.v1],
                                     vertex_colors[face_simplex.v2]};

      face_partition[f] = 3 * face_colors[f] - (v_colors[0] + v_colors[1] + v_colors[2]);
    }
  }

  /* Type 3 Partitions are redundant.
  *  It's just two prisms with partition type 0 stacked.
  *  By eliminating these we get more granularity when handling layer heights.
  *
  *  Example:
  *  If a face is partitioned as type 3 prism, it means its color due to neighbor reassignment is c+1
  *  when all the colors of its vertices are c.
  */
  void type_3_partition_elimination(){
    info("Test: Eliminating type 3 partitions");

    size_t num_buildings = _buildings.size();
    std::vector<int> min_building_color(num_buildings,num_layers);
    for (size_t i = 0; i < _ground_mesh.faces.size(); i++)
    {
      if (face_partition[i] == 3 ){
        face_colors[i]--;
        face_partition[i] = 0;
      }
    }

    //Check to ensure there are not layer height differences larger than 1.
    check_face_colors();
  }


  VolumeMesh layer_ground_mesh()
  {
    info("Mesh Layering Function.");

    VolumeMesh volume_mesh;

    const double min_layer_height = layer_heights[0];
    const double max_layer_height = layer_heights[num_layers - 1];
    info("Number of layers: " + str(num_layers));
    info("Min Layer Height: " + str(min_layer_height));
    info("Max Layer Height: " + str(max_layer_height));

    info("Domain height:" + str(domain_height));
    const double adjusted_domain_height =
        std::ceil(domain_height / max_layer_height) * max_layer_height;
    info("Domain height adjusted to fit chosen layer heights: " +
         str(adjusted_domain_height) + "m");


    // We could Start building mesh from minimum elevation to help smoother converge faster...
    const double min_elevation = 0.0; //_dem.min();

    //size_t volume_mesh_num_vertices = 0;
    for (size_t i = 0; i < num_ground_vertices; i++)
    {
      const double layer_h = layer_heights[vertex_colors[i]];
      const size_t col_count =
          static_cast<size_t>((adjusted_domain_height / layer_h)) + 1;
      for (size_t j = 0; j < col_count; j++)
      {
        Vector3D v(_ground_mesh.vertices[i].x,
                   _ground_mesh.vertices[i].y,
                   _ground_mesh.vertices[i].z + min_elevation + j * layer_h);
        _col_volume_mesh.vertices[i].push_back(v);
        //volume_mesh_num_vertices++;
      }
      _col_volume_mesh.vertices_offset[i + 1] = _col_volume_mesh.vertices_offset[i] + _col_volume_mesh.vertices[i].size();
    }
    //col_index_offset[num_ground_faces + 1] = col_index_offset[num_ground_faces] + col_vertices[num_ground_faces - 1].size();


    //std::size_t volume_mesh_num_cells = 0;
    for (size_t i = 0; i < num_ground_faces; i++)
    {
      const Simplex2D face_simplex = _ground_mesh.faces[i];
      std::array<size_t, 3> face = {face_simplex.v0, face_simplex.v1,
                                    face_simplex.v2};
      std::array<int, 3> v_colors = {vertex_colors[face_simplex.v0],
                                     vertex_colors[face_simplex.v1],
                                     vertex_colors[face_simplex.v2]};
      //color_sort(face, v_colors);

      // const std::array<size_t, 3> column_offsets = {col_index_offset[face[0]],
      //                                               col_index_offset[face[1]],
      //                                               col_index_offset[face[2]]};

      const std::array<size_t, 3> column_offsets = {0,0,0};
      const std::array<size_t, 3> column_len = {
          _col_volume_mesh.vertices_offset[face[0] + 1] - _col_volume_mesh.vertices_offset[face[0]],
          _col_volume_mesh.vertices_offset[face[1] + 1] - _col_volume_mesh.vertices_offset[face[1]],
          _col_volume_mesh.vertices_offset[face[2] + 1] - _col_volume_mesh.vertices_offset[face[2]]};
      const size_t num_prisms =
          (column_len.back() - 1) / (1 << (face_colors[i] - v_colors[2]));


      std::vector<std::array<size_t, 3>> prism_iterator(num_prisms);

      for (size_t j = 0; j < num_prisms; j++)
      {
        prism_iterator[j] = {column_offsets[0] + j * (1 << (face_colors[i] - v_colors[0])),
                             column_offsets[1] + j * (1 << (face_colors[i] - v_colors[1])),
                             column_offsets[2] + j * (1 << (face_colors[i] - v_colors[2]))};
      }
      switch (face_partition[i])
      {
      case 0:
      {
        for (const auto &ar : prism_iterator)
        {
          const size_t k = ar[0];
          const size_t l = ar[1];
          const size_t m = ar[2];

          ColumnarVertexIndex bot_triangle_0(face[0],k),
                              bot_triangle_1(face[1],l),
                              bot_triangle_2(face[2],m);
          ColumnarVertexIndex top_triangle_0(face[0],k+1),
                              top_triangle_1(face[1],l+1),
                              top_triangle_2(face[2],m+1);


          ColumnarSimplex3D K0(bot_triangle_0, bot_triangle_1, bot_triangle_2,
                       top_triangle_2);
          ColumnarSimplex3D K1(bot_triangle_0, top_triangle_1, bot_triangle_1,
                       top_triangle_2);
          ColumnarSimplex3D K2(bot_triangle_0, top_triangle_0, top_triangle_1,
                       top_triangle_2);

          // volume_mesh.cells.emplace_back(K0);
          // volume_mesh.cells.emplace_back(K1);
          // volume_mesh.cells.emplace_back(K2);

          // col_cells[i].emplace_back(K0);
          // col_cells[i].emplace_back(K1);
          // col_cells[i].emplace_back(K2);

          _col_volume_mesh.cells[i].emplace_back(K0);
          _col_volume_mesh.cells[i].emplace_back(K1);
          _col_volume_mesh.cells[i].emplace_back(K2);
          //volume_mesh_num_cells += 3;
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


          ColumnarVertexIndex bot_triangle_0(face[0],k),
                              bot_triangle_1(face[1],l),
                              bot_triangle_2(face[2],m);
          ColumnarVertexIndex mid_triangle_0(face[0],k+1);
          ColumnarVertexIndex top_triangle_0(face[0],k+2),
                              top_triangle_1(face[1],l+1),
                              top_triangle_2(face[2],m+1);


          // size_t bot_triangle_0 = k, bot_triangle_1 = l, bot_triangle_2 = m;
          // size_t mid_triangle_0 = k + 1;
          // size_t top_triangle_0 = k + 2, top_triangle_1 = l + 1,
          //        top_triangle_2 = m + 1;

          ColumnarSimplex3D K0(bot_triangle_0, bot_triangle_1, bot_triangle_2,
                       mid_triangle_0);
          ColumnarSimplex3D K1(bot_triangle_1, top_triangle_2, bot_triangle_2,
                       mid_triangle_0);
          ColumnarSimplex3D K2(bot_triangle_1, top_triangle_2, mid_triangle_0,
                       top_triangle_1);
          ColumnarSimplex3D K3(mid_triangle_0, top_triangle_0, top_triangle_1,
                       top_triangle_2);

          // volume_mesh.cells.emplace_back(K0);
          // volume_mesh.cells.emplace_back(K1);
          // volume_mesh.cells.emplace_back(K2);
          // volume_mesh.cells.emplace_back(K3);

          _col_volume_mesh.cells[i].emplace_back(K0);
          _col_volume_mesh.cells[i].emplace_back(K1);
          _col_volume_mesh.cells[i].emplace_back(K2);
          _col_volume_mesh.cells[i].emplace_back(K3);
          //volume_mesh_num_cells += 4;
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

          // size_t bot_triangle_0 = k, bot_triangle_1 = l, bot_triangle_2 = m;
          // size_t mid_triangle_0 = k + 1, mid_triangle_1 = l + 1;
          // size_t top_triangle_0 = k + 2, top_triangle_1 = l + 2,
          //        top_triangle_2 = m + 1;

          ColumnarVertexIndex bot_triangle_0(face[0],k),
                              bot_triangle_1(face[1],l),
                              bot_triangle_2(face[2],m);
          ColumnarVertexIndex mid_triangle_0(face[0],k+1),
                              mid_triangle_1(face[1],l+1);
          ColumnarVertexIndex top_triangle_0(face[0],k+2),
                              top_triangle_1(face[1],l+2),
                              top_triangle_2(face[2],m+1);

          ColumnarSimplex3D K0(bot_triangle_0, bot_triangle_1, bot_triangle_2,
                       mid_triangle_1);
          ColumnarSimplex3D K1(bot_triangle_0, bot_triangle_2, mid_triangle_0,
                       mid_triangle_1);
          ColumnarSimplex3D K2(bot_triangle_2, top_triangle_2, mid_triangle_0,
                       mid_triangle_1);
          ColumnarSimplex3D K3(top_triangle_0, top_triangle_2, top_triangle_1,
                       mid_triangle_0);
          ColumnarSimplex3D K4(top_triangle_1, top_triangle_2, mid_triangle_1,
                       mid_triangle_0);

          // volume_mesh.cells.emplace_back(K0);
          // volume_mesh.cells.emplace_back(K1);
          // volume_mesh.cells.emplace_back(K2);
          // volume_mesh.cells.emplace_back(K3);
          // volume_mesh.cells.emplace_back(K4);

          _col_volume_mesh.cells[i].emplace_back(K0);
          _col_volume_mesh.cells[i].emplace_back(K1);
          _col_volume_mesh.cells[i].emplace_back(K2);
          _col_volume_mesh.cells[i].emplace_back(K3);
          _col_volume_mesh.cells[i].emplace_back(K4);
          //volume_mesh_num_cells += 5;
        }
      }
      break;

      case 3:
      {
        info("Layer stage:" + str(i) + "Partition type 3");
        for (const auto &ar : prism_iterator)
        {
          const size_t k = ar[0];
          const size_t l = ar[1];
          const size_t m = ar[2];

          // size_t bot_triangle_0 = k, bot_triangle_1 = l, bot_triangle_2 = m;
          // size_t mid_triangle_0 = k + 1, mid_triangle_1 = l + 1,
          //        mid_triangle_2 = m + 1;
          // size_t top_triangle_0 = k + 2, top_triangle_1 = l + 2,
          //        top_triangle_2 = m + 2;

          ColumnarVertexIndex bot_triangle_0(face[0],k),
                              bot_triangle_1(face[1],l),
                              bot_triangle_2(face[2],m);
          ColumnarVertexIndex mid_triangle_0(face[0],k+1),
                              mid_triangle_1(face[1],l+1),
                              mid_triangle_2(face[2],m+1);  ;
          ColumnarVertexIndex top_triangle_0(face[0],k+2),
                              top_triangle_1(face[1],l+2),
                              top_triangle_2(face[2],m+2);

          ColumnarSimplex3D K0(bot_triangle_0, bot_triangle_1, bot_triangle_2,
                       mid_triangle_2);
          ColumnarSimplex3D K1(bot_triangle_0, mid_triangle_1, bot_triangle_1,
                       mid_triangle_2);
          ColumnarSimplex3D K2(bot_triangle_0, mid_triangle_0, mid_triangle_1,
                       mid_triangle_2);
          ColumnarSimplex3D K3(mid_triangle_0, mid_triangle_1, mid_triangle_2,
                       top_triangle_2);
          ColumnarSimplex3D K4(mid_triangle_0, top_triangle_1, mid_triangle_1,
                       top_triangle_2);
          ColumnarSimplex3D K5(mid_triangle_0, top_triangle_0, top_triangle_1,
                       top_triangle_2);

          // volume_mesh.cells.emplace_back(K0);
          // volume_mesh.cells.emplace_back(K1);
          // volume_mesh.cells.emplace_back(K2);
          // volume_mesh.cells.emplace_back(K3);
          // volume_mesh.cells.emplace_back(K4);
          // volume_mesh.cells.emplace_back(K5);

          _col_volume_mesh.cells[i].emplace_back(K0);
          _col_volume_mesh.cells[i].emplace_back(K1);
          _col_volume_mesh.cells[i].emplace_back(K2);
          _col_volume_mesh.cells[i].emplace_back(K3);
          _col_volume_mesh.cells[i].emplace_back(K4);
          _col_volume_mesh.cells[i].emplace_back(K5);
          //volume_mesh_num_cells += 6;
        }
        break;
      }
      default:
        error("Face Coloring Error: Large layer height difference: " +
              str(face_partition[i]) +
              "\nFace:" + str(i) + " color: " + str(face_colors[i])+ str(min_building_colors[_ground_mesh.markers[i]])+
              "\n v0 color: " + str(v_colors[0]) +
              "\n v1 color: " + str(v_colors[1]) +
              "\n v2 color: " + str(v_colors[2]));
        break;
      }
    }



    // Add Markers
    //col_markers.resize(num_ground_vertices);
    auto mesh_vertex_markers = face_to_vertex_markers();
    for (size_t i = 0; i < num_ground_vertices; i++)
    {
      std::vector<int> tmp(_col_volume_mesh.vertices[i].size(), -4);
      if (mesh_vertex_markers[i] == -1 || mesh_vertex_markers[i] == -2)
      {
        tmp.front() = mesh_vertex_markers[i];
      }
      else
      {
        tmp.front() = -1;
      }
      tmp.back() = -3;
      _col_volume_mesh.markers[i] = tmp;
    }

    return _col_volume_mesh.to_volume_mesh();
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
    vertex_markers.resize(num_ground_vertices, std::numeric_limits<int>::max());

    if (!_ground_mesh.markers.size())
    {
      error("Ground mesh has no face Markers. Treating all faces as ground");
      return vertex_markers;
    }

    for (size_t f = 0; f < num_ground_faces; f++)
    {
      if (_ground_mesh.markers[f] < -2)
      {
        info("Problem problem with marker:" + str(_ground_mesh.markers[f]));
      }

      const std::array<size_t, 3> I = {_ground_mesh.faces[f].v0,
                                       _ground_mesh.faces[f].v1,
                                       _ground_mesh.faces[f].v2};

      vertex_markers[I[0]] =
          std::min(vertex_markers[I[0]], _ground_mesh.markers[f]);
      vertex_markers[I[1]] =
          std::min(vertex_markers[I[1]], _ground_mesh.markers[f]);
      vertex_markers[I[2]] =
          std::min(vertex_markers[I[2]], _ground_mesh.markers[f]);
    }

    return vertex_markers;
  }

  /// Adjust the input building heights to the maximum layer height used by
  /// the vertices that belong to the building.
  std::vector<double> adjust_building_heights()
  {
    const size_t num_buildings = _buildings.size();
    max_building_colors.resize(num_buildings,0);
    // build map from buildings to cells in 2D mesh

    info("Building map from buildings to cells in 2D mesh");
    std::vector<std::vector<size_t>> building_faces(num_buildings);
    for (size_t face_index = 0; face_index < num_ground_faces; face_index++)
    {
      const int building_index = _ground_mesh.markers[face_index];
      if (building_index >= 0)
      {
        building_faces[building_index].push_back(face_index);
      }
    }

    std::vector<int> min_building_colors(num_buildings, 0);
    for (size_t building_index = 0; building_index < num_buildings;
         building_index++)
    {
      // Note: We could either work with face or vertex colors/layer-heights
      for (const size_t &f : building_faces[building_index])
      {
       // Working with face colors... This choice also changes the calculation for the number of prisms created.
       max_building_colors[building_index] =
           std::max({face_colors[f], max_building_colors[building_index]});

      //   const Simplex2D face = _ground_mesh.faces[f];
      //   const int max_v_color = std::max({vertex_colors[face.v0],
      //                                     vertex_colors[face.v1],
      //                                     vertex_colors[face.v2]});
      //   building_colors[building_index] =
      //       std::max({max_v_color, building_colors[building_index]});
      }
    }


    // Loop over buildings
    std::vector<double> adj_building_heights(num_buildings, 0);
    for (size_t building_index = 0; building_index < num_buildings;
         building_index++)
    {
      const double fitting_layer_height =
          layer_heights[max_building_colors[building_index]];

      double n = (_buildings[building_index].max_height() - _building_ground_height[building_index])/ fitting_layer_height; //CHANGE ASAP WROOONG

      adj_building_heights[building_index] = fitting_layer_height;
      if (n>1)
      {
        adj_building_heights[building_index] *= std::floor(n);
      }

    }

    std::cout << "Adjusted Building Heights Size: " << adj_building_heights.size() << std::endl;
    return adj_building_heights;
  }

  //   void ideal_building_colors()
  // {
  //   const size_t num_buildings = _buildings.size();

  //   building_colors.resize(num_buildings,0);
  //   for (size_t i = 0; i < num_buildings; i++)
  //   {
  //     const double building_height = _buildings[i].max_height() - _building_ground_height[i];


  //     double min_excess = 1.0;
  //     for (size_t j = layer_heights.size(); j < 0; --j)
  //     {
  //       double excess = fmod(building_height,layer_heights[j]);

  //       if (excess < min_excess)
  //       {
  //         building_colors[i] = j;
  //       }

  //     }
  //     // auto max_element_it = std::max_element(max_color_index.rbegin(),max_color_index.rend());
  //     // building_colors[i] = std::distance(max_color_index.begin(),max_element_it.base()-1);

  //     //std::cout << " Ideal color: " << building_colors[i] << " bh: "<< building_height << " lh: "<< layer_heights[building_colors[i]]<<  std::endl;

  //   }

  // }


  // TO BE REMOVED:
  // This function is not used.
  // It returns the starting index of the vertex that should be included in the trimmed column.
  std::vector<size_t> _get_building_trimming_indexes()
  {
    const size_t num_buildings = _buildings.size();
    std::vector<size_t> trimming_index(num_buildings,0);
    std::vector<size_t> trimming_index_j(num_buildings,0);

    std::vector<int> min_vertex_markers = face_to_vertex_markers_min();
    std::vector<double> max_layer_height(num_buildings,0.0);


    for (size_t i = 0; i < num_buildings; i++)
    {
      const double building_height = _buildings[i].max_height() - _building_ground_height[i];
      if (max_layer_height[i] > 0.0)
        trimming_index[i] = static_cast<size_t>(building_height/ max_layer_height[i]); //*( 1 << max_building_colors[i]);
      else
       trimming_index[i] = 1; //<< max_building_colors[i];

      // if (max_layer_height[i] > 0 )
      //   std::cout << i << ") Building h: " << building_height
      //             << " max_cell_height: " << max_layer_height[i]
      //             << " j: "<< trimming_index_j[i]
      //             << " trimming_index: " <<trimming_index[i]
      //             << " max layer h: " << layer_heights[ max_building_colors[i]]
      //             << " min layer h: " << layer_heights[ min_building_colors[i]]
      //             << std::endl;

    }

    return trimming_index;

  }

  std::vector<size_t> get_building_trimming_index(double buffer = 0.0)
  {
    const size_t num_buildings = _buildings.size();
    std::vector<size_t> trimming_index(num_buildings,0);

    //Looping over all Vertex columns of the columnar volume mesh.
    for (size_t i = 0; i < num_ground_vertices; i++)
    {
      int marker = vertex_markers[i];
      if (marker >=0 && vertex_colors[i] == max_building_colors[marker] )
      {
        size_t j = 0;
        while (_col_volume_mesh.vertices[i][j + 1].z <= (_buildings[marker].max_height() + buffer))
        {
          j++;
        }
        trimming_index[marker] = std::max(trimming_index[marker],static_cast<size_t>(j * (1<<vertex_colors[i])));

        //For debugging.
        //std::cout <<i << ") Building: "<<  marker << " trimming_index: " << j << " v_color: " << vertex_colors[i]<< std::endl;
      }

    }

    for (size_t i = 0; i < trimming_index.size(); i++)
    {
        if (trimming_index[i] == std::numeric_limits<size_t>::max())
        {
          trimming_index[i] = 1<<max_building_colors[i];
        }
        //trimming_index[i] *= (1<<max_building_colors[i]);
        // std::cout<<i << ") Trimming_index: " << trimming_index[i] << std::endl;
    }
    return trimming_index;
  }

  std::vector<double> get_adjusted_building_heights()
  {
    const size_t num_buildings = _buildings.size();
    auto trimming_index = get_building_trimming_index();

    std::vector<double> adj_b_heights(num_buildings, 0.0);
    for (size_t i = 0; i < num_buildings; i++)
    {
      const double building_height =
          _buildings[i].max_height() - _building_ground_height[i];
      adj_b_heights[i] = trimming_index[i] * layer_heights[max_building_colors[i]];

      if ((adj_b_heights[i] - building_height) > layer_heights[max_building_colors[i]])
        adj_b_heights[i] -= layer_heights[max_building_colors[i]];
      else if ((building_height - adj_b_heights[i]) >
               layer_heights[max_building_colors[i]])
        adj_b_heights[i] += layer_heights[max_building_colors[i]];
    }

    return adj_b_heights;
  }

  VolumeMesh trim_volume_mesh()
  {
    info("Trimming Volume Mesh..");


    // Adjust Building heights to the largest layer height assigned to buildings
    // faces.
    std::vector<double> adj_building_heights = adjust_building_heights();
    std::vector<int> min_vertex_markers = face_to_vertex_markers_min();
    std::cout<< "Short test to find error. Ajusted buisdfsdfflding Heights: " << std::endl;

    //adj_building_heights = get_adjusted_building_heights();

    auto trimming_index = get_building_trimming_index();
    // In this vector we store how many minimum layer vertices in each column
    // are skipped to trim mesh and create a building
    vert_skip.resize(num_ground_vertices, 0);
    //size_t volume_mesh_num_vertices = 0;
    for (size_t i = 0; i < num_ground_vertices; i++)
    {
      if (min_vertex_markers[i] >=0 )
      {

        const double building_adj_height =
            adj_building_heights[min_vertex_markers[i]];
        const double layer_h = layer_heights[vertex_colors[i]];


        size_t col_begin = static_cast<size_t>((building_adj_height / layer_h));
        vert_skip[i] = col_begin;

        // col_begin = (trimming_index[min_vertex_markers[i]] )/ (1 << vertex_colors[i]);
        //std::cout << "Vcol: " << i << " building_adjusted_height: " << building_adj_height << " fitting layer: " << layer_h << "col_skip: "<< col_begin << std::endl;
        std::vector<Vector3D> trimmed_column;
        const size_t col_end =  _col_volume_mesh.vertices[i].size();
        for (size_t j = col_begin ; j < col_end; j++)
        {
          Vector3D v = _col_volume_mesh.vertices[i][j];
          trimmed_column.push_back(v);
        }
        _col_volume_mesh.vertices[i] = trimmed_column;
      }
      _col_volume_mesh.vertices_offset[i + 1] = _col_volume_mesh.vertices_offset[i] + _col_volume_mesh.vertices[i].size();
      //volume_mesh_num_vertices += col_vertices[i].size();
    }

    std::cout<< "Short test to find error. Ajusted building Heights: " << adj_building_heights.size()<< std::endl;

    //std::vector<int> face_partition(num_ground_faces,-1);

    for (size_t i = 0; i < num_ground_faces; i++)
    {
      // No trimming needed for cell columns that are not marked as buildings.
      if (_ground_mesh.markers[i] <0 )
      {
       continue;
      }

      _col_volume_mesh.cells[i].clear();
      const Simplex2D face_simplex = _ground_mesh.faces[i];
      const int face_marker = _ground_mesh.markers[i];
      std::array<size_t, 3> face = {face_simplex.v0, face_simplex.v1,
                                    face_simplex.v2};

      std::array<int, 3> v_colors = {vertex_colors[face_simplex.v0],
                                     vertex_colors[face_simplex.v1],
                                     vertex_colors[face_simplex.v2]};

      const std::array<size_t, 3> column_offsets = {0,0,0};

      const std::array<size_t, 3> col_skip = {
          vert_skip[face[0]] * (1 << v_colors[0]),
          vert_skip[face[1]] * (1 << v_colors[1]),
          vert_skip[face[2]] * (1 << v_colors[2])};



      const std::array<size_t, 3> column_len = {
          _col_volume_mesh.vertices_offset[face[0] + 1] - _col_volume_mesh.vertices_offset[face[0]],
          _col_volume_mesh.vertices_offset[face[1] + 1] - _col_volume_mesh.vertices_offset[face[1]],
          _col_volume_mesh.vertices_offset[face[2] + 1] - _col_volume_mesh.vertices_offset[face[2]]};


      const std::array<size_t, 3> column_len_nrml = {
          (column_len[0] - 1) * (1 << v_colors[0]),
          (column_len[1] - 1) * (1 << v_colors[1]),
          (column_len[2] - 1) * (1 << v_colors[2])};

      size_t min_nrml_col_len = std::min(
          {column_len_nrml[0], column_len_nrml[1], column_len_nrml[2]});

      size_t max_col_skip = std::max({col_skip[0], col_skip[1], col_skip[2]});

      size_t face_col_skip = static_cast<size_t>(
          (adj_building_heights[face_marker] / layer_heights[0]));

      // size_t face_col_skip = static_cast<size_t>(
      //     ((_buildings[face_marker].max_height() - _building_ground_height[face_marker]) / layer_heights[0]));




      //We do that for the corner of the buildings. Although no vertices are trimmed we should not create cells because we elimenate the corners of the buildings.
      min_nrml_col_len -= face_col_skip - max_col_skip;
      max_col_skip = face_col_skip; //(trimming_index[face_marker])* 1 <<(max_building_colors[face_marker]);


      const std::array<size_t, 3> prism_start_indexes{
          static_cast<size_t>((max_col_skip - col_skip[0]) /
                              (1 << v_colors[0])),
          static_cast<size_t>((max_col_skip - col_skip[1]) /
                              (1 << v_colors[1])),
          static_cast<size_t>((max_col_skip - col_skip[2]) /
                              (1 << v_colors[2]))};

      const size_t num_prisms =
          static_cast<size_t>(min_nrml_col_len / (1 << face_colors[i]));

      // const size_t num_prisms =
      //     static_cast<size_t>(min_nrml_col_len / (1 << v_colors[2]));

      std::vector<std::array<size_t, 3>> prism_iterator(num_prisms);

      for (size_t j = 0; j < num_prisms; j++)
      {
        prism_iterator[j] = {column_offsets[0] + prism_start_indexes[0] +
                                 (j << (face_colors[i] - v_colors[0])),
                             column_offsets[1] + prism_start_indexes[1] +
                                 (j << (face_colors[i] - v_colors[1])),
                             column_offsets[2] + prism_start_indexes[2] +
                                 (j << (face_colors[i] - v_colors[2]))};

      }

      switch (face_partition[i])
      {
      case 0:
      {
        for (const auto &ar : prism_iterator)
        {
          const size_t k = ar[0];
          const size_t l = ar[1];
          const size_t m = ar[2];

          ColumnarVertexIndex bot_triangle_0(face[0],k),
                              bot_triangle_1(face[1],l),
                              bot_triangle_2(face[2],m);
          ColumnarVertexIndex top_triangle_0(face[0],k+1),
                              top_triangle_1(face[1],l+1),
                              top_triangle_2(face[2],m+1);


          ColumnarSimplex3D K0(bot_triangle_0, bot_triangle_1, bot_triangle_2,
                       top_triangle_2);
          ColumnarSimplex3D K1(bot_triangle_0, top_triangle_1, bot_triangle_1,
                       top_triangle_2);
          ColumnarSimplex3D K2(bot_triangle_0, top_triangle_0, top_triangle_1,
                       top_triangle_2);

          // volume_mesh.cells.emplace_back(K0);
          // volume_mesh.cells.emplace_back(K1);
          // volume_mesh.cells.emplace_back(K2);

          // col_cells[i].emplace_back(K0);
          // col_cells[i].emplace_back(K1);
          // col_cells[i].emplace_back(K2);

          _col_volume_mesh.cells[i].emplace_back(K0);
          _col_volume_mesh.cells[i].emplace_back(K1);
          _col_volume_mesh.cells[i].emplace_back(K2);
          //volume_mesh_num_cells += 3;
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


          ColumnarVertexIndex bot_triangle_0(face[0],k),
                              bot_triangle_1(face[1],l),
                              bot_triangle_2(face[2],m);
          ColumnarVertexIndex mid_triangle_0(face[0],k+1);
          ColumnarVertexIndex top_triangle_0(face[0],k+2),
                              top_triangle_1(face[1],l+1),
                              top_triangle_2(face[2],m+1);


          // size_t bot_triangle_0 = k, bot_triangle_1 = l, bot_triangle_2 = m;
          // size_t mid_triangle_0 = k + 1;
          // size_t top_triangle_0 = k + 2, top_triangle_1 = l + 1,
          //        top_triangle_2 = m + 1;

          ColumnarSimplex3D K0(bot_triangle_0, bot_triangle_1, bot_triangle_2,
                       mid_triangle_0);
          ColumnarSimplex3D K1(bot_triangle_1, top_triangle_2, bot_triangle_2,
                       mid_triangle_0);
          ColumnarSimplex3D K2(bot_triangle_1, top_triangle_2, mid_triangle_0,
                       top_triangle_1);
          ColumnarSimplex3D K3(mid_triangle_0, top_triangle_0, top_triangle_1,
                       top_triangle_2);

          // volume_mesh.cells.emplace_back(K0);
          // volume_mesh.cells.emplace_back(K1);
          // volume_mesh.cells.emplace_back(K2);
          // volume_mesh.cells.emplace_back(K3);

          _col_volume_mesh.cells[i].emplace_back(K0);
          _col_volume_mesh.cells[i].emplace_back(K1);
          _col_volume_mesh.cells[i].emplace_back(K2);
          _col_volume_mesh.cells[i].emplace_back(K3);
          //volume_mesh_num_cells += 4;
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

          // size_t bot_triangle_0 = k, bot_triangle_1 = l, bot_triangle_2 = m;
          // size_t mid_triangle_0 = k + 1, mid_triangle_1 = l + 1;
          // size_t top_triangle_0 = k + 2, top_triangle_1 = l + 2,
          //        top_triangle_2 = m + 1;

          ColumnarVertexIndex bot_triangle_0(face[0],k),
                              bot_triangle_1(face[1],l),
                              bot_triangle_2(face[2],m);
          ColumnarVertexIndex mid_triangle_0(face[0],k+1),
                              mid_triangle_1(face[1],l+1);
          ColumnarVertexIndex top_triangle_0(face[0],k+2),
                              top_triangle_1(face[1],l+2),
                              top_triangle_2(face[2],m+1);

          ColumnarSimplex3D K0(bot_triangle_0, bot_triangle_1, bot_triangle_2,
                       mid_triangle_1);
          ColumnarSimplex3D K1(bot_triangle_0, bot_triangle_2, mid_triangle_0,
                       mid_triangle_1);
          ColumnarSimplex3D K2(bot_triangle_2, top_triangle_2, mid_triangle_0,
                       mid_triangle_1);
          ColumnarSimplex3D K3(top_triangle_0, top_triangle_2, top_triangle_1,
                       mid_triangle_0);
          ColumnarSimplex3D K4(top_triangle_1, top_triangle_2, mid_triangle_1,
                       mid_triangle_0);

          // volume_mesh.cells.emplace_back(K0);
          // volume_mesh.cells.emplace_back(K1);
          // volume_mesh.cells.emplace_back(K2);
          // volume_mesh.cells.emplace_back(K3);
          // volume_mesh.cells.emplace_back(K4);

          _col_volume_mesh.cells[i].emplace_back(K0);
          _col_volume_mesh.cells[i].emplace_back(K1);
          _col_volume_mesh.cells[i].emplace_back(K2);
          _col_volume_mesh.cells[i].emplace_back(K3);
          _col_volume_mesh.cells[i].emplace_back(K4);
          //volume_mesh_num_cells += 5;
        }
      }
      break;

      case 3:
      {
        info("Trim stage:" + str(i) + "Partition type 3");
        for (const auto &ar : prism_iterator)
        {
          const size_t k = ar[0];
          const size_t l = ar[1];
          const size_t m = ar[2];

          // size_t bot_triangle_0 = k, bot_triangle_1 = l, bot_triangle_2 = m;
          // size_t mid_triangle_0 = k + 1, mid_triangle_1 = l + 1,
          //        mid_triangle_2 = m + 1;
          // size_t top_triangle_0 = k + 2, top_triangle_1 = l + 2,
          //        top_triangle_2 = m + 2;

          ColumnarVertexIndex bot_triangle_0(face[0],k),
                              bot_triangle_1(face[1],l),
                              bot_triangle_2(face[2],m);
          ColumnarVertexIndex mid_triangle_0(face[0],k+1),
                              mid_triangle_1(face[1],l+1),
                              mid_triangle_2(face[2],m+1);  ;
          ColumnarVertexIndex top_triangle_0(face[0],k+2),
                              top_triangle_1(face[1],l+2),
                              top_triangle_2(face[2],m+2);

          ColumnarSimplex3D K0(bot_triangle_0, bot_triangle_1, bot_triangle_2,
                       mid_triangle_2);
          ColumnarSimplex3D K1(bot_triangle_0, mid_triangle_1, bot_triangle_1,
                       mid_triangle_2);
          ColumnarSimplex3D K2(bot_triangle_0, mid_triangle_0, mid_triangle_1,
                       mid_triangle_2);
          ColumnarSimplex3D K3(mid_triangle_0, mid_triangle_1, mid_triangle_2,
                       top_triangle_2);
          ColumnarSimplex3D K4(mid_triangle_0, top_triangle_1, mid_triangle_1,
                       top_triangle_2);
          ColumnarSimplex3D K5(mid_triangle_0, top_triangle_0, top_triangle_1,
                       top_triangle_2);

          // volume_mesh.cells.emplace_back(K0);
          // volume_mesh.cells.emplace_back(K1);
          // volume_mesh.cells.emplace_back(K2);
          // volume_mesh.cells.emplace_back(K3);
          // volume_mesh.cells.emplace_back(K4);
          // volume_mesh.cells.emplace_back(K5);

          _col_volume_mesh.cells[i].emplace_back(K0);
          _col_volume_mesh.cells[i].emplace_back(K1);
          _col_volume_mesh.cells[i].emplace_back(K2);
          _col_volume_mesh.cells[i].emplace_back(K3);
          _col_volume_mesh.cells[i].emplace_back(K4);
          _col_volume_mesh.cells[i].emplace_back(K5);
          //volume_mesh_num_cells += 6;
        }
        break;
      }
      default:
        error("Face Color Error: Large layer height difference");
        break;
      }
    }

    const auto mesh_vertex_markers = face_to_vertex_markers();

    // col_markers.resize(num_ground_vertices);
    for (size_t i = 0; i < num_ground_vertices; i++)
    {
      const int marker = mesh_vertex_markers[i];
      std::vector<int> tmp(_col_volume_mesh.vertices[i].size(), -4);
      if (marker < 0)
      {
        tmp.front() = marker;
      }
      else
      {
        const double layer_h = layer_heights[vertex_colors[i]];
        const size_t idx =
            static_cast<size_t>((adj_building_heights[marker] / layer_h)) -
            vert_skip[i];
        if (idx == 0)
        {
          tmp.front() = marker;
        }
        else
        {
          tmp.front() = -1;
          tmp[idx] = marker;
          //for (size_t j = 1; j < idx; j++) tmp[j] = - 5 - marker;

        }
        tmp.back() = -3;
      }
      _col_volume_mesh.markers[i] = tmp;
    }

    return _col_volume_mesh.to_volume_mesh();
  }


VolumeMesh add_domain_padding(double padding_height, double max_scale = 3.0)
{
  //Max Layer height.
  const double max_layer_height = layer_heights[num_layers - 1];

  //Number of Max height layers needed to cover padding height.
  const size_t n = static_cast<size_t>(2 * padding_height/(max_layer_height*(max_scale + 1)));


  info("Padding Domain with "+ str(n) + " max layers scaled from 1 to "+ str(max_scale));
  std::vector<size_t> offset_before_padding = _col_volume_mesh.vertices_offset;

  for (size_t i = 0; i < num_ground_vertices; i++)
  {
    //Z coordinate of the last veertex of each column.
    const double top_vertex_z = _col_volume_mesh.vertices[i].back().z ;

    const size_t n_k = static_cast<size_t>(max_layer_height / layer_heights[vertex_colors[i]]);
    double layer_h = 0.0;
    const double s_1 = n>1 ? (max_scale - 1)/(n-1): 0 ;
    for (size_t j = 1; j < n; j++)
      {
        // const double s = 1.0 + j*(max_scale - 1)/(n-1);
        const double s = 1.0 + j*s_1;

        for (size_t k = 0; k < n_k; k++)
        {
          layer_h += s * layer_heights[vertex_colors[i]];
          Vector3D v(_ground_mesh.vertices[i].x,
                    _ground_mesh.vertices[i].y,
                    top_vertex_z + layer_h);
          _col_volume_mesh.vertices[i].push_back(v);
        }
      }
    _col_volume_mesh.vertices_offset[i+1] =  _col_volume_mesh.vertices_offset[i] +  _col_volume_mesh.vertices[i].size();
  }

  size_t volume_mesh_num_cells = 0;
  for (size_t i = 0; i < num_ground_faces; i++)
    {
      const Simplex2D face_simplex = _ground_mesh.faces[i];
      std::array<size_t, 3> face = {face_simplex.v0, face_simplex.v1,
                                    face_simplex.v2};
      std::array<int, 3> v_colors = {vertex_colors[face_simplex.v0],
                                     vertex_colors[face_simplex.v1],
                                     vertex_colors[face_simplex.v2]};
      //color_sort(face, v_colors);

      // const std::array<size_t, 3> column_offsets = {col_index_offset[face[0]],
      //                                               col_index_offset[face[1]],
      //                                               col_index_offset[face[2]]};

      const std::array<size_t, 3> column_offsets = {
                offset_before_padding[face[0]+1]- offset_before_padding[face[0]] - 1,
                offset_before_padding[face[1]+1]- offset_before_padding[face[1]] - 1,
                offset_before_padding[face[2]+1]- offset_before_padding[face[2]] - 1,
      };
      const std::array<size_t, 3> column_len = {
          _col_volume_mesh.vertices_offset[face[0] + 1] - _col_volume_mesh.vertices_offset[face[0]],
          _col_volume_mesh.vertices_offset[face[1] + 1] - _col_volume_mesh.vertices_offset[face[1]],
          _col_volume_mesh.vertices_offset[face[2] + 1] - _col_volume_mesh.vertices_offset[face[2]]};
      const size_t num_prisms =
          (column_len.back() - column_offsets.back() - 1) / (1 << (face_colors[i] - v_colors[2]));




      std::vector<std::array<size_t, 3>> prism_iterator(num_prisms);

      for (size_t j = 0; j < num_prisms; j++)
      {
        prism_iterator[j] = {column_offsets[0] + j * (1 << (face_colors[i] - v_colors[0])),
                             column_offsets[1] + j * (1 << (face_colors[i] - v_colors[1])),
                             column_offsets[2] + j * (1 << (face_colors[i] - v_colors[2]))};
      }
      switch (face_partition[i])
      {
      case 0:
      {
        for (const auto &ar : prism_iterator)
        {
          const size_t k = ar[0];
          const size_t l = ar[1];
          const size_t m = ar[2];

          ColumnarVertexIndex bot_triangle_0(face[0],k),
                              bot_triangle_1(face[1],l),
                              bot_triangle_2(face[2],m);
          ColumnarVertexIndex top_triangle_0(face[0],k+1),
                              top_triangle_1(face[1],l+1),
                              top_triangle_2(face[2],m+1);


          ColumnarSimplex3D K0(bot_triangle_0, bot_triangle_1, bot_triangle_2,
                       top_triangle_2);
          ColumnarSimplex3D K1(bot_triangle_0, top_triangle_1, bot_triangle_1,
                       top_triangle_2);
          ColumnarSimplex3D K2(bot_triangle_0, top_triangle_0, top_triangle_1,
                       top_triangle_2);

          // volume_mesh.cells.emplace_back(K0);
          // volume_mesh.cells.emplace_back(K1);
          // volume_mesh.cells.emplace_back(K2);

          // col_cells[i].emplace_back(K0);
          // col_cells[i].emplace_back(K1);
          // col_cells[i].emplace_back(K2);

          _col_volume_mesh.cells[i].emplace_back(K0);
          _col_volume_mesh.cells[i].emplace_back(K1);
          _col_volume_mesh.cells[i].emplace_back(K2);
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


          ColumnarVertexIndex bot_triangle_0(face[0],k),
                              bot_triangle_1(face[1],l),
                              bot_triangle_2(face[2],m);
          ColumnarVertexIndex mid_triangle_0(face[0],k+1);
          ColumnarVertexIndex top_triangle_0(face[0],k+2),
                              top_triangle_1(face[1],l+1),
                              top_triangle_2(face[2],m+1);


          // size_t bot_triangle_0 = k, bot_triangle_1 = l, bot_triangle_2 = m;
          // size_t mid_triangle_0 = k + 1;
          // size_t top_triangle_0 = k + 2, top_triangle_1 = l + 1,
          //        top_triangle_2 = m + 1;

          ColumnarSimplex3D K0(bot_triangle_0, bot_triangle_1, bot_triangle_2,
                       mid_triangle_0);
          ColumnarSimplex3D K1(bot_triangle_1, top_triangle_2, bot_triangle_2,
                       mid_triangle_0);
          ColumnarSimplex3D K2(bot_triangle_1, top_triangle_2, mid_triangle_0,
                       top_triangle_1);
          ColumnarSimplex3D K3(mid_triangle_0, top_triangle_0, top_triangle_1,
                       top_triangle_2);

          // volume_mesh.cells.emplace_back(K0);
          // volume_mesh.cells.emplace_back(K1);
          // volume_mesh.cells.emplace_back(K2);
          // volume_mesh.cells.emplace_back(K3);

          _col_volume_mesh.cells[i].emplace_back(K0);
          _col_volume_mesh.cells[i].emplace_back(K1);
          _col_volume_mesh.cells[i].emplace_back(K2);
          _col_volume_mesh.cells[i].emplace_back(K3);
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

          // size_t bot_triangle_0 = k, bot_triangle_1 = l, bot_triangle_2 = m;
          // size_t mid_triangle_0 = k + 1, mid_triangle_1 = l + 1;
          // size_t top_triangle_0 = k + 2, top_triangle_1 = l + 2,
          //        top_triangle_2 = m + 1;

          ColumnarVertexIndex bot_triangle_0(face[0],k),
                              bot_triangle_1(face[1],l),
                              bot_triangle_2(face[2],m);
          ColumnarVertexIndex mid_triangle_0(face[0],k+1),
                              mid_triangle_1(face[1],l+1);
          ColumnarVertexIndex top_triangle_0(face[0],k+2),
                              top_triangle_1(face[1],l+2),
                              top_triangle_2(face[2],m+1);

          ColumnarSimplex3D K0(bot_triangle_0, bot_triangle_1, bot_triangle_2,
                       mid_triangle_1);
          ColumnarSimplex3D K1(bot_triangle_0, bot_triangle_2, mid_triangle_0,
                       mid_triangle_1);
          ColumnarSimplex3D K2(bot_triangle_2, top_triangle_2, mid_triangle_0,
                       mid_triangle_1);
          ColumnarSimplex3D K3(top_triangle_0, top_triangle_2, top_triangle_1,
                       mid_triangle_0);
          ColumnarSimplex3D K4(top_triangle_1, top_triangle_2, mid_triangle_1,
                       mid_triangle_0);

          // volume_mesh.cells.emplace_back(K0);
          // volume_mesh.cells.emplace_back(K1);
          // volume_mesh.cells.emplace_back(K2);
          // volume_mesh.cells.emplace_back(K3);
          // volume_mesh.cells.emplace_back(K4);

          _col_volume_mesh.cells[i].emplace_back(K0);
          _col_volume_mesh.cells[i].emplace_back(K1);
          _col_volume_mesh.cells[i].emplace_back(K2);
          _col_volume_mesh.cells[i].emplace_back(K3);
          _col_volume_mesh.cells[i].emplace_back(K4);
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

          // size_t bot_triangle_0 = k, bot_triangle_1 = l, bot_triangle_2 = m;
          // size_t mid_triangle_0 = k + 1, mid_triangle_1 = l + 1,
          //        mid_triangle_2 = m + 1;
          // size_t top_triangle_0 = k + 2, top_triangle_1 = l + 2,
          //        top_triangle_2 = m + 2;

          ColumnarVertexIndex bot_triangle_0(face[0],k),
                              bot_triangle_1(face[1],l),
                              bot_triangle_2(face[2],m);
          ColumnarVertexIndex mid_triangle_0(face[0],k+1),
                              mid_triangle_1(face[1],l+1),
                              mid_triangle_2(face[2],m+1);
          ColumnarVertexIndex top_triangle_0(face[0],k+2),
                              top_triangle_1(face[1],l+2),
                              top_triangle_2(face[2],m+2);

          ColumnarSimplex3D K0(bot_triangle_0, bot_triangle_1, bot_triangle_2,
                       mid_triangle_2);
          ColumnarSimplex3D K1(bot_triangle_0, mid_triangle_1, bot_triangle_1,
                       mid_triangle_2);
          ColumnarSimplex3D K2(bot_triangle_0, mid_triangle_0, mid_triangle_1,
                       mid_triangle_2);
          ColumnarSimplex3D K3(mid_triangle_0, mid_triangle_1, mid_triangle_2,
                       top_triangle_2);
          ColumnarSimplex3D K4(mid_triangle_0, top_triangle_1, mid_triangle_1,
                       top_triangle_2);
          ColumnarSimplex3D K5(mid_triangle_0, top_triangle_0, top_triangle_1,
                       top_triangle_2);

          // volume_mesh.cells.emplace_back(K0);
          // volume_mesh.cells.emplace_back(K1);
          // volume_mesh.cells.emplace_back(K2);
          // volume_mesh.cells.emplace_back(K3);
          // volume_mesh.cells.emplace_back(K4);
          // volume_mesh.cells.emplace_back(K5);

          _col_volume_mesh.cells[i].emplace_back(K0);
          _col_volume_mesh.cells[i].emplace_back(K1);
          _col_volume_mesh.cells[i].emplace_back(K2);
          _col_volume_mesh.cells[i].emplace_back(K3);
          _col_volume_mesh.cells[i].emplace_back(K4);
          _col_volume_mesh.cells[i].emplace_back(K5);
          volume_mesh_num_cells += 6;
        }
        break;
      }
      default:
        error("Face Coloring Error: Large layer height difference");
        break;
      }
    }
  info("Cells added to the Volume Mesh during domain padding: " + str(volume_mesh_num_cells));
  return _col_volume_mesh.to_volume_mesh();

}

};

} // namespace DTCC_BUILDER

#endif
