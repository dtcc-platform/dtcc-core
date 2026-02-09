// Copyright (C) 2018 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_MESH_BUILDER_H
#define DTCC_MESH_BUILDER_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <stack>
#include <tuple>
#include <vector>

#include "Geometry.h"
#include "Logging.h"
#include "MeshProcessor.h"
#include "Timer.h"
#include "Triangulate.h"
#include "VertexSmoother.h"
#include "model/GridField.h"
#include "model/Mesh.h"
#include "model/Surface.h"
#include "model/Vector.h"
// #include "spade/triangulate.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace DTCC_BUILDER
{

class MeshBuilder
{
public:
  static Mesh build_terrain_surface_mesh(const std::vector<Polygon> &subdomains,
                                         const std::vector<Polygon> &holes,
                                         const std::vector<double> &subdomain_triangle_size,
                                         const GridField &dtm, double max_mesh_size,
                                         double min_mesh_angle, size_t smooth_ground = 0,
                                         bool sort_triangles = false)
  {

    // Get bounding box
    const BoundingBox2D &bbox = dtm.grid.bounding_box;
    // build boundary
    Mesh ground_mesh =
        build_city_flat_mesh(subdomains, holes, subdomain_triangle_size, bbox.P.x, bbox.P.y,
                             bbox.Q.x, bbox.Q.y, max_mesh_size, min_mesh_angle, sort_triangles);
    // Displace ground surface. Fill all points with maximum height. This is
    // used to always choose the smallest height for each point since each point
    // may be visited multiple times.
    const double z_max = dtm.max();
    for (size_t i = 0; i < ground_mesh.vertices.size(); i++)
      ground_mesh.vertices[i].z = z_max;

    // If ground is not float, iterate over the triangles
    for (size_t i = 0; i < ground_mesh.faces.size(); i++)
    {
      // Get cell marker
      const int cell_marker = ground_mesh.markers[i];

      // Get triangle
      const Simplex2D &T = ground_mesh.faces[i];

      // Check cell marker
      if (cell_marker != -2) // not ground
      {
        // Compute minimum height of vertices
        double z_min = std::numeric_limits<double>::max();
        z_min = std::min(z_min, dtm(ground_mesh.vertices[T.v0]));
        z_min = std::min(z_min, dtm(ground_mesh.vertices[T.v1]));
        z_min = std::min(z_min, dtm(ground_mesh.vertices[T.v2]));

        // Set minimum height for all vertices
        set_min(ground_mesh.vertices[T.v0].z, z_min);
        set_min(ground_mesh.vertices[T.v1].z, z_min);
        set_min(ground_mesh.vertices[T.v2].z, z_min);
      }
      else
      {
        // Sample height map at vertex position for all vertices
        set_min(ground_mesh.vertices[T.v0].z, dtm(ground_mesh.vertices[T.v0]));
        set_min(ground_mesh.vertices[T.v1].z, dtm(ground_mesh.vertices[T.v1]));
        set_min(ground_mesh.vertices[T.v2].z, dtm(ground_mesh.vertices[T.v2]));
      }
    }
    info("smooth ground...");
    if (smooth_ground > 0)
      VertexSmoother::smooth_mesh(ground_mesh, smooth_ground, true, true);

    // for (size_t i = 0; i < ground_mesh.faces.size(); i++)
    // {
    //   auto normal = Geometry::face_normal(ground_mesh.faces[i], ground_mesh);
    //   if (normal.z < 0)
    //     ground_mesh.faces[i].flip();
    // }

    info("ground mesh done");
    return ground_mesh;
  }

  // Build ground mesh for city.
  //
  // The mesh is a triangular mesh of the rectangular region
  // defined by (xmin, xmax) x (ymin, ymax). The edges of the mesh respect
  // the boundaries of the buildings.
  //
  // markers:
  //
  // -2: ground (cells outside buildings and halos)
  // -1: halos (cells close to buildings)
  //  0: building 0 (cells inside building 0)
  //  1: building 1 (cells inside building 1)
  //  etc (non-negative integers mark cells inside buildings)
  static Mesh build_city_flat_mesh(const std::vector<Polygon> &subdomains,
                                   const std::vector<Polygon> &holes,
                                   const std::vector<double> &subdomain_triangle_size, double xmin,
                                   double ymin, double xmax, double ymax, double max_mesh_size,
                                   double min_mesh_angle, bool sort_triangles = false)
  {
    info("Building city flat mesh...");
    Timer timer("build_city_flat_mesh");

    const BoundingBox2D bounding_box(Vector2D(xmin, ymin), Vector2D(xmax, ymax));
    const size_t nx = static_cast<size_t>((bounding_box.Q.x - bounding_box.P.x) / max_mesh_size);
    const size_t ny = static_cast<size_t>((bounding_box.Q.y - bounding_box.P.y) / max_mesh_size);
    const size_t n = nx * ny;
    info("Bounds: " + str(bounding_box));
    info("Max mesh size: " + str(max_mesh_size));
    info("Estimated number of faces: " + str(n));
    info("Number of subdomains (buildings): " + str(subdomains.size()));
    info("Number of explicit holes: " + str(holes.size()));

    std::vector<std::vector<Vector2D>> triangle_sub_domains;
    triangle_sub_domains.reserve(subdomains.size());
    for (const auto &sd : subdomains)
    {
      if (!sd.vertices.empty())
        triangle_sub_domains.push_back(sd.vertices);
      for (const auto &hole : sd.holes)
      {
        if (!hole.empty())
          triangle_sub_domains.push_back(hole);
      }
    }
    info("Number of subdomains (buildings + building holes): " + str(triangle_sub_domains.size()));

    std::vector<std::vector<Vector2D>> triangle_holes;
    triangle_holes.reserve(holes.size());
    for (const auto &hole_polygon : holes)
    {
      if (!hole_polygon.vertices.empty())
        triangle_holes.push_back(hole_polygon.vertices);
      for (const auto &nested : hole_polygon.holes)
      {
        if (!nested.empty())
          triangle_holes.push_back(nested);
      }
    }
    info("Number of explicit hole loops (including nested): " + str(triangle_holes.size()));

    std::vector<Vector2D> boundary{};
    boundary.push_back(bounding_box.P);
    boundary.push_back(Vector2D(bounding_box.Q.x, bounding_box.P.y));
    boundary.push_back(bounding_box.Q);
    boundary.push_back(Vector2D(bounding_box.P.x, bounding_box.Q.y));

    Mesh mesh;
#ifdef DTCC_HAVE_TRIANGLE
    info("Triangulation Backend: Triangle");
    Triangulate::call_triangle(mesh, boundary, triangle_sub_domains, triangle_holes,
                               subdomain_triangle_size, max_mesh_size, min_mesh_angle,
                               sort_triangles);
#else
    info("Triangulation Backend: Spade");
    double effective_maxh = max_mesh_size;
    if (!subdomain_triangle_size.empty())
    {
      double min_subdomain = std::numeric_limits<double>::max();
      for (double h : subdomain_triangle_size)
      {
        if (h > 0.0)
          min_subdomain = std::min(min_subdomain, h);
      }
      if (min_subdomain < std::numeric_limits<double>::max())
      {
        if (effective_maxh > 0.0)
          effective_maxh = std::min(effective_maxh, min_subdomain);
        else
          effective_maxh = min_subdomain;
      }
    }

    Triangulate::call_spade(mesh, boundary, triangle_holes, triangle_sub_domains, effective_maxh,
                            min_mesh_angle, sort_triangles);
#endif

    MeshProcessor::compute_mesh_domain_markers(mesh, subdomains);

    return mesh;
  }

  // Layer ground mesh to create a volume mesh.
  //
  // The volume mesh is a tetrahedral mesh constructed
  // extruding the 2D mesh in the vertical (z) direction.
  //
  // markers:
  //
  // -2: ground (cells outside buildings and halos)
  // -1: halos (cells close to buildings)
  //  0: building 0 (cells inside building 0)
  //  1: building 1 (cells inside building 1)
  //  etc (non-negative integers mark cells inside buildings)
  //
  // Note that the markers are just propagatated upward from the
  // ground mesh, meaning that the markers will be the same in each
  // column of the ground mesh.
  //
  // It is assumed that the ground mesh is sorted, which is the case if the
  // mesh has been built by calling build_mesh().
  static VolumeMesh layer_ground_mesh(const Mesh &ground_mesh,
                                      const std::vector<double> layer_heights)
  {
    Timer timer("build_volume_mesh");

    const size_t num_layers = layer_heights.size();

    const size_t layer_size = ground_mesh.vertices.size();

    info("Building volume mesh with " + str(num_layers) + " layers...");

    // Initialize volume mesh
    VolumeMesh volume_mesh;
    volume_mesh.vertices.resize((num_layers + 1) * ground_mesh.vertices.size());
    volume_mesh.cells.resize(num_layers * 3 * ground_mesh.faces.size());
    volume_mesh.markers.resize(volume_mesh.cells.size());
    volume_mesh.num_layers = num_layers;

    // Add vertices
    {
      size_t k = 0;
      double z = 0;

      // set zero layer
      for (const auto &p_2d : ground_mesh.vertices)
        volume_mesh.vertices[k++] = Vector3D(p_2d.x, p_2d.y, z);
      for (const auto &h : layer_heights)
      {
        z += h;
        for (const auto &p_2d : ground_mesh.vertices)
          volume_mesh.vertices[k++] = Vector3D(p_2d.x, p_2d.y, z);
      }
    }

    // Add cells
    {
      size_t k = 0;
      size_t offset = 0;
      for (size_t layer = 0; layer < num_layers; layer++)
      {
        // Iterate over triangles in layer
        for (const auto &T : ground_mesh.faces)
        {
          // Get sorted vertex indices for bottom layer
          const size_t u0 = T.v0 + offset;
          const size_t u1 = T.v1 + offset;
          const size_t u2 = T.v2 + offset;

          // Get sorted vertices for top layer
          const size_t v0 = u0 + layer_size;
          const size_t v1 = u1 + layer_size;
          const size_t v2 = u2 + layer_size;

          // Create three tetrahedra by connecting the first vertex
          // of each edge in the bottom layer with the second
          // vertex of the corresponding edge in the top layer.
          volume_mesh.cells[k++] = Simplex3D(u0, u1, u2, v2);
          volume_mesh.cells[k++] = Simplex3D(u0, v1, u1, v2);
          volume_mesh.cells[k++] = Simplex3D(u0, v0, v1, v2);
        }

        // Add to offset
        offset += layer_size;
      }
    }

    // Add domain markers
    {
      size_t k = 0;
      for (size_t layer = 0; layer < num_layers; layer++)
      {
        for (const auto &marker : ground_mesh.markers)
        {
          int m = 0;

          // Top layer marked as -3
          if (layer == num_layers - 1)
            m = -3;

          // Halo and ground only marked for bottom layer
          else if (marker == -1 || marker == -2)
          {
            if (layer == 0)
              m = marker;
            else
              m = -4;
          }

          // buildings marked for all layers (except top layer).
          // Later adjusted to -4 above buildings in trim_volume_mesh.
          else
          {
            m = marker;
          }

          // Set markers
          volume_mesh.markers[k++] = m;
          volume_mesh.markers[k++] = m;
          volume_mesh.markers[k++] = m;
        }
      }
    }

    return volume_mesh;
  }

  // Trim volume mesh by removing cells inside buildings.
  //
  // markers:
  //
  // -4: other (not top, ground, halo, or building)
  // -3: top (cells in *top layer*)
  // -2: ground (cells in *bottom layer* outside buildings and halos)
  // -1: halos (cells in *bottom layer* close to buildings)
  //  0: building 0 (cells *first layer* inside building 0)
  //  1: building 1 (cells *first layer* inside building 1)
  //  etc (non-negative integers mark cells inside buildings)
  //
  // Note that the markers are adjusted in the following way:
  //
  // - Only cells in bottom layer marked as ground (-2) or halo (-1)
  // - Only cells in first layer above a building marked as building
  // - cells in top layer marked as top (-3)
  // - All other cells (in between) marked as other (-4)
  static VolumeMesh trim_volume_mesh(const VolumeMesh &volume_mesh, const Mesh &mesh,
                                     const std::vector<Surface> &buildings)
  {
    info("Trimming volume mesh...");
    Timer timer("trim_volume_mesh");

    // Get sizes
    const size_t num_buildings = buildings.size();
    const size_t num_cells_2d = mesh.faces.size();
    const size_t num_cells_3d = volume_mesh.cells.size();
    const size_t layer_size = 3 * mesh.faces.size();

    // Phase 1: Determine which cells should be trimmed
    // ------------------------------------------------

    // build map from buildings to cells in 2D mesh
    std::vector<std::vector<size_t>> building_cells_2d(num_buildings);
    for (size_t cell_index_2d = 0; cell_index_2d < num_cells_2d; cell_index_2d++)
    {
      const int building_index = mesh.markers[cell_index_2d];
      if (building_index >= 0)
        building_cells_2d[building_index].push_back(cell_index_2d);
    }

    // Create markers for cells to be trimmed (keep by default)
    std::vector<bool> trim_cell(num_cells_3d);
    std::fill(trim_cell.begin(), trim_cell.end(), false);

    // Keep track of first layer for each building
    std::vector<size_t> first_layer(num_buildings);
    std::fill(first_layer.begin(), first_layer.end(), 0);

    // Iterate over buildings
    for (size_t building_index = 0; building_index < num_buildings; building_index++)
    {
      // Iterate over layers
      for (size_t layer = 0; layer < volume_mesh.num_layers; layer++)
      {
        // build list of 3D cells for building in current layer
        std::vector<size_t> cells_3d;
        for (const auto &cell_index_2d : building_cells_2d[building_index])
        {
          for (size_t j = 0; j < 3; j++)
            cells_3d.push_back(index_3d(layer, layer_size, cell_index_2d, j));
        }

        // Trim layer if any cell midpoint is below building height
        bool trim_layer = false;
        for (const auto &cell_index_3d : cells_3d)
        {
          const double z = volume_mesh.mid_point(cell_index_3d).z;
          const double h = buildings[building_index].max_height();
          if (z < h)
          {
            trim_layer = true;
            break;
          }
        }

        // Check if layer should be trimmed
        if (trim_layer)
        {
          // Mark cells for trimming
          for (const auto &cell_index_3d : cells_3d)
            trim_cell[cell_index_3d] = true;
        }
        else
        {
          // If layer should be kept, no need to check more layers
          first_layer[building_index] = layer;
          break;
        }
      }
    }

    // Phase 2: Adjust markers
    // -----------------------

    // Create copy of markeres
    std::vector<int> markers{volume_mesh.markers};

    // Mark cells between bottom and top layer as -4
    for (size_t layer = 1; layer < volume_mesh.num_layers - 1; layer++)
    {
      for (size_t cell_index_2d = 0; cell_index_2d < num_cells_2d; cell_index_2d++)
        for (size_t j = 0; j < 3; j++)
          markers[index_3d(layer, layer_size, cell_index_2d, j)] = -4;
    }

    // Mark cells in top layer as -3
    for (size_t cell_index_2d = 0; cell_index_2d < num_cells_2d; cell_index_2d++)
    {
      for (size_t j = 0; j < 3; j++)
        markers[index_3d(volume_mesh.num_layers - 1, layer_size, cell_index_2d, j)] = -3;
    }

    // Mark cells in first layer above each building:
    //
    // 0, 1, 2, ... if building is not covered by bottom layer (normal case)
    // -1           if building is covered by bottom layer (modify to halo)
    for (size_t building_index = 0; building_index < num_buildings; building_index++)
    {
      const size_t layer = first_layer[building_index];
      size_t marker = building_index;
      if (layer == 0)
      {
        warning("Building " + str(building_index) + " is covered by bottom layer");
        marker = -1;
      }
      for (const auto &cell_index_2d : building_cells_2d[building_index])
      {
        for (size_t j = 0; j < 3; j++)
          markers[index_3d(layer, layer_size, cell_index_2d, j)] = marker;
      }
    }

    // Phase 3: Extract new mesh for all cells that should be kept
    // -----------------------------------------------------------

    // Renumber vertices and cells
    std::unordered_map<size_t, size_t> vertex_map;
    std::unordered_map<size_t, size_t> cell_map;
    size_t k = 0;
    size_t l = 0;
    for (size_t cell_index_3d = 0; cell_index_3d < volume_mesh.cells.size(); cell_index_3d++)
    {
      if (!trim_cell[cell_index_3d])
      {
        // Get cell
        const Simplex3D &T = volume_mesh.cells[cell_index_3d];

        // Renumbers vertices
        if (vertex_map.find(T.v0) == vertex_map.end())
          vertex_map[T.v0] = k++;
        if (vertex_map.find(T.v1) == vertex_map.end())
          vertex_map[T.v1] = k++;
        if (vertex_map.find(T.v2) == vertex_map.end())
          vertex_map[T.v2] = k++;
        if (vertex_map.find(T.v3) == vertex_map.end())
          vertex_map[T.v3] = k++;

        // Renumber cells
        cell_map[cell_index_3d] = l++;
      }
    }

    // Initialize new mesh data
    const size_t num_vertices = vertex_map.size();
    const size_t num_cells = cell_map.size();
    std::vector<Vector3D> _vertices(num_vertices);
    std::vector<Simplex3D> _cells(num_cells);
    std::vector<int> _markers(num_cells);

    // Set new mesh data
    for (const auto v : vertex_map)
      _vertices[v.second] = volume_mesh.vertices[v.first];
    for (const auto c : cell_map)
    {
      _cells[c.second].v0 = vertex_map[volume_mesh.cells[c.first].v0];
      _cells[c.second].v1 = vertex_map[volume_mesh.cells[c.first].v1];
      _cells[c.second].v2 = vertex_map[volume_mesh.cells[c.first].v2];
      _cells[c.second].v3 = vertex_map[volume_mesh.cells[c.first].v3];
      _markers[c.second] = markers[c.first];
    }

    // Create new mesh and assign data
    VolumeMesh _volume_mesh;
    _volume_mesh.vertices = _vertices;
    _volume_mesh.cells = _cells;
    _volume_mesh.markers = _markers;

    return _volume_mesh;
  }

  static std::vector<Mesh>
  build_city_surface_mesh(const std::vector<Surface> &buildings, const std::vector<Surface> &holes,
                          const std::vector<int> meshing_directive,
                          const std::vector<double> &subdomain_triangle_size, const GridField &dtm,
                          double max_mesh_size, double min_mesh_angle, size_t smooth_ground = 0,
                          bool merge_meshes = true, bool sort_triangles = false)
  {
    auto build_city_surface_t = Timer("build_city_surface_mesh");
    auto terrain_time = Timer("build_city_surface_mesh: step 1 terrain");
    const size_t num_buildings = buildings.size();
    std::vector<Polygon> subdomains;
    subdomains.reserve(buildings.size());
    for (const auto &b : buildings)
    {
      subdomains.push_back(b.to_polygon());
    }
    std::vector<Polygon> hole_domains;
    hole_domains.reserve(holes.size());
    for (const auto &h : holes)
    {
      hole_domains.push_back(h.to_polygon());
    }
    Mesh terrain_mesh =
        build_terrain_surface_mesh(subdomains, hole_domains, subdomain_triangle_size, dtm,
                                   max_mesh_size, min_mesh_angle, smooth_ground, sort_triangles);
    terrain_time.stop();

    std::vector<Mesh> city_mesh;
    std::vector<Mesh> building_meshes;

    std::map<size_t, std::vector<Simplex2D>> building_faces;
    std::vector<size_t> building_indices;

    std::map<size_t, std::vector<Simplex2D>> platform_faces;
    std::map<size_t, double> platform_min_z;

    info("finding markes");
    auto find_markers_t = Timer("build_city_surface_mesh: step 2 find markers");
    for (size_t i = 0; i < terrain_mesh.markers.size(); i++)
    {
      auto marker = terrain_mesh.markers[i];
      if (marker < 0)
        continue;

      const auto &face = terrain_mesh.faces[i];

      if (meshing_directive[marker] > 0)
      {
        building_faces[marker].push_back(face);
        building_indices.push_back(i);
      }
      else
      {
        platform_faces[marker].push_back(face);

        const auto &v0 = terrain_mesh.vertices[face.v0];
        const auto &v1 = terrain_mesh.vertices[face.v1];
        const auto &v2 = terrain_mesh.vertices[face.v2];

        auto [it, _] = platform_min_z.try_emplace(marker, std::numeric_limits<double>::infinity());
        double &minz = it->second;
        minz = std::min({minz, v0.z, v1.z, v2.z});
      }
    }

    find_markers_t.stop();

    info("build_city_surface_mesh: step 3 flatten platforms");
    for (const auto &kv : platform_faces)
    {
      const size_t marker = kv.first;
      const auto &faces = kv.second;
      auto it = platform_min_z.find(marker);
      if (it == platform_min_z.end())
        continue; // no faces? skip
      const double zflat = it->second;

      // Collect unique vertex indices used by these faces
      std::unordered_set<size_t> vset;
      vset.reserve(faces.size() * 3);
      for (const auto &f : faces)
      {
        vset.insert(static_cast<size_t>(f.v0));
        vset.insert(static_cast<size_t>(f.v1));
        vset.insert(static_cast<size_t>(f.v2));
      }
      // Set their z to zflat
      for (size_t vi : vset)
      {
        terrain_mesh.vertices[vi].z = zflat;
      }
    }

    if (smooth_ground)
    {
      VertexSmoother::smooth_mesh(terrain_mesh, smooth_ground, false, true);
    }

    info("building meshes");
    auto building_meshes_t = Timer("build_city_surface_mesh: step 3 building meshes");
    for (const auto &kv : building_faces)
    {
      const auto marker = kv.first;
      const auto &faces = kv.second;
      const auto &building = buildings[marker];
      const auto roof_height = building.max_height();

      auto naked_edges = MeshProcessor::find_naked_edges(faces);

      Mesh building_mesh;
      building_mesh.vertices.reserve(faces.size() * 3 + naked_edges.size() * 4);
      building_mesh.faces.reserve(faces.size() + naked_edges.size() * 2);
      building_mesh.markers.reserve(faces.size() + naked_edges.size() * 2);

      // add roofs
      for (const auto &face : faces)
      {
        const auto v0 = terrain_mesh.vertices[face.v0];
        const auto v1 = terrain_mesh.vertices[face.v1];
        const auto v2 = terrain_mesh.vertices[face.v2];
        const auto v3 = Vector3D(v0.x, v0.y, roof_height);
        const auto v4 = Vector3D(v1.x, v1.y, roof_height);
        const auto v5 = Vector3D(v2.x, v2.y, roof_height);

        const auto num_vertices = building_mesh.vertices.size();
        building_mesh.vertices.push_back(v3);
        building_mesh.vertices.push_back(v4);
        building_mesh.vertices.push_back(v5);
        // info("vertices: " + str(v3) + " " + str(v4) + " " + str(v5));
        // info("faces: " + str(face.v0) + " " + str(face.v1) + " " +
        //     str(face.v2));
        building_mesh.faces.push_back(Simplex2D(num_vertices, num_vertices + 1, num_vertices + 2));
        building_mesh.markers.push_back(static_cast<int>(num_buildings + marker));
      }

      // add walls
      for (const auto &edge_faces : naked_edges)
      {
        const Simplex1D edge = edge_faces.first;
        const Simplex2D edge_face = edge_faces.second;
        const auto face_center = Geometry::face_center(edge_face, terrain_mesh);

        const auto ground_v0 = terrain_mesh.vertices[edge.v0];
        const auto ground_v1 = terrain_mesh.vertices[edge.v1];
        const auto roof_v0 = Vector3D(ground_v0.x, ground_v0.y, roof_height);
        const auto roof_v1 = Vector3D(ground_v1.x, ground_v1.y, roof_height);

        const auto base = building_mesh.vertices.size();
        building_mesh.vertices.push_back(ground_v0);
        building_mesh.vertices.push_back(ground_v1);
        building_mesh.vertices.push_back(roof_v0);
        building_mesh.vertices.push_back(roof_v1);

        const auto wall_normal = Geometry::triangle_normal(ground_v0, ground_v1, roof_v1);
        if (Geometry::dot_3d(wall_normal, face_center - ground_v0) > 0)
        {
          building_mesh.faces.push_back(Simplex2D(base, base + 3, base + 1));
          building_mesh.faces.push_back(Simplex2D(base, base + 2, base + 3));
          building_mesh.markers.push_back(static_cast<int>(marker));
          building_mesh.markers.push_back(static_cast<int>(marker));
        }
        else
        {
          building_mesh.faces.push_back(Simplex2D(base, base + 1, base + 3));
          building_mesh.faces.push_back(Simplex2D(base, base + 3, base + 2));
          building_mesh.markers.push_back(static_cast<int>(marker));
          building_mesh.markers.push_back(static_cast<int>(marker));
        }
      }

      building_meshes.push_back(MeshProcessor::weld_mesh(building_mesh));
    }
    building_meshes_t.stop();

    // remove triangles inside houses from terrain
    auto remove_inside_t = Timer("build_city_surface_mesh: step 4 remove inside");
    std::sort(building_indices.begin(), building_indices.end());

    // keep faces/markers/normals aligned while removing building triangles
    std::vector<Simplex2D> filtered_faces;
    std::vector<int> filtered_markers;
    std::vector<Vector3D> filtered_normals;
    const bool copy_normals = terrain_mesh.normals.size() == terrain_mesh.faces.size();
    filtered_faces.reserve(terrain_mesh.faces.size());
    filtered_markers.reserve(terrain_mesh.markers.size());
    if (copy_normals)
      filtered_normals.reserve(terrain_mesh.normals.size());

    size_t bpos = 0;
    for (size_t i = 0; i < terrain_mesh.faces.size(); ++i)
    {
      while (bpos < building_indices.size() && building_indices[bpos] < i)
        ++bpos;
      if (bpos < building_indices.size() && building_indices[bpos] == i)
      {
        ++bpos;
        continue;
      }
      filtered_faces.push_back(terrain_mesh.faces[i]);
      filtered_markers.push_back(terrain_mesh.markers[i]);
      if (copy_normals)
        filtered_normals.push_back(terrain_mesh.normals[i]);
    }

    terrain_mesh.faces.swap(filtered_faces);
    terrain_mesh.markers.swap(filtered_markers);
    if (copy_normals)
      terrain_mesh.normals.swap(filtered_normals);
    remove_inside_t.stop();

    auto final_merger_t = Timer("build_city_surface_mesh: step 5 final merge");
    city_mesh.push_back(terrain_mesh);
    city_mesh.insert(city_mesh.end(), building_meshes.begin(), building_meshes.end());
    if (merge_meshes)
    {
      auto merged_mesh = MeshProcessor::merge_meshes(city_mesh, true);
      city_mesh = {merged_mesh};
    }
    final_merger_t.stop();
    build_city_surface_t.stop();
    // Timer::report("city surface");
    return city_mesh;
  }

  static Mesh mesh_surface(const Surface &surface,

                           double max_triangle_area_size = -1, double min_mesh_angle = 25)
  // Convert 3D Surface to triangle Mesh.
  // - If max_triangle_area_size < 0: uses fast_mesh (earcut/fan triangulation)
  // - If max_triangle_area_size >= 0:
  //   - With DTCC_HAVE_TRIANGLE: uses Triangle library (preferred)
  //   - Without DTCC_HAVE_TRIANGLE: uses SPADE library (default)
  {
    Mesh mesh;
    if (surface.vertices.size() < 3)
      return mesh;
    if (max_triangle_area_size < 0)
    {
      Triangulate::fast_mesh(mesh, surface);
    }
    else
    {
#ifdef DTCC_HAVE_TRIANGLE
      Triangulate::call_triangle(mesh, surface, max_triangle_area_size, min_mesh_angle);
#else
      Triangulate::call_spade(mesh, surface, max_triangle_area_size, min_mesh_angle);
#endif
    }
    return mesh;
  }

  static Mesh mesh_multisurface(const MultiSurface &multi_surface,
                                double max_triangle_area_size = -1, double min_mesh_angle = 25,
                                bool weld = false, double snap = 0)
  {
    std::vector<Mesh> multimesh(multi_surface.surfaces.size());
    //    info("meshing multisurface with " + str(multi_surface.surfaces.size())
    //    +
    //         " surfaces");
    // #pragma omp parallel for
    for (size_t i = 0; i < multi_surface.surfaces.size(); i++)
    {
      multimesh[i] =
          mesh_surface(multi_surface.surfaces[i], max_triangle_area_size, min_mesh_angle);
    }
    // for (const auto &surface : multi_surface.surfaces)
    // {
    //   auto surface_mesh =
    //       mesh_surface(surface, max_triangle_area_size, min_mesh_angle);
    //   multimesh.push_back(surface_mesh);
    // }
    auto mesh = MeshProcessor::merge_meshes(multimesh, weld, snap);
    mesh.normalize_normal_direction();
    return mesh;
  }

  static std::vector<Mesh> mesh_multisurfaces(const std::vector<MultiSurface> &multi_surfaces,
                                              double max_triangle_area_size = -1,
                                              double min_mesh_angle = 25, bool weld = false)
  {
    int n = multi_surfaces.size();
    std::vector<Mesh> meshes(n);
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
      auto mesh =
          mesh_multisurface(multi_surfaces[i], max_triangle_area_size, min_mesh_angle, weld);
      meshes[i] = mesh;
    }

    return meshes;
  }

private:
  // Map from 2D cell index to 3D cell indices
  static size_t index_3d(size_t layer, size_t layer_size, size_t cell_index_2d, size_t j)
  {
    return layer * layer_size + 3 * cell_index_2d + j;
  }

  // Compute domain markers for subdomains
  static void compute_domain_markers(Mesh &mesh, const std::vector<Polygon> &subdomains)
  {
    info("Computing domain markers...");
    Timer timer("compute_domain_markers");

    // build search tree for subdomains

    auto search_tree = BoundingBoxTree2D();
    std::vector<BoundingBox2D> bounding_boxes;
    for (const auto &subdomain : subdomains)
    {
      bounding_boxes.push_back(BoundingBox2D(subdomain));
    }
    search_tree.build(bounding_boxes);

    // Initialize domain markers and set all markers to -2 (ground)
    mesh.markers.resize(mesh.faces.size());
    std::fill(mesh.markers.begin(), mesh.markers.end(), -2);

    // Initialize markers for vertices belonging to a building
    std::vector<bool> is_building_vertex(mesh.vertices.size());
    std::fill(is_building_vertex.begin(), is_building_vertex.end(), false);

    // Iterate over cells to mark buildings
    if (subdomains.size() > 0)
    {
      for (size_t i = 0; i < mesh.faces.size(); i++)
      {
        // find building containg midpoint of cell (if any)
        const Vector3D c_3d = mesh.mid_point(i);
        const Vector2D c_2d(c_3d.x, c_3d.y);
        std::vector<size_t> indices = search_tree.find(Vector2D(c_2d));

        if (indices.size() > 0)
        {
          for (const auto &index : indices)
          {
            if (Geometry::polygon_contains_2d(subdomains[index], c_2d))
            {
              mesh.markers[i] = index;
              const Simplex2D &T = mesh.faces[i];
              // Mark all cell vertices as belonging to a building
              is_building_vertex[T.v0] = true;
              is_building_vertex[T.v1] = true;
              is_building_vertex[T.v2] = true;

              // // Check if individual vertices are inside a building
              // // (not only midpoint). Necessary for when building
              // // visualization meshes that are not boundary-fitted.
              // if (search_tree.find(mesh.vertices[T.v0]).size() == 0)
              //   is_building_vertex[T.v0] = false;
              // if (search_tree.find(mesh.vertices[T.v1]).size() == 0)
              //   is_building_vertex[T.v1] = false;
              // if (search_tree.find(mesh.vertices[T.v2]).size() == 0)
              //   is_building_vertex[T.v2] = false;

              break;
            }
          }
        }
      }

      // Iterate over cells to mark building halos
      for (size_t i = 0; i < mesh.faces.size(); i++)
      {
        // Check if any of the cell vertices belongs to a building
        const Simplex2D &T = mesh.faces[i];
        const bool touches_building =
            (is_building_vertex[T.v0] || is_building_vertex[T.v1] || is_building_vertex[T.v2]);

        // Mark as halo (-1) if the cell touches a building but is not
        // itself inside footprint (not marked in the previous step)
        if (touches_building && mesh.markers[i] == -2)
          mesh.markers[i] = -1;
      }
    }
  }

  // Set x = min(x, y)
  static void set_min(double &x, double y)
  {
    if (y < x)
      x = y;
  }
};

} // namespace DTCC_BUILDER

#endif
