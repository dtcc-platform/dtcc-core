// Copyright (C) 2018 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_MESH_BUILDER_H
#define DTCC_MESH_BUILDER_H

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <stack>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
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
    return build_terrain_surface_mesh_from_ground_mesh(ground_mesh, dtm, smooth_ground);
  }

  static Mesh build_terrain_surface_mesh_from_ground_mesh(
      Mesh ground_mesh,
      const GridField &dtm,
      size_t smooth_ground = 0)
  {
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
    info("Surface mesh | smoothing ground surface");
    if (smooth_ground > 0)
      VertexSmoother::smooth_mesh(ground_mesh, smooth_ground, true, true);

    // for (size_t i = 0; i < ground_mesh.faces.size(); i++)
    // {
    //   auto normal = Geometry::face_normal(ground_mesh.faces[i], ground_mesh);
    //   if (normal.z < 0)
    //     ground_mesh.faces[i].flip();
    // }

    info("Surface mesh | ground surface ready");
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
                                   double min_mesh_angle, bool sort_triangles = false,
                                   const std::string &backend = "auto")
  {
    info("Building city flat mesh...");
    Timer timer("build_city_flat_mesh");

    const BoundingBox2D bounding_box(Vector2D(xmin, ymin), Vector2D(xmax, ymax));
    const bool has_global_mesh_limit = max_mesh_size > 0.0;
    const size_t nx = has_global_mesh_limit ?
        static_cast<size_t>((bounding_box.Q.x - bounding_box.P.x) / max_mesh_size) : 0;
    const size_t ny = has_global_mesh_limit ?
        static_cast<size_t>((bounding_box.Q.y - bounding_box.P.y) / max_mesh_size) : 0;
    const size_t n = nx * ny;
    info("Bounds: " + str(bounding_box));
    info("Max mesh size: " + (has_global_mesh_limit ? str(max_mesh_size) : std::string("unrestricted")));
    if (has_global_mesh_limit)
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
    resolve_2d_backend(backend);
    info("Triangulation backend: triangle");
    Triangulate::call_triangle(mesh, boundary, triangle_sub_domains, triangle_holes,
                               subdomain_triangle_size, max_mesh_size, min_mesh_angle,
                               sort_triangles);

    MeshProcessor::compute_mesh_domain_markers(mesh, subdomains);

    return mesh;
  }

  static std::vector<Mesh>
  build_city_surface_mesh_from_terrain_mesh(
      const std::vector<Surface> &buildings,
      const std::vector<int> meshing_directive,
      Mesh terrain_mesh,
      size_t smooth_ground = 0,
      bool merge_meshes = true)
  {
    auto build_city_surface_t = Timer("build_city_surface_mesh");
    const size_t num_buildings = buildings.size();

    if (meshing_directive.size() != num_buildings)
    {
      throw std::invalid_argument(
          "build_city_surface_mesh_from_terrain_mesh requires one meshing directive per building.");
    }

    std::vector<Mesh> city_mesh;
    std::vector<Mesh> building_meshes;

    std::map<size_t, std::vector<Simplex2D>> building_faces;
    std::vector<size_t> building_indices;

    std::map<size_t, std::vector<Simplex2D>> platform_faces;
    std::map<size_t, double> platform_min_z;
    std::map<size_t, double> building_min_z;

    info("Surface mesh | assigning region markers");
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

        const auto &v0 = terrain_mesh.vertices[face.v0];
        const auto &v1 = terrain_mesh.vertices[face.v1];
        const auto &v2 = terrain_mesh.vertices[face.v2];

        auto [it, _] = building_min_z.try_emplace(marker, std::numeric_limits<double>::infinity());
        double &minz = it->second;
        minz = std::min({minz, v0.z, v1.z, v2.z});
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

    // Buildings whose footprints share a vertex form a block. A block has one
    // base level, and the walls between its buildings are step walls.
    const auto blocks = building_blocks(building_faces);
    for (const auto &block : blocks)
    {
      double base = std::numeric_limits<double>::infinity();
      for (const auto marker : block)
        base = std::min(base, building_min_z[marker]);
      for (const auto marker : block)
        building_min_z[marker] = base;
    }

    info("Surface mesh | flattening supported regions");
    std::vector<char> platform_freeze(terrain_mesh.vertices.size(), 0);
    const auto flatten_region_faces =
        [&](const std::map<size_t, std::vector<Simplex2D>> &region_faces,
            const std::map<size_t, double> &region_min_z)
    {
      for (const auto &kv : region_faces)
      {
        const size_t marker = kv.first;
        const auto &faces = kv.second;
        auto it = region_min_z.find(marker);
        if (it == region_min_z.end())
          continue;
        const double zflat = it->second;

        std::unordered_set<size_t> vset;
        vset.reserve(faces.size() * 3);
        for (const auto &f : faces)
        {
          vset.insert(static_cast<size_t>(f.v0));
          vset.insert(static_cast<size_t>(f.v1));
          vset.insert(static_cast<size_t>(f.v2));
        }

        for (size_t vi : vset)
        {
          terrain_mesh.vertices[vi].z = zflat;
          platform_freeze[vi] = 1;
        }
      }
    };
    flatten_region_faces(platform_faces, platform_min_z);
    flatten_region_faces(building_faces, building_min_z);

    if (smooth_ground)
    {
      VertexSmoother::smooth_mesh(terrain_mesh, smooth_ground, platform_freeze, true);
    }

    info("Surface mesh | building shell meshes");
    auto building_meshes_t = Timer("build_city_surface_mesh: step 3 building meshes");

    // Roof level of each building. Roofs in a block that differ by less than
    // a tolerance are unified, so they meet without a sliver step.
    std::map<size_t, double> roof_z;
    for (const auto &block : blocks)
    {
      std::vector<std::pair<double, size_t>> roofs;
      for (const auto marker : block)
        roofs.emplace_back(buildings[marker].max_height(), marker);
      std::sort(roofs.begin(), roofs.end());
      size_t start = 0;
      for (size_t i = 1; i <= roofs.size(); ++i)
      {
        if (i < roofs.size() &&
            roofs[i].first - roofs[start].first <= roof_tolerance(roofs[start].first))
          continue;
        for (size_t j = start; j < i; ++j)
          roof_z[roofs[j].second] = roofs[i - 1].first;
        start = i;
      }
    }

    // The building markers on each edge of a building face.
    std::unordered_map<Simplex1D, std::vector<size_t>, Simplex1DHash> edge_markers;
    for (const auto &kv : building_faces)
    {
      for (const auto &face : kv.second)
      {
        for (int i = 0; i < 3; ++i)
        {
          auto &markers = edge_markers[Simplex1D(face[i], face[(i + 1) % 3], true)];
          if (std::find(markers.begin(), markers.end(), kv.first) == markers.end())
            markers.push_back(kv.first);
        }
      }
    }

    // Walls. An edge on the outline of a block gets a wall from the base to
    // the roof. An edge between two buildings of a block gets a step wall from
    // the lower roof to the higher, owned by the taller building, and none
    // when the roofs are equal. All walls of a block are split at the same
    // levels (the base, every roof, and strips between them), so walls and
    // roofs meeting at a vertical line or a roof edge share their vertices.
    std::map<size_t, std::vector<BuildingWall>> walls;
    std::map<size_t, std::vector<double>> wall_levels;
    for (const auto &block : blocks)
    {
      const double base = building_min_z[block.front()];
      std::vector<double> levels = {base};
      std::vector<double> horizontal_lengths;
      for (const auto marker : block)
      {
        levels.push_back(roof_z[marker]);
        for (const auto &edge_face : MeshProcessor::find_naked_edges(building_faces[marker]))
        {
          const Simplex1D edge = edge_face.first;
          double z0 = base;
          bool shared = false;
          for (const auto other : edge_markers[edge])
          {
            if (other == marker)
              continue;
            shared = true;
            z0 = roof_z[other];
          }
          if (shared && z0 >= roof_z[marker])
            continue;
          walls[marker].push_back({edge, edge_face.second, z0, roof_z[marker]});
          const auto &v0 = terrain_mesh.vertices[edge.v0];
          const auto &v1 = terrain_mesh.vertices[edge.v1];
          const double horizontal_length =
              Geometry::distance_2d(Vector2D(v0.x, v0.y), Vector2D(v1.x, v1.y));
          if (horizontal_length > Constants::epsilon)
            horizontal_lengths.push_back(horizontal_length);
        }
      }
      const auto block_levels = subdivided_wall_levels(levels, horizontal_lengths);
      for (const auto marker : block)
        wall_levels[marker] = block_levels;
    }

    for (const auto &kv : building_faces)
    {
      const auto marker = kv.first;
      const auto &faces = kv.second;
      const auto roof_height = roof_z[marker];
      const auto &building_walls = walls[marker];
      const auto &levels = wall_levels[marker];

      Mesh building_mesh;
      building_mesh.vertices.reserve(faces.size() * 3 + building_walls.size() * 4 * levels.size());
      building_mesh.faces.reserve(faces.size() + building_walls.size() * 2 * levels.size());
      building_mesh.markers.reserve(faces.size() + building_walls.size() * 2 * levels.size());

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
        building_mesh.faces.push_back(Simplex2D(num_vertices, num_vertices + 1, num_vertices + 2));
        building_mesh.markers.push_back(static_cast<int>(num_buildings + marker));
      }

      // add walls
      for (const auto &wall : building_walls)
      {
        const auto face_center = Geometry::face_center(wall.face, terrain_mesh);
        const auto ground_v0 = terrain_mesh.vertices[wall.edge.v0];
        const auto ground_v1 = terrain_mesh.vertices[wall.edge.v1];
        const auto bottom_v0 = Vector3D(ground_v0.x, ground_v0.y, wall.z0);
        const auto bottom_v1 = Vector3D(ground_v1.x, ground_v1.y, wall.z0);
        const auto top_v1 = Vector3D(ground_v1.x, ground_v1.y, wall.z1);

        const auto wall_normal = Geometry::triangle_normal(bottom_v0, bottom_v1, top_v1);
        const bool flip_orientation =
            Geometry::dot_3d(wall_normal, face_center - bottom_v0) > 0;
        append_vertical_quad_levels(
            building_mesh,
            ground_v0,
            ground_v1,
            wall.z0,
            wall.z1,
            levels,
            static_cast<int>(marker),
            flip_orientation);
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
    terrain_mesh = MeshProcessor::compact_mesh(terrain_mesh);
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

                           double max_triangle_area_size = -1, double min_mesh_angle = 25,
                           const std::string &backend = "auto")
  // Convert 3D Surface to triangle Mesh.
  // - If max_triangle_area_size < 0: uses fast_mesh (earcut/fan triangulation)
  // - If max_triangle_area_size >= 0: uses the Triangle library (requires
  //   DTCC_HAVE_TRIANGLE).
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
      resolve_2d_backend(backend);
      Triangulate::call_triangle(mesh, surface, max_triangle_area_size, min_mesh_angle);
    }
    return mesh;
  }

  static Mesh mesh_multisurface(const MultiSurface &multi_surface,
                                double max_triangle_area_size = -1, double min_mesh_angle = 25,
                                bool weld = false, double snap = 0,
                                const std::string &backend = "auto")
  {
    std::vector<Mesh> multimesh(multi_surface.surfaces.size());
    //    info("meshing multisurface with " + str(multi_surface.surfaces.size())
    //    +
    //         " surfaces");
    // #pragma omp parallel for
    for (size_t i = 0; i < multi_surface.surfaces.size(); i++)
    {
      multimesh[i] =
          mesh_surface(multi_surface.surfaces[i], max_triangle_area_size, min_mesh_angle, backend);
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
                                              double min_mesh_angle = 25, bool weld = false,
                                              const std::string &backend = "auto")
  {
    int n = multi_surfaces.size();
    std::vector<Mesh> meshes(n);
    // Validate the backend before entering the parallel region: an exception
    // escaping an OpenMP loop calls std::terminate instead of propagating.
    if (max_triangle_area_size >= 0)
      resolve_2d_backend(backend);
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
      auto mesh =
          mesh_multisurface(multi_surfaces[i], max_triangle_area_size, min_mesh_angle, weld, 0.0,
                            backend);
      meshes[i] = mesh;
    }

    return meshes;
  }

private:
  enum class TriangulationBackend
  {
    Triangle,
  };

  static std::string normalize_backend_name(const std::string &backend)
  {
    std::string normalized = backend;
    std::transform(
        normalized.begin(), normalized.end(), normalized.begin(),
        [](unsigned char c)
        { return static_cast<char>(std::tolower(c)); });
    return normalized;
  }

  static TriangulationBackend resolve_2d_backend(const std::string &backend)
  {
    const std::string normalized = normalize_backend_name(backend);
    if (normalized.empty() || normalized == "auto")
    {
#ifdef DTCC_HAVE_TRIANGLE
      return TriangulationBackend::Triangle;
#else
      throw std::runtime_error("No triangulation backend is available in this build.");
#endif
    }

    if (normalized == "triangle")
    {
#ifdef DTCC_HAVE_TRIANGLE
      return TriangulationBackend::Triangle;
#else
      throw std::runtime_error(
          "Triangle support not built; reinstall dtcc-core with DTCC_USE_TRIANGLE=ON.");
#endif
    }

    throw std::invalid_argument(
        "Unsupported 2D mesher '" + backend + "'. Expected one of: auto, triangle.");
  }

  static size_t vertical_quad_strip_count(double horizontal_length,
                                          double vertical_height,
                                          double target_aspect_ratio = 15.0)
  {
    if (horizontal_length <= Constants::epsilon || vertical_height <= Constants::epsilon)
      return 1;

    const double target_height = horizontal_length * target_aspect_ratio;
    if (target_height <= Constants::epsilon)
      return 1;

    return std::max<size_t>(
        1, static_cast<size_t>(std::ceil(vertical_height / target_height)));
  }

  // A wall of a building shell: the footprint edge, a face of the building
  // on that edge (to orient the wall outwards) and its bottom and top level.
  struct BuildingWall
  {
    Simplex1D edge;
    Simplex2D face;
    double z0;
    double z1;
  };

  // Roofs closer than this are the same level.
  static double roof_tolerance(double z)
  {
    return 1e-9 * std::max(1.0, std::abs(z));
  }

  // Groups of building markers whose faces share a vertex, each sorted.
  static std::vector<std::vector<size_t>>
  building_blocks(const std::map<size_t, std::vector<Simplex2D>> &building_faces)
  {
    std::map<size_t, size_t> parent;
    for (const auto &kv : building_faces)
      parent[kv.first] = kv.first;
    const auto find = [&parent](size_t marker)
    {
      size_t root = marker;
      while (parent[root] != root)
        root = parent[root];
      while (parent[marker] != root)
      {
        const size_t next = parent[marker];
        parent[marker] = root;
        marker = next;
      }
      return root;
    };

    std::unordered_map<size_t, size_t> vertex_owner;
    for (const auto &kv : building_faces)
    {
      for (const auto &face : kv.second)
      {
        for (int i = 0; i < 3; ++i)
        {
          const auto inserted = vertex_owner.emplace(face[i], kv.first);
          if (inserted.second)
            continue;
          const size_t a = find(inserted.first->second);
          const size_t b = find(kv.first);
          if (a != b)
            parent[std::max(a, b)] = std::min(a, b);
        }
      }
    }

    std::map<size_t, std::vector<size_t>> groups;
    for (const auto &kv : building_faces)
      groups[find(kv.first)].push_back(kv.first);
    std::vector<std::vector<size_t>> blocks;
    blocks.reserve(groups.size());
    for (auto &kv : groups)
      blocks.push_back(std::move(kv.second));
    return blocks;
  }

  // The levels at which the walls of one block are split: the given levels
  // (base and roofs), with every interval between consecutive levels split
  // into strips. The shortest wall edge sets one strip height for the block,
  // so the most slender wall panels are subdivided by construction (TetGen
  // quality is dominated by preserved wall panels).
  static std::vector<double>
  subdivided_wall_levels(std::vector<double> levels,
                         const std::vector<double> &horizontal_lengths)
  {
    std::sort(levels.begin(), levels.end());
    levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
    const double shortest =
        horizontal_lengths.empty()
            ? 0.0
            : *std::min_element(horizontal_lengths.begin(), horizontal_lengths.end());

    std::vector<double> subdivided;
    if (levels.empty())
      return subdivided;
    subdivided.push_back(levels.front());
    for (size_t i = 1; i < levels.size(); ++i)
    {
      const double lo = levels[i - 1];
      const double hi = levels[i];
      const size_t strips =
          horizontal_lengths.empty() ? 1 : vertical_quad_strip_count(shortest, hi - lo, 5.0);
      for (size_t step = 1; step < strips; ++step)
      {
        const double t = static_cast<double>(step) / static_cast<double>(strips);
        subdivided.push_back(lo + t * (hi - lo));
      }
      subdivided.push_back(hi);
    }
    return subdivided;
  }

  // Append a vertical wall over the footprint edge (v0, v1) from level z0 to
  // level z1, as one quad (two triangles) between each pair of consecutive
  // levels in that range.
  static void append_vertical_quad_levels(Mesh &mesh,
                                          const Vector3D &v0,
                                          const Vector3D &v1,
                                          double z0,
                                          double z1,
                                          const std::vector<double> &levels,
                                          int marker,
                                          bool flip_orientation)
  {
    std::vector<double> zs;
    for (const double z : levels)
    {
      if (z >= std::min(z0, z1) && z <= std::max(z0, z1))
        zs.push_back(z);
    }
    if (z0 > z1)
      std::reverse(zs.begin(), zs.end());

    for (size_t step = 0; step + 1 < zs.size(); ++step)
    {
      const auto base = mesh.vertices.size();
      mesh.vertices.push_back(Vector3D(v0.x, v0.y, zs[step]));
      mesh.vertices.push_back(Vector3D(v1.x, v1.y, zs[step]));
      mesh.vertices.push_back(Vector3D(v0.x, v0.y, zs[step + 1]));
      mesh.vertices.push_back(Vector3D(v1.x, v1.y, zs[step + 1]));

      if (flip_orientation)
      {
        mesh.faces.push_back(Simplex2D(base, base + 3, base + 1));
        mesh.faces.push_back(Simplex2D(base, base + 2, base + 3));
      }
      else
      {
        mesh.faces.push_back(Simplex2D(base, base + 1, base + 3));
        mesh.faces.push_back(Simplex2D(base, base + 3, base + 2));
      }

      mesh.markers.push_back(marker);
      mesh.markers.push_back(marker);
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
