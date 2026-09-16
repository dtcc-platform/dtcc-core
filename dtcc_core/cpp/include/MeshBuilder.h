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
    for (const auto &kv : building_faces)
    {
      const auto marker = kv.first;
      const auto &faces = kv.second;
      const auto &building = buildings[marker];
      const auto roof_height = building.max_height();

      auto naked_edges = MeshProcessor::find_naked_edges(faces);
      const auto wall_slices = building_wall_strip_count(naked_edges, terrain_mesh, roof_height);

      Mesh building_mesh;
      building_mesh.vertices.reserve(faces.size() * 3 + naked_edges.size() * 4 * wall_slices);
      building_mesh.faces.reserve(faces.size() + naked_edges.size() * 2 * wall_slices);
      building_mesh.markers.reserve(faces.size() + naked_edges.size() * 2 * wall_slices);

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
        const bool flip_orientation =
            Geometry::dot_3d(wall_normal, face_center - ground_v0) > 0;
        building_mesh.vertices.resize(base);
        append_vertical_quad_strips(
            building_mesh,
            ground_v0,
            ground_v1,
            roof_v0,
            roof_v1,
            static_cast<int>(marker),
            flip_orientation,
            wall_slices);
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

  static size_t building_wall_strip_count(
      const std::vector<std::pair<Simplex1D, Simplex2D>> &naked_edges,
      const Mesh &terrain_mesh,
      double roof_height)
  {
    std::vector<double> horizontal_lengths;
    horizontal_lengths.reserve(naked_edges.size());
    double max_wall_height = 0.0;

    for (const auto &edge_faces : naked_edges)
    {
      const Simplex1D edge = edge_faces.first;
      const auto &ground_v0 = terrain_mesh.vertices[edge.v0];
      const auto &ground_v1 = terrain_mesh.vertices[edge.v1];
      const double horizontal_length = Geometry::distance_2d(
          Vector2D(ground_v0.x, ground_v0.y), Vector2D(ground_v1.x, ground_v1.y));
      if (horizontal_length > Constants::epsilon)
        horizontal_lengths.push_back(horizontal_length);

      max_wall_height = std::max(max_wall_height, std::abs(roof_height - ground_v0.z));
      max_wall_height = std::max(max_wall_height, std::abs(roof_height - ground_v1.z));
    }

    if (horizontal_lengths.empty())
      return 1;

    std::sort(horizontal_lengths.begin(), horizontal_lengths.end());
    // Use the shortest supported wall edge for one conforming strip count per
    // building shell. This keeps the wall mesh watertight while ensuring that
    // the most slender preserved wall panels are subdivided by construction.
    // A tighter target than the historical default is deliberate here: TetGen
    // quality is dominated by preserved wall panels, so the shell builder must
    // hand over a genuinely well-proportioned wall mesh instead of relying on
    // a later repair pass.
    const double shortest_horizontal_length = horizontal_lengths.front();
    return vertical_quad_strip_count(
        shortest_horizontal_length, max_wall_height, 5.0);
  }

  static void append_vertical_quad_strips(Mesh &mesh,
                                          const Vector3D &bottom_v0,
                                          const Vector3D &bottom_v1,
                                          const Vector3D &top_v0,
                                          const Vector3D &top_v1,
                                          int marker,
                                          bool flip_orientation,
                                          size_t strip_count)
  {
    const size_t strips = std::max<size_t>(1, strip_count);
    mesh.vertices.reserve(mesh.vertices.size() + 4 * strips);
    mesh.faces.reserve(mesh.faces.size() + 2 * strips);
    mesh.markers.reserve(mesh.markers.size() + 2 * strips);

    const auto interpolate = [](const Vector3D &a, const Vector3D &b, double t)
    {
      if (t <= 0.0)
        return a;
      if (t >= 1.0)
        return b;
      return Vector3D(a.x + t * (b.x - a.x),
                      a.y + t * (b.y - a.y),
                      a.z + t * (b.z - a.z));
    };

    for (size_t step = 0; step < strips; ++step)
    {
      const double t0 = static_cast<double>(step) / static_cast<double>(strips);
      const double t1 = static_cast<double>(step + 1) / static_cast<double>(strips);

      const auto lower_v0 = interpolate(bottom_v0, top_v0, t0);
      const auto lower_v1 = interpolate(bottom_v1, top_v1, t0);
      const auto upper_v0 = interpolate(bottom_v0, top_v0, t1);
      const auto upper_v1 = interpolate(bottom_v1, top_v1, t1);

      const auto base = mesh.vertices.size();
      mesh.vertices.push_back(lower_v0);
      mesh.vertices.push_back(lower_v1);
      mesh.vertices.push_back(upper_v0);
      mesh.vertices.push_back(upper_v1);

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
