#ifndef DTCC_CORE_SPADE_TRIANGULATE_HPP
#define DTCC_CORE_SPADE_TRIANGULATE_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <spade_wrapper.h>

#include "BoundingBox.h"
#include "BoundingBoxTree.h"
#include "Geometry.h"
#include "Logging.h"
#include "Timer.h"
#include "model/Mesh.h"
#include "model/Polygon.h"
#include "model/Vector.h"

namespace DTCC_BUILDER
{

// Simple smoke test to verify SPADE integration is available at runtime.
inline void spade_wrapper_smoke_test()
{
  std::vector<spade::Point> outer = {
      {0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 0.0}};

  std::vector<std::vector<spade::Point>> holes;

  auto result = spade::triangulate(outer, holes, 0.0, spade::Quality::Default, true);

  if (result.points.size() < 3 || result.triangles.empty())
  {
    throw std::runtime_error("SPADE triangulation returned an unexpected result");
  }
}

namespace detail
{
inline void compute_spade_domain_markers(Mesh &mesh, const std::vector<Polygon> &subdomains)
{
  info("Computing domain markers...");
  Timer timer("compute_domain_markers");

  BoundingBoxTree2D search_tree;
  std::vector<BoundingBox2D> bounding_boxes;
  bounding_boxes.reserve(subdomains.size());
  for (const auto &subdomain : subdomains)
  {
    bounding_boxes.emplace_back(subdomain);
  }
  search_tree.build(bounding_boxes);

  mesh.markers.resize(mesh.faces.size());
  std::fill(mesh.markers.begin(), mesh.markers.end(), -2);

  std::vector<bool> is_building_vertex(mesh.vertices.size());
  std::fill(is_building_vertex.begin(), is_building_vertex.end(), false);

  if (!subdomains.empty())
  {
    for (size_t i = 0; i < mesh.faces.size(); i++)
    {
      const Vector3D c_3d = mesh.mid_point(i);
      const Vector2D c_2d(c_3d.x, c_3d.y);
      std::vector<size_t> indices = search_tree.find(Vector2D(c_2d));

      if (!indices.empty())
      {
        for (const auto &index : indices)
        {
          if (Geometry::polygon_contains_2d(subdomains[index], c_2d))
          {
            mesh.markers[i] = index;
            const Simplex2D &T = mesh.faces[i];
            is_building_vertex[T.v0] = true;
            is_building_vertex[T.v1] = true;
            is_building_vertex[T.v2] = true;
            break;
          }
        }
      }
    }

    for (size_t i = 0; i < mesh.faces.size(); i++)
    {
      const Simplex2D &T = mesh.faces[i];
      const bool touches_building =
          (is_building_vertex[T.v0] || is_building_vertex[T.v1] || is_building_vertex[T.v2]);

      if (touches_building && mesh.markers[i] == -2)
        mesh.markers[i] = -1;
    }
  }
}
} // namespace detail

inline Mesh spade_build_ground_mesh(const std::vector<Polygon> &subdomains,
                                    const std::vector<Polygon> &holes,
                                    const std::vector<double> &subdomain_triangle_size,
                                    double xmin,
                                    double ymin,
                                    double xmax,
                                    double ymax,
                                    double max_mesh_size,
                                    double min_mesh_angle,
                                    bool sort_triangles = false)
{

  info("SPADE Triangulation: Building ground mesh...");
  Timer timer("build_ground_mesh");

  const BoundingBox2D bounding_box(Vector2D(xmin, ymin), Vector2D(xmax, ymax));
  const size_t nx =
      static_cast<size_t>((bounding_box.Q.x - bounding_box.P.x) / max_mesh_size);
  const size_t ny =
      static_cast<size_t>((bounding_box.Q.y - bounding_box.P.y) / max_mesh_size);
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

  auto make_loop = [](const std::vector<Vector2D> &polygon) -> std::vector<spade::Point> {
    std::vector<spade::Point> loop;
    if (polygon.empty())
      return loop;
    loop.reserve(polygon.size() + 1);
    for (const auto &p : polygon)
      loop.push_back({p.x, p.y, 0.0});
    const Vector2D &first = polygon.front();
    const Vector2D &last = polygon.back();
    if (std::abs(first.x - last.x) > 1e-10 || std::abs(first.y - last.y) > 1e-10)
      loop.push_back({first.x, first.y, 0.0});
    return loop;
  };

  std::vector<spade::Point> outer_loop = make_loop(boundary);
  std::vector<std::vector<spade::Point>> inner_loops;
  inner_loops.reserve(triangle_sub_domains.size());
  for (const auto &polygon : triangle_sub_domains)
  {
    auto loop = make_loop(polygon);
    if (loop.size() >= 3)
      inner_loops.push_back(std::move(loop));
  }

  std::vector<std::vector<spade::Point>> hole_loops;
  hole_loops.reserve(triangle_holes.size());
  for (const auto &polygon : triangle_holes)
  {
    auto loop = make_loop(polygon);
    if (loop.size() >= 3)
      hole_loops.push_back(std::move(loop));
  }

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

  spade::Quality quality =
      (min_mesh_angle >= 25.0) ? spade::Quality::Moderate : spade::Quality::Default;

  auto result = spade::triangulate(outer_loop, hole_loops, inner_loops, effective_maxh, quality, true);

  if (result.points.size() < 3 || result.triangles.empty())
  {
    throw std::runtime_error("SPADE triangulation returned an unexpected result");
  }

  Mesh mesh;
  mesh.vertices.reserve(result.points.size());
  for (const auto &p : result.points)
  {
    mesh.vertices.emplace_back(p.x, p.y, p.z);
  }

  mesh.faces.reserve(result.triangles.size());
  for (const auto &tri : result.triangles)
  {
    mesh.faces.emplace_back(tri.v0, tri.v1, tri.v2, sort_triangles);
  }

  detail::compute_spade_domain_markers(mesh, subdomains);

  return mesh;
}

} // namespace DTCC_BUILDER

#endif // DTCC_CORE_SPADE_TRIANGULATE_HPP
