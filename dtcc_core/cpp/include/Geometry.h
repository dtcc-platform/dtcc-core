// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_GEOMETRY_H
#define DTCC_GEOMETRY_H

#include <algorithm>
#include <cmath>
#include <iso646.h>
#include <limits>
#include <stack>
#include <tuple>
#include <vector>

#include "Eigen/Eigen"
#include "Eigen/Geometry"

#include "BoundingBox.h"
#include "Constants.h"
#include "model/Mesh.h"
#include "model/Polygon.h"
#include "model/Simplices.h"
#include "model/Surface.h"
#include "model/Vector.h"
#include "model/VolumeMesh.h"

namespace DTCC_BUILDER
{

class Geometry
{
public:
  // FIXME: This needs to be reworked in light of the introduction of two
  // different classes Point and Vector. Currently somewhat inconsistent.

  // Compute squared norm (2D)
  static double squared_norm_2d(const Vector2D &v) { return dot_2d(v, v); }

  // Compute squared norm (3D)
  static double squared_norm_3d(const Vector3D &v) { return dot_3d(v, v); }

  // Compute norm (2D)
  static double norm_2d(const Vector2D &v) { return std::sqrt(squared_norm_2d(v)); }

  // Compute norm (3D)
  static double norm_3d(const Vector3D &v) { return std::sqrt(squared_norm_3d(v)); }

  // Compute dot product (2D)
  static double dot_2d(const Vector2D &u, const Vector2D &v) { return u.x * v.x + u.y * v.y; }

  // Compute dot product (3D)
  static double dot_3d(const Vector3D &u, const Vector3D &v)
  {
    return u.x * v.x + u.y * v.y + u.z * v.z;
  }

  // Compute cross product (3D)
  static Vector3D cross_3d(const Vector3D &u, const Vector3D &v)
  {
    return Vector3D(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);
  }

  // Compute distance between points (2D)
  static double distance_2d(const Vector2D &p, const Vector2D &q)
  {
    return std::sqrt(squared_distance_2d(p, q));
  }

  // Compute distance between segment (p0, p1) and point q (2D)
  static double distance_2d(const Vector2D &p0, const Vector2D &p1, const Vector2D &q)
  {
    return std::sqrt(squared_distance_2d(p0, p1, q));
  }

  // Compute distance between polygon and point (2D)
  static double distance_2d(const Polygon &polygon, const Vector2D &p)
  {
    return std::sqrt(squared_distance_2d(polygon, p));
  }

  // Compute distance between polygons (2D)
  static double distance_2d(const Polygon &polygon0, const Polygon &polygon1)
  {
    return std::sqrt(squared_distance_2d(polygon0, polygon1));
  }

  // Compute distance between points (3D)
  static double distance_3d(const Vector3D &p, const Vector3D &q)
  {
    return std::sqrt(squared_distance_3d(p, q));
  }

  // Compute squared distance between points (2D)
  static double squared_distance_2d(const Vector2D &p, const Vector2D &q)
  {
    const double dx = p.x - q.x;
    const double dy = p.y - q.y;
    return dx * dx + dy * dy;
  }

  // Compute squared distance between segment (p0, p1) and point q (2D)
  static double squared_distance_2d(const Vector2D &p0, const Vector2D &p1, const Vector2D &q)
  {
    // Project point to line
    const Vector2D u(p0, q);
    const Vector2D v(p0, p1);
    const Vector2D p = p0 + v * (dot_2d(u, v) / v.squared_magnitude());

    // Check whether projected point is inside segment. Check either
    // x or y coordinates depending on which is largest (most stable)
    const bool inside = std::abs(v.x) > std::abs(v.y)
                            ? std::min(p0.x, p1.x) <= p.x && p.x <= std::max(p0.x, p1.x)
                            : std::min(p0.y, p1.y) <= p.y && p.y <= std::max(p0.y, p1.y);

    // Use distance to projection if inside
    if (inside)
      return squared_distance_2d(p, q);

    // Otherwise use distance to closest end point
    const double d0 = squared_distance_2d(p0, q);
    const double d1 = squared_distance_2d(p1, q);
    return std::min(d0, d1);
  }

  // Compute squared distance between polygon and point(2D)
  static double squared_distance_2d(const Polygon &polygon, const Vector2D &p)
  {
    // Check if point is contained in polygon
    if (polygon_contains_2d(polygon, p))
      return 0.0;

    // If not, compute minimal squared distance to all segments
    double d_to_min = std::numeric_limits<double>::max();
    for (size_t i = 0; i < polygon.vertices.size(); i++)
    {
      Vector2D p0 = polygon.vertices[i];
      Vector2D p1 = polygon.vertices[(i + 1) % polygon.vertices.size()];
      d_to_min = std::min(d_to_min, squared_distance_2d(p0, p1, p));
    }

    return d_to_min;
  }

  // Compute squared distance between polygons (2D)
  static double squared_distance_2d(const Polygon &polygon0, const Polygon &polygon1)
  {
    //    std::cout << "Checking polygon with " << polygon0.vertices.size()
    //              << " vertices" << std::endl;
    //    std::cout << "Checking polygon with " << polygon1.vertices.size()
    //             << " vertices" << std::endl;
    double d_to_min = std::numeric_limits<double>::max();

    // Check all vertices in first polygon
    for (auto const &p : polygon0.vertices)
    {
      //     std::cout << "Point of polygon0 " << str(p) << std::endl;
      d_to_min = std::min(d_to_min, squared_distance_2d(polygon1, p));
      //     std::cout << "D2min " << d_to_min << std::endl;
    }
    // Check all vertices in second polygon
    for (auto const &p : polygon1.vertices)
    {
      //      std::cout << "Point of polygon1 " << str(p) << std::endl;
      d_to_min = std::min(d_to_min, squared_distance_2d(polygon0, p));
      //      std::cout << "D2min " << d_to_min << std::endl;
    }
    return d_to_min;
  }

  // Compute squared distance between points (3D)
  static double squared_distance_3d(const Vector3D &p, const Vector3D &q)
  {
    const double dx = p.x - q.x;
    const double dy = p.y - q.y;
    const double dz = p.z - q.z;
    return dx * dx + dy * dy + dz * dz;
  }

  // Compute orientation of point q relative to edge (p0, p1) (2D)
  static double orient_2d(const Vector2D &p0, const Vector2D &p1, const Vector2D &q)
  {
    const Vector2D u(p0, p1);
    const Vector2D v(p0, q);
    return u.x * v.y - u.y * v.x;
  }

  // Compute sign of point q relative to edge (p0, p1) (2D).
  // -1 --- (p0) --- 0 --- (p1) --- +1
  static int edge_sign_2d(const Vector2D &p0, const Vector2D &p1, const Vector2D &q)
  {
    double l{}, d0{}, d1{};
    if (std::abs(p0.x - p1.x) > std::abs(p0.y - p1.y))
    {
      l = std::abs(p0.x - p1.x);
      d0 = std::abs(p0.x - q.x);
      d1 = std::abs(p1.x - q.x);
    }
    else
    {
      l = std::abs(p0.y - p1.y);
      d0 = std::abs(p0.y - q.y);
      d1 = std::abs(p1.y - q.y);
    }
    if (d0 > l - Constants::epsilon && d0 > d1)
      return 1;
    else if (d1 > l - Constants::epsilon && d1 > d0)
      return -1;
    else
      return 0;
  }

  // Compute strictly increasing function [-pi, pi] -> [-2, 2] of angle of v
  // relative to u. This is a cheap alternative compared to working with asin,
  // acos.
  static double vector_angle_2d(const Vector2D &u, const Vector2D &v)
  {
    const double u2 = u.x * u.x + u.y * u.y;
    const double v2 = v.x * v.x + v.y * v.y;
    const double sin = u.x * v.y - u.y * v.x;
    const double cos = u.x * v.x + u.y * v.y;
    const double a = sin * sin / (u2 * v2);
    return (sin > 0.0 ? (cos > 0.0 ? a : 2.0 - a) : (cos > 0.0 ? -a : a - 2.0));
  }

  // Compute quadrant angle of point p relative to polygon (2D)
  static int quadrant_angle_2d(const Vector2D &p, const std::vector<Vector2D> &polygon)
  {
    // Compute angle to first vertex
    Vector2D q0 = polygon[0];
    int v0 = quadrant_angle_2d(q0, p);

    // Sum up total angle
    int total_angle = 0;
    for (size_t i = 1; i < polygon.size() + 1; i++)
    {
      // Compute angle increment
      Vector2D q1 = polygon[i % polygon.size()];
      int v1 = quadrant_angle_2d(q1, p);
      int dv = v1 - v0;

      // Adjust angle increment for wrap-around
      if (dv == 3)
        dv = -1;
      else if (dv == -3)
        dv = 1;
      else if (dv == 2 || dv == -2)
      {
        double xx = q1.x - ((q1.y - p.y) * ((q0.x - q1.x) / (q0.y - q1.y)));
        if (xx > p.x)
          dv = -dv;
      }

      // Add to total angle and update
      total_angle += dv;
      q0 = q1;
      v0 = v1;
    }

    return total_angle;
  }

  // Compute quadrant angle of point p relative to point q (2D)
  static int quadrant_angle_2d(const Vector2D &p, const Vector2D &q)
  {
    return ((p.x > q.x) ? ((p.y > q.y) ? 0 : 3) : ((p.y > q.y) ? 1 : 2));
  }

  // Compute face normal
  static Vector3D face_normal_3d(const Simplex2D &face, const VolumeMesh &mesh_3d)
  {
    const Vector3D p0{mesh_3d.vertices[face.v0]};
    const Vector3D p1{mesh_3d.vertices[face.v1]};
    const Vector3D p2{mesh_3d.vertices[face.v2]};
    return triangle_normal(p0, p1, p2);
  }

  static Vector3D triangle_normal(const Vector3D &p0, const Vector3D &p1, const Vector3D &p2)
  {
    const Vector3D u = p1 - p0;
    const Vector3D v = p2 - p0;
    Vector3D n = cross_3d(u, v);
    n /= Geometry::norm_3d(n);
    return n;
  }

  static Vector3D surface_normal(const Surface &surface)
  {
    const Vector3D p0{surface.vertices[0]};
    const Vector3D p1{surface.vertices[1]};
    const Vector3D p2{surface.vertices[2]};
    return triangle_normal(p0, p1, p2);
  }

  static Vector3D surface_centroid(const Surface &surface)
  {
    Vector3D c{};
    for (auto const &p : surface.vertices)
      c += p;
    c /= static_cast<double>(surface.vertices.size());
    return c;
  }

  static bool is_convex(const Polygon &polygon)
  {
    int orientation = polygon_orientation_2d(polygon);
    if (orientation == 0)
      orientation = -1;
    size_t num_vertices = polygon.vertices.size();
    for (size_t i = 0; i < num_vertices; i++)
    {
      Vector2D p0 = polygon.vertices[i];
      Vector2D p1 = polygon.vertices[(i + 1) % num_vertices];
      Vector2D p2 = polygon.vertices[(i + 2) % num_vertices];
      Vector2D v0(p0, p1);
      Vector2D v1(p1, p2);
      double cross = v0.x * v1.y - v0.y * v1.x;
      if (cross * orientation < 0)
        return false;
    }
    return true;
  }

  static Polygon project_surface(const Surface &surface)
  {
    const auto z_normal = Eigen::Vector3d(0, 0, 1);
    auto normal = Geometry::surface_normal(surface);
    auto e_norm = Eigen::Vector3d(normal.x, normal.y, normal.z);
    auto transform = Eigen::Transform<double, 3, Eigen::Isometry>();
    transform = Eigen::Quaterniond::FromTwoVectors(e_norm, z_normal);
    Polygon projected_polygon;
    for (const auto &v : surface.vertices)
    {
      auto e_v = Eigen::Vector3d(v.x, v.y, v.z);
      auto e_v_prime = transform * e_v;
      projected_polygon.vertices.push_back(Vector2D(e_v_prime.x(), e_v_prime.y()));
    }
    return projected_polygon;
  }
  static bool is_convex(const Surface &surface)
  {
    auto poly = project_surface(surface);
    return is_convex(poly);
  }

  static Vector3D face_normal(const Simplex2D &face, const Mesh &mesh)
  {
    const Vector3D p0{mesh.vertices[face.v0]};
    const Vector3D p1{mesh.vertices[face.v1]};
    const Vector3D p2{mesh.vertices[face.v2]};
    return triangle_normal(p0, p1, p2);
  }

  // Compute face center
  static Vector3D face_center(const Simplex2D &face, const Mesh &mesh)
  {
    Vector3D c = mesh.vertices[face.v0];
    c += mesh.vertices[face.v1];
    c += mesh.vertices[face.v2];
    c /= 3.0;
    return c;
  }

  // Compute cell center
  static Vector3D cell_center_3d(const Simplex3D &cell, const VolumeMesh &mesh_3d)
  {
    Vector3D c{};
    c += Vector3D(mesh_3d.vertices[cell.v0]);
    c += Vector3D(mesh_3d.vertices[cell.v1]);
    c += Vector3D(mesh_3d.vertices[cell.v2]);
    c += Vector3D(mesh_3d.vertices[cell.v3]);
    c /= 4.0;
    return c;
  }

  // Compute signed determinant of polygon (2D)
  static double polygon_determinant_2d(const Polygon &polygon)
  {
    double sum = 0.0;
    for (size_t i = 0; i < polygon.vertices.size(); i++)
    {
      Vector2D p0 = polygon.vertices[i];
      Vector2D p1 = polygon.vertices[(i + 1) % polygon.vertices.size()];
      sum += (p1.x - p0.x) * (p1.y + p0.y);
    }
    return sum;
  }

  // Compute orientation of polygon (0 = counter-clockwise, 1 = clockwise)
  static size_t polygon_orientation_2d(const Polygon &polygon)
  {
    return polygon_determinant_2d(polygon) < 0 ? 0 : 1;
  }

  // Compute area of polygon (2D)
  static double polygon_area(const Polygon &polygon)
  {
    return 0.5 * std::abs(polygon_determinant_2d(polygon));
  }

  static bool self_intersects(const Polygon &polygon)
  {
    for (size_t i = 0; i < polygon.vertices.size(); i++)
    {
      Vector2D p0 = polygon.vertices[i];
      Vector2D p1 = polygon.vertices[(i + 1) % polygon.vertices.size()];
      for (size_t j = i + 2; j < polygon.vertices.size(); j++)
      {
        Vector2D q0 = polygon.vertices[j];
        Vector2D q1 = polygon.vertices[(j + 1) % polygon.vertices.size()];
        if (intersects_2d(p0, p1, q0, q1, true))
          return true;
      }
    }
    return false;
  }

  static double surface_area(const Surface &surface)
  {
    return polygon_area(project_surface(surface));
  }

  // Return a 'random' point inside the polygon
  static Vector2D point_inside_polygon_2d(const Polygon &polygon)
  {
    bool found = false;
    Vector2D pc;
    for (size_t idx = 0; idx < polygon.vertices.size() - 1; idx++)
    {
      auto v0 = polygon.vertices[idx];
      auto v1 = polygon.vertices[idx + 1];
      auto v2 = polygon.vertices[(idx + 2) % polygon.vertices.size()];
      auto u = v0 - v1;
      auto v = v2 - v1;
      auto c = u.x * v.y - u.y * v.x;
      if (abs(c) < 1e-6)
        continue; // colinear
      pc = (v0 + v1 + v2) / 3;
      if (polygon_contains_2d(polygon, pc))
      {
        found = true;
        break;
      }
    }
    if (!found)
      error("failed to find point inside polygon");
    return pc;
  }

  // Compute area of a triangle (3D)
  static double triangle_area(const Vector3D &p0, const Vector3D &p1, const Vector3D &p2)
  {
    double side_a = (p0 - p1).magnitude();
    double side_b = (p1 - p2).magnitude();
    double side_c = (p2 - p0).magnitude();

    double s = 0.5 * (side_a + side_b + side_c);
    return std::sqrt(s * (s - side_a) * (s - side_b) * (s - side_c));
  }

  // Computer Perimeter of polygon (2D)
  static double polygon_perimeter_2d(const Polygon &polygon)
  {
    double sum = 0.0;
    for (size_t i = 0; i < polygon.vertices.size(); i++)
    {
      Vector2D p0 = polygon.vertices[i];
      Vector2D p1 = polygon.vertices[(i + 1) % polygon.vertices.size()];
      sum += Geometry::distance_2d(p0, p1);
    }
    return sum;
  }

  // Compute center of polygon (2D)
  static Vector2D polygon_center_2d(const Polygon &polygon)
  {
    Vector2D o{};
    Vector2D c{};
    for (auto const &p : polygon.vertices)
      c += p;
    c /= static_cast<double>(polygon.vertices.size());
    return c;
  }

  // Compute radius of polygon relative to center (2D)
  static double polygon_radius_2d(const Polygon &polygon, const Vector2D &center)
  {
    double r_to_max = 0.0;
    for (auto const &p : polygon.vertices)
    {
      const double r2 = squared_distance_2d(p, center);
      if (r2 > r_to_max)
        r_to_max = r2;
    }
    return std::sqrt(r_to_max);
  }

  // FIXME: Cleanup Contains vs Intersects vs Collide. What should we name it?

  // Check whether edge (p0, p1) contains point q. It is assumed that the
  // point is located on the line defined by the edge.
  static bool edge_contains_2d(const Vector2D &p0, const Vector2D &p1, const Vector2D &q,
                               double tol = 0.0)
  {
    const Vector2D v(p0, p1);
    if (std::abs(v.x) > std::abs(v.y))
      return std::min(p0.x, p1.x) - tol < q.x and std::max(p0.x, p1.x) + tol > q.x;
    else
      return std::min(p0.y, p1.y) - tol < q.y and std::max(p0.y, p1.y) + tol > q.y;
  }

  // Check whether polygon contains point (2D)
  static bool polygon_contains_2d(const Polygon &polygon, const Vector2D &p)
  {
    // Compute total quadrant relative to polygon. If the point
    // is inside the polygon, the angle should be 4 (or -4).
    bool inside = Geometry::quadrant_angle_2d(p, polygon.vertices) != 0;

    if (inside && polygon.holes.size() > 0) // inside shell, check holes
    {
      for (auto const &hole : polygon.holes)
      {
        bool inside_hole = Geometry::quadrant_angle_2d(p, hole) != 0;
        if (inside_hole)
          return false;
      }
    }

    return inside;
  }

  // Check whether bounding box contains point (2D)
  static bool bounding_box_contains_2d(const BoundingBox2D &bbox, const Vector2D &p,
                                       double margin = 0.0)
  {
    return (bbox.P.x + margin <= p.x && p.x + margin <= bbox.Q.x && bbox.P.y + margin <= p.y &&
            p.y + margin <= bbox.Q.y);
  }

  // Check whether bounding box contains point (3D)
  static bool bounding_box_contains_3d(const BoundingBox3D &bbox, const Vector3D &p,
                                       double margin = 0.0)
  {
    return (bbox.P.x + margin <= p.x && p.x + margin <= bbox.Q.x && bbox.P.y + margin <= p.y &&
            p.y + margin <= bbox.Q.y && bbox.P.z + margin <= p.z && p.z + margin <= bbox.Q.z);
  }

  // Check whether bounding box contains polygon (2D)
  static bool bounding_box_contains_2d(const BoundingBox2D &bbox, const Polygon &polygon,
                                       double margin = 0.0)
  {
    for (const auto &p : polygon.vertices)
      if (!bounding_box_contains_2d(bbox, p, margin))
        return false;
    return true;
  }

  // Check whether edges (p0, p1) and (q0, q1) intersect
  static bool intersects_2d(const Vector2D &p0, const Vector2D &p1, const Vector2D &q0,
                            const Vector2D &q1, bool strict = false)

  {

    bool crosses = (orient_2d(p0, p1, q0) * orient_2d(p0, p1, q1) <= 0.0 &&
                    orient_2d(q0, q1, p0) * orient_2d(q0, q1, p1) <= 0.0);
    if (crosses && strict)
    {
      if (p0 == q0 || p0 == q1 || p1 == q0 || p1 == q1)
        return false;
    }
    return crosses;
  }

  // Check whether bounding boxes intersect (2D)
  static bool intersect_2d(const BoundingBox2D &bbox_a, const BoundingBox2D &bbox_b)
  {
    return (bbox_a.P.x <= bbox_b.Q.x && bbox_b.P.x <= bbox_a.Q.x && bbox_a.P.y <= bbox_b.Q.y &&
            bbox_b.P.y <= bbox_a.Q.y);
  }

  // Check whether polygon intersects with polygon (2D)
  static bool intersects_2d(const Polygon &polygon_a, const Polygon &polygon_b)
  {
    // Check if bounding boxes intersect
    if (!intersect_2d(BoundingBox2D(polygon_a), BoundingBox2D(polygon_b)))
      return false;

    // Check if any edge of polygon_a intersects with any edge of polygon_b
    for (const auto &p0 : polygon_a.vertices)
      for (const auto &p1 : polygon_a.vertices)
        for (const auto &q0 : polygon_b.vertices)
          for (const auto &q1 : polygon_b.vertices)
            if (intersects_2d(p0, p1, q0, q1))
              return true;

    // Check if one polygon contains the other.
    if (polygon_contains_2d(polygon_a, polygon_b.vertices[0]) ||
        polygon_contains_2d(polygon_b, polygon_a.vertices[0]))
      return true;

    return false;
  }

  // Compute intersection between edges p0 - p1 and q0 - q1 (2D)
  static Vector2D edge_intersection_2d(const Vector2D &p0, const Vector2D &p1, const Vector2D &q0,
                                       const Vector2D &q1)
  {
    // Solve for intersection: p0 + k*(p1 - p0) = q0 + l*(q1 - q0)

    // Compute vectors
    const Vector2D u(p0, q0);
    const Vector2D v(p0, p1);
    const Vector2D w(q0, q1);

    // Create linear system
    const double a = v.x;
    const double b = -w.x;
    const double c = v.y;
    const double d = -w.y;
    const double e = u.x;
    const double f = u.y;

    // Compute determinant
    const double det = a * d - b * c;

    // Check if close to parallel
    // if (std::(det) < parallelTolerance)
    // {
    //     throw std::runtime_error("Segments are parallel.");
    // }

    // Solve linear system
    const double k = (d * e - b * f) / det;
    Vector2D p = p0 + v * k;

    return p;
  }

  // Compute convex hull of point set (2D)
  static Polygon convex_hull_2d(const std::vector<Vector2D> &points)
  {
    // The convex hull is computed by doing a Graham scan: select an
    // extreme base point, sort remaining points by angle and then
    // add points that create a left turn around the perimeter.

    // find point with smallest y-coordinate. If y-coordinate is
    // the same, sort by smallest x-coordinate.
    double x_min = points[0].x;
    double y_min = points[0].y;
    size_t i_min = 0;
    const size_t num_points = points.size();
    for (size_t i = 1; i < num_points; i++)
    {
      const double x = points[i].x;
      const double y = points[i].y;
      if (y < y_min || (y == y_min && x < x_min))
      {
        x_min = x;
        y_min = y;
        i_min = i;
      }
    }

    // Set base point
    const size_t base_index = i_min;
    Vector2D base_point = points[base_index];

    // Compute angles and distances relative to base point
    std::vector<std::tuple<double, double, size_t>> angles(num_points - 1);
    size_t k = 0;
    for (size_t i = 0; i < num_points; i++)
    {
      // Skip base point
      if (i == base_index)
        continue;

      // Compute angle (negative cosine) and distance
      const Vector2D &p = points[i];
      const Vector2D v(base_point, p);
      const double distance = v.magnitude();
      const double angle = (distance > Constants::epsilon ? -v.x / distance : 0.0);

      // Store angle and distance along with index (for sorting)
      angles[k++] = std::make_tuple(angle, distance, i);
    }

    // Sort by angles (primary) and distance (secondary) to base point
    std::sort(angles.begin(), angles.end());

    // Filter out points with unique angles, keeping only furthest point
    std::vector<size_t> filtered_indices;
    double last_angle = 2.0; // no angle has this value
    for (size_t i = 0; i < num_points - 1; i++)
    {
      // Get data for current point
      const double current_angle = std::get<0>(angles[i]);
      const size_t current_index = std::get<2>(angles[i]);

      // Add point or replace last point
      if (std::abs(current_angle - last_angle) > Constants::epsilon)
        filtered_indices.push_back(current_index);
      else
        filtered_indices[filtered_indices.size() - 1] = current_index;

      // Update last index
      last_angle = current_angle;
    }

    // Create stack of points and add first three candidates
    std::stack<size_t> convex_hull;
    convex_hull.push(base_index);
    convex_hull.push(filtered_indices[0]);
    convex_hull.push(filtered_indices[1]);

    // Graham-Scan: Push candidates to stack and pop until
    // we have a left turn
    for (size_t i = 2; i < filtered_indices.size(); i++)
    {
      // Get next point
      const size_t i2 = filtered_indices[i];
      const Vector2D &p2 = points[i2];

      // Keep popping from stack until we see a left turn
      while (true)
      {
        // Get last two points from stack
        const size_t i1 = convex_hull.top();
        convex_hull.pop();
        const size_t i0 = convex_hull.top();
        const Vector2D &p0 = points[i0];
        const Vector2D &p1 = points[i1];

        // Check orientation, keep p1 if orientation is positive
        if (orient_2d(p0, p1, p2) > Constants::epsilon)
        {
          convex_hull.push(i1);
          break;
        }
      }

      // Push next candidate to stack
      convex_hull.push(i2);
    }

    // Extract polygon points from stack
    Polygon polygon;
    while (!convex_hull.empty())
    {
      polygon.vertices.push_back(points[convex_hull.top()]);
      convex_hull.pop();
    }

    // Reverse polygon to make it counter-clockwise
    std::reverse(polygon.vertices.begin(), polygon.vertices.end());

    return polygon;
  }

  // Compute tetrahedron volume (3D)
  static double tetrahedron_volume(const Vector3D &v0, const Vector3D &v1, const Vector3D &v2,
                                   const Vector3D &v3)
  {
    return std::abs((v1 - v0).dot((v2 - v0).cross(v3 - v0))) / 6.0;
  }

  // Compute tetrahedron circumcenter (3D)
  static Vector3D tetrahedron_circumcenter(const Vector3D &v0, const Vector3D &v1,
                                           const Vector3D &v2, const Vector3D &v3)
  {
    // Translate so v0 is at origin
    const Vector3D &B = v1 - v0;
    const Vector3D &C = v2 - v0;
    const Vector3D &D = v3 - v0;

    // Compute squared magnitudes
    const double B2 = B.squared_magnitude();
    const double C2 = C.squared_magnitude();
    const double D2 = D.squared_magnitude();

    // Create linear system
    double M[3][3] = {{B.x, B.y, B.z}, {C.x, C.y, C.z}, {D.x, D.y, D.z}};
    double R[3] = {0.5 * B2, 0.5 * C2, 0.5 * D2};

    // Compute determinant of M
    const double det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) -
                       M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
                       M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);

    // Cramer's rule for each variable

    auto cramer = [&](int col)
    {
      double _M[3][3];
      for (int i = 0; i < 3; i++)
      {
        for (int j = 0; j < 3; j++)
        {
          _M[i][j] = M[i][j];
        }
      }
      for (int i = 0; i < 3; i++)
      {
        _M[i][col] = R[i];
      }
      double d = _M[0][0] * (_M[1][1] * _M[2][2] - _M[1][2] * _M[2][1]) -
                 _M[0][1] * (_M[1][0] * _M[2][2] - _M[1][2] * _M[2][0]) +
                 _M[0][2] * (_M[1][0] * _M[2][1] - _M[1][1] * _M[2][0]);
      return d;
    };

    // Solve for each variable
    const double detX = cramer(0);
    const double detY = cramer(1);
    const double detZ = cramer(2);
    const double x = detX / det;
    const double y = detY / det;
    const double z = detZ / det;

    // Return circumcenter
    return Vector3D(v0.x + x, v0.y + y, v0.z + z);
  }

  // Compute face area (3D)
  static double face_area(const Vector3D &v0, const Vector3D &v1, const Vector3D &v2)
  {
    return 0.5 * (v1 - v0).cross(v2 - v0).magnitude();
  }

  // Compute normalized aspect ratio of tetrahedron (3D)
  static double aspect_ratio(const Vector3D &v0, const Vector3D &v1, const Vector3D &v2,
                             const Vector3D &v3)
  {
    // Compute volume
    const double V = tetrahedron_volume(v0, v1, v2, v3);

    // Handle degenerate tetrahedra (zero volume)
    if (V == 0)
      return 0.0;

    // Compute radius of inscribed
    const double A0 = face_area(v1, v2, v3);
    const double A1 = face_area(v0, v2, v3);
    const double A2 = face_area(v0, v1, v3);
    const double A3 = face_area(v0, v1, v2);
    const double A = A0 + A1 + A2 + A3;
    const double r = 3.0 * V / A;

    // Compute radius of circumscribed sphere
    const Vector3D c = tetrahedron_circumcenter(v0, v1, v2, v3);
    const double R = (c - v0).magnitude();

    // Compute aspect ratio (normalized by regular tetrahedron)
    const double ar = (R / r) / 3.0;

    return ar;
  }

  // Compute normalized aspect ratios of volume mesh
  static std::vector<double> aspect_ratios(const VolumeMesh &mesh)
  {
    // Compute aspect ratios of all tetrahedra
    std::vector<double> aspect_ratios;
    aspect_ratios.reserve(mesh.cells.size());
    for (const auto &cell : mesh.cells)
    {
      const Vector3D &v0 = mesh.vertices[cell.v0];
      const Vector3D &v1 = mesh.vertices[cell.v1];
      const Vector3D &v2 = mesh.vertices[cell.v2];
      const Vector3D &v3 = mesh.vertices[cell.v3];
      const double ar = aspect_ratio(v0, v1, v2, v3);
      aspect_ratios.push_back(ar);
    }

    return aspect_ratios;
  }

  // Compute normalized aspect ratio of volume mesh: min, max, median (3D)
  static std::tuple<double, double, double> aspect_ratio(const VolumeMesh &mesh)
  {
    // Compute aspect ratios of all tetrahedra
    std::vector<double> _aspect_ratios = aspect_ratios(mesh);

    // Compute min, max and median aspect ratios
    double min_ar = *std::min_element(_aspect_ratios.begin(), _aspect_ratios.end());
    double max_ar = *std::max_element(_aspect_ratios.begin(), _aspect_ratios.end());
    std::sort(_aspect_ratios.begin(), _aspect_ratios.end());
    double median_ar = _aspect_ratios[_aspect_ratios.size() / 2];

    return std::make_tuple(min_ar, max_ar, median_ar);
  }
};

} // namespace DTCC_BUILDER

#endif
