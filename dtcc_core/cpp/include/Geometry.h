// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_GEOMETRY_H
#define DTCC_GEOMETRY_H

#include <algorithm>
#include <cmath>
#include <iso646.h>
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

namespace DTCC_BUILDER
{

class Geometry
{
public:
  // Compute squared norm (3D)
  static double squared_norm_3d(const Vector3D &v) { return dot_3d(v, v); }

  // Compute norm (3D)
  static double norm_3d(const Vector3D &v) { return std::sqrt(squared_norm_3d(v)); }

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

  // Compute squared distance between points (2D)
  static double squared_distance_2d(const Vector2D &p, const Vector2D &q)
  {
    const double dx = p.x - q.x;
    const double dy = p.y - q.y;
    return dx * dx + dy * dy;
  }

  // Compute orientation of point q relative to edge (p0, p1) (2D)
  static double orient_2d(const Vector2D &p0, const Vector2D &p1, const Vector2D &q)
  {
    const Vector2D u(p0, p1);
    const Vector2D v(p0, q);
    return u.x * v.y - u.y * v.x;
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


};

} // namespace DTCC_BUILDER

#endif
