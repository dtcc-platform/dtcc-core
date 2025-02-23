// Copyright (C) 2024 Anders Logg, Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_TRIANGULATE_BUILDER_H
#define DTCC_TRIANGULATE_BUILDER_H

#include <cmath>
#include <iostream>
#include <map>
#include <stack>
#include <tuple>
#include <vector>

extern "C"
{
#include <triangle.h>
}
#include <earcut.hpp>

#include "Eigen/Eigen"
#include "Eigen/Geometry"

#include "Logging.h"
#include "Timer.h"

#include "model/Mesh.h"
#include "model/Surface.h"
#include "model/Vector.h"



namespace DTCC_BUILDER
{

class Triangulate
{
public:

  static void call_earcut(Mesh &mesh, const Surface &surface, bool sort_triangles = true) {
    auto area = Geometry::surface_area(surface);
    //    info("surface area " + str(area) + " m^2");
    if (isnan(area) || area < 1e-3)
      return;

    auto [projected_polygon, transform_inv] = project_surface(surface);
    std::vector<std::vector<std::array<double,2>>> earcut_polygon;
    std::vector<std::array<double,2>> shell;
    for (auto const &v : projected_polygon.vertices)
    {
      shell.push_back({v.x, v.y});
    }
    earcut_polygon.push_back(shell);
    for (auto const &hole : projected_polygon.holes)
    {
      std::vector<std::array<double,2>> hole_;
      hole_.reserve(hole.size());
      for (auto const &v : hole)
      {
        hole_.push_back({v.x, v.y});
      }
      earcut_polygon.push_back(hole_);
    }
    auto tri_indices = mapbox::earcut<std::size_t>(earcut_polygon);
    mesh.faces.reserve(tri_indices.size() / 3);
    for (size_t i = 0; i < tri_indices.size(); i += 3)
    {
      mesh.faces.emplace_back(tri_indices[i], tri_indices[i + 1], tri_indices[i + 2], false);
    }
    for (auto const &linestrings: earcut_polygon)
    {
      for (auto const &v : linestrings)
      {
        auto e_v = Eigen::Vector3d(v[0], v[1], 0);
        auto e_v_prime = transform_inv * e_v;
        mesh.vertices.emplace_back(e_v_prime.x(), e_v_prime.y(), e_v_prime.z());
      }
    }
    mesh.calc_normals();
    const auto normal = Geometry::surface_normal(surface);
    for (size_t i = 0; i < mesh.normals.size(); i++)
    {
      if (!mesh.normals[i].close_to(normal))
      {
        std::swap(mesh.faces[i].v1, mesh.faces[i].v2);
        mesh.normals[i] = -mesh.normals[i];
      }
    }


  }
  static void call_triangle(Mesh &mesh, const Surface &surface, double max_mesh_size,
                            double min_mesh_angle, bool sort_triangles = true)
  {
    std::vector<std::vector<Vector2D>> sd;

    auto area = Geometry::surface_area(surface);
    //    info("surface area " + str(area) + " m^2");
    if (isnan(area) || area < 1e-3)
      return;

    auto [projected_polygon, transform_inv] = project_surface(surface);



    std::vector<double> sd_size;
    call_triangle(mesh, projected_polygon.vertices, projected_polygon.holes, sd_size, max_mesh_size, min_mesh_angle,
                  sort_triangles);
    for (auto &v : mesh.vertices)
    {
      auto e_v = Eigen::Vector3d(v.x, v.y, 0);
      auto e_v_prime = transform_inv * e_v;
      v.x = e_v_prime.x();
      v.y = e_v_prime.y();
      v.z = e_v_prime.z();
    }

    mesh.calc_normals();
    auto normal = Geometry::surface_normal(surface);
    for (size_t i = 0; i < mesh.normals.size(); i++)
    {
      if (!mesh.normals[i].close_to(normal))
      {
        std::swap(mesh.faces[i].v1, mesh.faces[i].v2);
        mesh.normals[i] = -mesh.normals[i];
      }
    }
  }

  static void fast_mesh(Mesh &mesh, const Surface &surface)
  {
    bool has_holes = surface.holes.size() > 0;
    if (surface.vertices.size() < 3)
      return;
    if (surface.vertices.size() == 3 && !has_holes)
    {
      mesh.vertices = surface.vertices;
      mesh.faces.push_back(Simplex2D(0, 1, 2, false));
      mesh.sort_vertices();
      return;
    }

    if (Geometry::is_convex(surface) && !has_holes)
    {
      mesh.vertices = surface.vertices;
      for (size_t i = 1; i < surface.vertices.size() - 1; i++)
      {
        mesh.faces.push_back(Simplex2D(0, i, i + 1, false));
      }
      mesh.sort_vertices();
    }
    else
    {
      call_earcut(mesh, surface);
    }
  }

  // Call Triangle to compute 2D mesh
  static void call_triangle(Mesh &mesh, const std::vector<Vector2D> &boundary,
                            const std::vector<std::vector<Vector2D>> &sub_domains,
                            const std::vector<double> &subdomain_triangle_size,
                            double max_mesh_size, double min_mesh_angle, bool sort_triangles)
  {
    Timer timer("call_triangle");

    // Set area constraint to control mesh size
    const double max_area = 0.5 * max_mesh_size * max_mesh_size;

    // Set input switches for Triangle
    std::string triswitches = "zQp";
    if (min_mesh_angle > 0)
      triswitches += "q" + str(min_mesh_angle, 3);
    if (max_mesh_size > 0 && subdomain_triangle_size.size() == 0)
      triswitches += "a" + str(max_area, 3);
    if (subdomain_triangle_size.size() > 0)
      triswitches += "a";

    debug("Triangle switches: " + triswitches);

    // Convert to C-style string
    char *triswitches_c = new char[triswitches.length() + 1];
    std::strcpy(triswitches_c, triswitches.c_str());

    // z = use zero-based numbering
    // p = use polygon input (segments)
    // q = control mesh quality
    //
    // Note that the minimum angle (here 25) should be
    // as large as possible for high quality meshes but
    // it should be less than 28.6 degrees to guarantee
    // that Triangle terminates. Default is 20 degrees.

    // Create input data structure for Triangle
    struct triangulateio in = create_triangle_io();

    // Check for duplicate points
    // size_t duplicate_vertices = 0;
    for (auto polygon : sub_domains)
    {
      auto first = polygon.front();
      auto last = polygon.back();
      if (first.close_to(last))
      {
        polygon.pop_back();
        // duplicate_vertices++;
      }
    }
//    if (duplicate_vertices > 0)
//      debug("Removed " + str(duplicate_vertices) + " duplicate vertices");

    // Set number of points
    size_t num_points = boundary.size();
    //    info("triangluate with " + str(num_points) + " points");
    for (auto const &innerPolygon : sub_domains)
    {
      num_points += innerPolygon.size();
    }
    in.numberofpoints = num_points;

    // Set points
    in.pointlist = new double[2 * num_points];
    {
      size_t k = 0;
      for (auto const &p : boundary)
      {
        in.pointlist[k++] = p.x;
        in.pointlist[k++] = p.y;
      }
      for (auto const &innerPolygon : sub_domains)
      {
        for (auto const &p : innerPolygon)
        {
          in.pointlist[k++] = p.x;
          in.pointlist[k++] = p.y;
        }
      }
    }

    // Set number of segments
    const size_t num_segments = num_points;
    in.numberofsegments = num_segments;

    // Set segments
    in.segmentlist = new int[2 * num_segments];
    {
      size_t k = 0;
      size_t n = 0;
      for (size_t j = 0; j < boundary.size(); j++)
      {
        const size_t j0 = j;
        const size_t j1 = (j + 1) % boundary.size();
        in.segmentlist[k++] = n + j0;
        in.segmentlist[k++] = n + j1;
      }
      n += boundary.size();
      for (size_t i = 0; i < sub_domains.size(); i++)
      {
        for (size_t j = 0; j < sub_domains[i].size(); j++)
        {
          const size_t j0 = j;
          const size_t j1 = (j + 1) % sub_domains[i].size();
          in.segmentlist[k++] = n + j0;
          in.segmentlist[k++] = n + j1;
        }
        n += sub_domains[i].size();
      }
    }

    if (subdomain_triangle_size.size() > 0)
    {
      in.regionlist = new double[4 * (1 + sub_domains.size())];
      auto boundary_center = Geometry::polygon_center_2d(Polygon(boundary));
      auto corner_to_center = boundary_center - boundary[0];
      corner_to_center.normalize();
      auto inside_boundary = boundary[0] + (corner_to_center * 0.2);

      in.regionlist[0] = inside_boundary.x;
      in.regionlist[1] = inside_boundary.y;
      in.regionlist[2] = 1;
      in.regionlist[3] = max_area;

      auto k = 4;

      for (size_t i = 0; i < sub_domains.size(); i++)
      {
        auto inner_polygon = sub_domains[i];
        auto c = Geometry::point_inside_polygon_2d(Polygon(inner_polygon));
        double max_inner_area = subdomain_triangle_size[i] * subdomain_triangle_size[i] * 0.5;
        in.regionlist[k++] = c.x;
        in.regionlist[k++] = c.y;
        in.regionlist[k++] = 1;
        in.regionlist[k++] = max_inner_area;
      }
      in.numberofregions = 1 + sub_domains.size();
    }

    // Note: This is how set holes but it's not used here since we
    // need the triangles for the interior *above* the buildings.

    /*
    // Set number of holes
    const size_t numHoles = SubDomains.size();
    in.numberofholes = numHoles;

    // Set holes. Note that we assume that we can get an
    // interior point of each hole (inner polygon) by computing
    // its center of mass.
    in.holelist = new double[2 * numHoles];
    {
    size_t k = 0;
    Vector2D c;
    for (auto const & InnerPolygon : SubDomains)
    {
    for (auto const & p : InnerPolygon)
    {
    c += p;
    }
    c /= InnerPolygon.size();
    in.holelist[k++] = c.x;
    in.holelist[k++] = c.y;
    }
    }
    */

    // Prepare output data for Triangl;e
    struct triangulateio out = create_triangle_io();
    struct triangulateio vorout = create_triangle_io();
    //    print_triangle_io(in);
    // Call Triangle
    triangulate(triswitches_c, &in, &out, &vorout);
    delete[] triswitches_c;

    // Uncomment for debugging
    // print_triangle_io(out);
    // print_triangle_io(vorout);

    // Extract points
    mesh.vertices.reserve(out.numberofpoints);
    for (int i = 0; i < out.numberofpoints; i++)
    {
      Vector3D p(out.pointlist[2 * i], out.pointlist[2 * i + 1], 0.0);
      mesh.vertices.push_back(p);
    }

    // Extract triangles
    mesh.faces.reserve(out.numberoftriangles);
    for (int i = 0; i < out.numberoftriangles; i++)
    {
      // Note the importance of creating a sorted simplex here!
      Simplex2D t(out.trianglelist[3 * i], out.trianglelist[3 * i + 1], out.trianglelist[3 * i + 2],
                  sort_triangles);
      mesh.faces.push_back(t);
    }

    // Free memory
    // trifree(&out); // causes segfault
    delete[] in.pointlist;
    delete[] in.segmentlist;
    delete[] in.holelist;
  }

private:

  static std::tuple<Polygon, Eigen::Transform<double,3,1>> project_surface(const Surface &surface)
  {
    const auto z_normal = Eigen::Vector3d(0, 0, 1);
    auto normal = Geometry::surface_normal(surface);
    auto e_norm = Eigen::Vector3d(normal.x, normal.y, normal.z);
    auto centroid = Geometry::surface_centroid(surface);
    auto e_centroid = Eigen::Vector3d(centroid.x, centroid.y, centroid.z);

    auto transform = Eigen::Transform<double, 3, Eigen::Isometry>();
    auto translation = Eigen::Translation3d(-e_centroid);
    auto rotation = Eigen::Quaterniond::FromTwoVectors(e_norm, z_normal);

    transform = rotation * translation;
    auto transform_inv = transform.inverse();
    // std::cout << "trans_matrix " << trans_matrix << std::endl;

    // std::cout << "transform " << transform << std::endl;
    // auto transform_inv = transform.inverse();
    // std::cout << "transform inv" << transform_inv << std::endl;
    Polygon projected_polygon;
    for (const auto &v : surface.vertices)
    {
      auto e_v = Eigen::Vector3d(v.x, v.y, v.z);
      auto e_v_prime = transform * e_v;
      auto projected_v = Vector2D(e_v_prime.x(), e_v_prime.y());
      if (projected_polygon.vertices.empty() || (!projected_v.close_to(projected_polygon.vertices.front()) &&
                                            !projected_v.close_to(projected_polygon.vertices.back())))
        projected_polygon.vertices.push_back(Vector2D(e_v_prime.x(), e_v_prime.y()));
    }
//    size_t removed_vertices = surface.vertices.size() - projected_polygon.vertices.size();
//    if (removed_vertices > 0)
//      info("Removed " + str(removed_vertices) + " duplicate vertices");

    for (const auto &hole : surface.holes)
    {
      std::vector<Vector2D> projected_hole;
      for (const auto &v : hole)
      {
        auto e_v = Eigen::Vector3d(v.x, v.y, v.z);
        auto e_v_prime = transform * e_v;
        auto projected_v = Vector2D(e_v_prime.x(), e_v_prime.y());
        if (projected_hole.size() == 0 || (!projected_v.close_to(projected_hole.front()) &&
                                           !projected_v.close_to(projected_hole.back())))
          projected_hole.push_back(Vector2D(e_v_prime.x(), e_v_prime.y()));
      }
      projected_polygon.holes.push_back(projected_hole);
    }
    return std::make_tuple(projected_polygon, transform_inv);
  }

  // Create and reset Triangle I/O data structure
  static struct triangulateio create_triangle_io()
  {
    struct triangulateio io;

    io.pointlist = nullptr;
    io.pointmarkerlist = nullptr;
    io.pointmarkerlist = nullptr;
    io.numberofpoints = 0;
    io.numberofpointattributes = 0;
    io.trianglelist = nullptr;
    io.triangleattributelist = nullptr;
    io.trianglearealist = nullptr;
    io.neighborlist = nullptr;
    io.numberoftriangles = 0;
    io.numberofcorners = 0;
    io.numberoftriangleattributes = 0;
    io.segmentlist = nullptr;
    io.segmentmarkerlist = nullptr;
    io.numberofsegments = 0;
    io.holelist = nullptr;
    io.numberofholes = 0;
    io.regionlist = nullptr;
    io.numberofregions = 0;
    io.edgelist = nullptr;
    io.edgemarkerlist = nullptr;
    io.normlist = nullptr;
    io.numberofedges = 0;

    return io;
  }

  // print triangle I/O data
  static void print_triangle_io(const struct triangulateio &io)
  {
    info("Triangle I/O data: ");
    info("  pointlist = " + str(reinterpret_cast<std::uintptr_t>(io.pointlist)));
    info("  pointmarkerlist = " + str(reinterpret_cast<std::uintptr_t>(io.pointmarkerlist)));
    if (io.pointmarkerlist)
    {
      std::stringstream string_builder{};
      string_builder << "   ";
      for (int i = 0; i < io.numberofpoints; i++)
        string_builder << " " << io.pointmarkerlist[i];
      string_builder << std::endl;
      info(string_builder.str());
    }
    info("  numberofpoints = " + str(io.numberofpoints));
    info("  numberofpointattributes = " + str(io.numberofpointattributes));
    info("  trianglelist = " + str(reinterpret_cast<std::uintptr_t>(io.trianglelist)));
    info("  triangleattributelist = " +
         str(reinterpret_cast<std::uintptr_t>(io.triangleattributelist)));
    info("  trianglearealist = " + str(reinterpret_cast<std::uintptr_t>(io.trianglearealist)));
    info("  neighborlist = " + str(reinterpret_cast<std::uintptr_t>(io.neighborlist)));
    info("  numberoftriangles = " + str(io.numberoftriangles));
    info("  numberofcorners = " + str(io.numberofcorners));
    info("  numberoftriangleattributes = " + str(io.numberoftriangleattributes));
    info("  segmentlist = " + str(reinterpret_cast<std::uintptr_t>(io.segmentlist)));
    info("  segmentmarkerlist = " + str(reinterpret_cast<std::uintptr_t>(io.segmentmarkerlist)));
    if (io.segmentmarkerlist)
    {
      std::stringstream string_builder{};
      string_builder << "   ";
      for (int i = 0; i < io.numberofsegments; i++)
        string_builder << " " << io.segmentmarkerlist[i];
      string_builder << std::endl;
      info(string_builder.str());
    }
    info("  numberofsegments = " + str(io.numberofsegments));
    info("  holelist = " + str(reinterpret_cast<std::uintptr_t>(io.holelist)));
    info("  numberofholes = " + str(io.numberofholes));
    info("  regionlist = " + str(reinterpret_cast<std::uintptr_t>(io.regionlist)));
    info("  numberofregions = " + str(io.numberofregions));
    info("  edgelist = " + str(reinterpret_cast<std::uintptr_t>(io.edgelist)));
    info("  edgemarkerlist = " + str(reinterpret_cast<std::uintptr_t>(io.edgemarkerlist)));
    info("  normlist = " + str(reinterpret_cast<std::uintptr_t>(io.normlist)));
    info("  numberofedges = " + str(io.numberofedges));
  }
};

} // namespace DTCC_BUILDER

#endif