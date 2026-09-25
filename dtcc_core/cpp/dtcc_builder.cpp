// Copyright (C) 2023 Dag Wästberg
// Licensed under the MIT License
//
// Modified by Anders Logg 2023

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <type_traits>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "BuildingProcessor.h"
#include "Intersection.h"
#include "MeshBuilder.h"
#include "MeshProcessor.h"
#include "VertexSmoother.h"
#include "model/GridField.h"
#include "model/Mesh.h"
#include "model/Polygon.h"
#include "model/Simplices.h"
#include "model/Vector.h"
#include "model/VolumeMesh.h"
#include "terrain-mesher/zemlya.hpp"

namespace py = pybind11;

namespace DTCC_BUILDER
{

Polygon create_polygon(py::list vertices, py::list holes)
{
  Polygon poly;
  for (size_t i = 0; i < vertices.size(); i++)
  {
    auto pt = vertices[i].cast<py::tuple>();
    poly.vertices.push_back(Vector2D(pt[0].cast<double>(), pt[1].cast<double>()));
  }
  if (holes.size() > 0)
  {
    for (size_t i = 0; i < holes.size(); i++)
    {
      auto hl = holes[i].cast<py::list>();
      std::vector<Vector2D> hole;
      for (size_t j = 0; j < hl.size(); j++)
      {
        auto pt = hl[j].cast<py::tuple>();
        hole.push_back(Vector2D(pt[0].cast<double>(), pt[1].cast<double>()));
      }
      poly.holes.push_back(hole);
    }
  }
  return poly;
}

namespace
{

py::ssize_t mesh_array_rows(const py::array &array, py::ssize_t width, const char *name)
{
  if (array.ndim() == 1 && array.shape(0) == 0)
    return 0;
  if (array.ndim() != 2 || array.shape(1) != width)
    throw py::value_error(std::string(name) + " must have shape (N, " +
                          std::to_string(width) + ") or (0,)");
  return array.shape(0);
}

// NumPy can expose unaligned buffers. memcpy keeps scalar reads safe even when
// a contiguous typed conversion can reuse the original buffer.
template <typename T>
T mesh_array_value(const T *data, py::ssize_t index)
{
  T value;
  std::memcpy(&value, data + index, sizeof(T));
  return value;
}

template <typename Copy>
void copy_integer_mesh_array(const py::array &array, const char *name, Copy copy)
{
  const char kind = array.dtype().kind();
  // Model defaults use empty floating arrays for cells and markers.
  if (array.size() == 0 && kind == 'f')
    return;
  if ((kind != 'i' && kind != 'u') || array.itemsize() > sizeof(std::uint64_t))
    throw py::value_error(std::string(name) + " must contain integer values");
  if (array.size() == 0)
    return;
  // Widen before checking so neither unsigned overflow nor negative wrapping
  // can hide invalid input. Typed conversion handles byte order and strides.
  if (kind == 'i')
    copy(py::array_t<std::int64_t, py::array::c_style | py::array::forcecast>(array));
  else
    copy(py::array_t<std::uint64_t, py::array::c_style | py::array::forcecast>(array));
}

std::vector<Vector3D> copy_mesh_vectors(const py::array &array, const char *name,
                                       py::ssize_t expected_count = -1)
{
  const auto count = mesh_array_rows(array, 3, name);
  if (count != 0 && expected_count >= 0 && count != expected_count)
    throw py::value_error(std::string(name) + " must be empty or contain one vector per face");
  const char kind = array.dtype().kind();
  if (kind != 'f' && kind != 'i' && kind != 'u')
    throw py::value_error(std::string(name) + " must contain finite real values");
  py::array_t<double, py::array::c_style | py::array::forcecast> values(array);
  std::vector<Vector3D> result;
  result.reserve(count);
  for (py::ssize_t i = 0; i < count; ++i)
  {
    const double x = mesh_array_value(values.data(), 3 * i);
    const double y = mesh_array_value(values.data(), 3 * i + 1);
    const double z = mesh_array_value(values.data(), 3 * i + 2);
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
      throw py::value_error(std::string(name) + " must contain finite real values");
    result.emplace_back(x, y, z);
  }
  return result;
}

template <std::size_t Width, typename Simplex>
std::vector<Simplex> copy_mesh_connectivity(const py::array &array,
                                           std::size_t num_vertices, const char *name)
{
  const auto count = mesh_array_rows(array, Width, name);
  std::vector<Simplex> result;
  result.reserve(count);
  copy_integer_mesh_array(array, name, [&](const auto &values)
  {
    for (py::ssize_t i = 0; i < count; ++i)
    {
      std::size_t indices[Width];
      for (std::size_t j = 0; j < Width; ++j)
      {
        const auto value = mesh_array_value(values.data(), Width * i + j);
        if constexpr (std::is_signed_v<decltype(value)>)
        {
          if (value < 0)
            throw py::value_error(std::string(name) + " contains an out-of-range vertex index");
        }
        if (static_cast<std::uint64_t>(value) >= num_vertices)
          throw py::value_error(std::string(name) + " contains an out-of-range vertex index");
        indices[j] = static_cast<std::size_t>(value);
      }
      if constexpr (Width == 3)
        result.emplace_back(indices[0], indices[1], indices[2]);
      else
        result.emplace_back(indices[0], indices[1], indices[2], indices[3]);
    }
  });
  return result;
}

std::vector<int> copy_mesh_markers(const py::array &markers, std::size_t num_elements)
{
  if (markers.ndim() != 1 || (markers.size() != 0 &&
                              static_cast<std::size_t>(markers.size()) != num_elements))
    throw py::value_error("markers must be empty or a vector with one integer per face/cell");
  std::vector<int> result;
  result.reserve(markers.size());
  copy_integer_mesh_array(markers, "markers", [&](const auto &values)
  {
    for (py::ssize_t i = 0; i < values.size(); ++i)
    {
      const auto value = mesh_array_value(values.data(), i);
      if constexpr (std::is_signed_v<decltype(value)>)
      {
        if (value < std::numeric_limits<int>::min())
          throw py::value_error("markers must fit in a native int");
      }
      if (value > std::numeric_limits<int>::max())
        throw py::value_error("markers must fit in a native int");
      result.push_back(static_cast<int>(value));
    }
  });
  return result;
}

} // namespace

Mesh create_mesh(py::array vertices, py::array faces, py::array markers, py::array normals)
{
  Mesh mesh;
  mesh.vertices = copy_mesh_vectors(vertices, "vertices");
  mesh.faces = copy_mesh_connectivity<3, Simplex2D>(faces, mesh.vertices.size(), "faces");
  mesh.markers = copy_mesh_markers(markers, mesh.faces.size());
  mesh.normals = copy_mesh_vectors(normals, "normals", mesh.faces.size());
  return mesh;
}

VolumeMesh create_volume_mesh(py::array vertices, py::array cells, py::array markers)
{
  VolumeMesh mesh;
  mesh.vertices = copy_mesh_vectors(vertices, "vertices");
  mesh.cells = copy_mesh_connectivity<4, Simplex3D>(cells, mesh.vertices.size(), "cells");
  mesh.markers = copy_mesh_markers(markers, mesh.cells.size());
  return mesh;
}

py::tuple mesh_as_arrays(const Mesh &mesh)
{
  py::array_t<double> vertices(mesh.vertices.size() * 3);
  py::array_t<size_t> faces(mesh.faces.size() * 3);
  py::array_t<int> markers(mesh.markers.size());
  py::array_t<double> normals(mesh.normals.size() * 3);
  auto *n = normals.mutable_data();
  auto *v = vertices.mutable_data();
  auto *f = faces.mutable_data();
  auto *m = markers.mutable_data();
  for (size_t i = 0; i < mesh.vertices.size(); ++i)
  {
    v[3 * i] = mesh.vertices[i].x;
    v[3 * i + 1] = mesh.vertices[i].y;
    v[3 * i + 2] = mesh.vertices[i].z;
  }
  for (size_t i = 0; i < mesh.faces.size(); ++i)
  {
    f[3 * i] = mesh.faces[i].v0;
    f[3 * i + 1] = mesh.faces[i].v1;
    f[3 * i + 2] = mesh.faces[i].v2;
  }
  for (size_t i = 0; i < mesh.markers.size(); ++i)
    m[i] = mesh.markers[i];
  for (size_t i = 0; i < mesh.normals.size(); ++i)
  {
    n[3 * i] = mesh.normals[i].x;
    n[3 * i + 1] = mesh.normals[i].y;
    n[3 * i + 2] = mesh.normals[i].z;
  }
  return py::make_tuple(vertices, faces, markers, normals);
}

py::tuple volume_mesh_as_arrays(const VolumeMesh &mesh)
{
  py::array_t<double> vertices(mesh.vertices.size() * 3);
  py::array_t<size_t> cells(mesh.cells.size() * 4);
  py::array_t<int> markers(mesh.markers.size());
  auto *v = vertices.mutable_data();
  auto *c = cells.mutable_data();
  auto *m = markers.mutable_data();
  for (size_t i = 0; i < mesh.vertices.size(); ++i)
  {
    v[3 * i] = mesh.vertices[i].x;
    v[3 * i + 1] = mesh.vertices[i].y;
    v[3 * i + 2] = mesh.vertices[i].z;
  }
  for (size_t i = 0; i < mesh.cells.size(); ++i)
  {
    c[4 * i] = mesh.cells[i].v0;
    c[4 * i + 1] = mesh.cells[i].v1;
    c[4 * i + 2] = mesh.cells[i].v2;
    c[4 * i + 3] = mesh.cells[i].v3;
  }
  for (size_t i = 0; i < mesh.markers.size(); ++i)
    m[i] = mesh.markers[i];
  return py::make_tuple(vertices, cells, markers);
}

Surface create_surface(py::array_t<double> vertices, py::list holes)
{
  Surface surface;
  auto verts_r = vertices.unchecked<2>();
  size_t num_vertices = verts_r.shape(0);
  for (size_t i = 0; i < num_vertices; i++)
  {
    surface.vertices.push_back(Vector3D(verts_r(i, 0), verts_r(i, 1), verts_r(i, 2)));
  }
  for (size_t i = 0; i < holes.size(); i++)
  {
    auto hole = holes[i].cast<py::array_t<double>>();
    auto hole_r = hole.unchecked<2>();
    std::vector<Vector3D> hole_vertices;
    size_t num_hole_vertices = hole_r.shape(0);
    for (size_t j = 0; j < num_hole_vertices; j++)
    {
      hole_vertices.push_back(Vector3D(hole_r(j, 0), hole_r(j, 1), hole_r(j, 2)));
    }
    surface.holes.push_back(hole_vertices);
  }
  return surface;
}

py::list points_in_polygons(const py::array &pts, const std::vector<Polygon> &polygons)
{
  auto pc = copy_mesh_vectors(pts, "points");
  auto pips = PointCloudProcessor::points_in_polygons(pc, polygons);
  py::list in_polygons_list;
  for (auto const &pip : pips)
  {
    py::array_t<size_t> indices(pip.size());
    auto *data = indices.mutable_data();
    for (size_t i = 0; i < pip.size(); i++)
    {
      data[i] = pip[i];
    }
    in_polygons_list.append(indices);
  }

  return in_polygons_list;
}

py::list extract_building_points(std::vector<Polygon> &buildings, const py::array_t<double> &pts,
                                 bool statistical_outlier_remover, size_t neighbors,
                                 double outlier_margin)
{
  py::list roof_points;
  auto pts_r = pts.unchecked<2>();
  size_t pt_count = pts_r.shape(0);
  std::vector<Vector3D> pc;
  for (size_t i = 0; i < pt_count; i++)
  {
    pc.push_back(Vector3D(pts_r(i, 0), pts_r(i, 1), pts_r(i, 2)));
  }

  auto _roof_points = BuildingProcessor::extract_building_points(buildings, pc);
  if (statistical_outlier_remover)
  {
    for (auto &rp : _roof_points)
    {
      PointCloudProcessor::statistical_outlier_remover(rp, neighbors, outlier_margin);
    }
  }
  for (auto const &rp : _roof_points)
  {
    py::array_t<double> pts(rp.size() * 3);
    for (size_t i = 0; i < rp.size(); i++)
    {
      pts.mutable_at(i * 3) = rp[i].x;
      pts.mutable_at(i * 3 + 1) = rp[i].y;
      pts.mutable_at(i * 3 + 2) = rp[i].z;
    }
    pts = pts.reshape(std::vector<long>{static_cast<long>(rp.size()), 3});
    roof_points.append(pts);
  }
  return roof_points;
}

MultiSurface create_multisurface(py::list surfaces)
{
  MultiSurface multi_surface;
  for (size_t i = 0; i < surfaces.size(); i++)
  {
    auto surface = surfaces[i].cast<Surface>();
    multi_surface.surfaces.push_back(surface);
  }
  return multi_surface;
}

GridField create_gridfield(py::array data, py::tuple bounds, py::ssize_t xsize, py::ssize_t ysize)
{
  if (xsize <= 0 || ysize <= 0)
    throw py::value_error("Grid dimensions must be positive");
  if (data.ndim() != 1 || data.size() / xsize != ysize || data.size() % xsize != 0)
    throw py::value_error("Grid data must be a flat array matching the grid dimensions");
  const char kind = data.dtype().kind();
  if (kind != 'f' && kind != 'i' && kind != 'u')
    throw py::value_error("Grid data must contain finite real values; fill missing data before meshing");
  if (bounds.size() != 4)
    throw py::value_error("Grid bounds must contain xmin, ymin, xmax, ymax");
  GridField grid_field;
  double px = bounds[0].cast<double>();
  double py = bounds[1].cast<double>();
  double qx = bounds[2].cast<double>();
  double qy = bounds[3].cast<double>();
  if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(qx) || !std::isfinite(qy) ||
      qx <= px || qy <= py)
    throw py::value_error("Grid bounds must be finite with positive width and height");
  auto bbox = BoundingBox2D(Vector2D(px, py), Vector2D(qx, qy));

  grid_field.grid.bounding_box = bbox;
  grid_field.grid.xstep = (qx - px) / xsize;
  grid_field.grid.ystep = (qy - py) / ysize;
  if (!std::isfinite(grid_field.grid.xstep) || !std::isfinite(grid_field.grid.ystep) ||
      grid_field.grid.xstep <= 0 || grid_field.grid.ystep <= 0)
    throw py::value_error("Grid pixel spacing must be finite and positive");
  grid_field.grid.cell_centered = true;

  grid_field.grid.xsize = xsize;
  grid_field.grid.ysize = ysize;

  py::array_t<double, py::array::c_style | py::array::forcecast> values(data);
  grid_field.values.reserve(values.size());
  for (py::ssize_t i = 0; i < values.size(); i++)
  {
    const double value = mesh_array_value(values.data(), i);
    if (!std::isfinite(value))
      throw py::value_error("Grid data must contain finite real values; fill missing data before meshing");
    grid_field.values.push_back(value);
  }

  return grid_field;
}

terrain_mesher::core::RasterDouble gridfield_to_raster(const GridField &grid_field)
{
  if (grid_field.grid.xsize == 0 || grid_field.grid.ysize == 0)
    error("build_terrain_mesh_zemlya: GridField has empty grid dimensions");

  const size_t expected_size = grid_field.grid.xsize * grid_field.grid.ysize;
  if (grid_field.values.size() != expected_size)
  {
    error("build_terrain_mesh_zemlya: GridField values size (" + str(grid_field.values.size()) +
          ") does not match grid dimensions (" + str(expected_size) + ")");
  }

  const double xstep = grid_field.grid.xstep;
  const double ystep = grid_field.grid.ystep;
  const double tolerance = std::max({1.0, std::fabs(xstep), std::fabs(ystep)}) * 1e-12;
  if (std::fabs(xstep - ystep) > tolerance)
  {
    error("build_terrain_mesh_zemlya: Zemlya requires equal x/y spacing, got xstep=" +
          str(xstep) + " and ystep=" + str(ystep));
  }

  terrain_mesher::core::RasterDouble raster(grid_field.grid.xsize, grid_field.grid.ysize);
  raster.set_cell_size(xstep);
  const auto first_sample = grid_field.grid.index_to_point(0);
  raster.set_pos_x(first_sample.x - 0.5 * xstep);
  raster.set_pos_y(first_sample.y - 0.5 * ystep);

  for (size_t row = 0; row < grid_field.grid.ysize; row++)
  {
    const size_t source_row = grid_field.grid.ysize - 1 - row;
    for (size_t col = 0; col < grid_field.grid.xsize; col++)
    {
      raster.value(row, col) = grid_field.values[source_row * grid_field.grid.xsize + col];
    }
  }

  return raster;
}

Mesh build_terrain_mesh_zemlya(const GridField &grid_field, double max_error,size_t smoothing_iterations)
{
  auto raster = gridfield_to_raster(grid_field);
  return terrain_mesher::core::generate_zemlya_mesh(std::move(raster), max_error, smoothing_iterations);
}

py::array_t<double> ray_surface_intersection(const Surface &surface,
                                             const py::array_t<double> &py_ray_origin,
                                             const py::array_t<double> &py_ray_vector)
{

  Vector3D ray_origin(py_ray_origin.at(0), py_ray_origin.at(1), py_ray_origin.at(2));
  Vector3D ray_vector(py_ray_vector.at(0), py_ray_vector.at(1), py_ray_vector.at(2));

  auto intersection = Intersection::ray_surface_intersection(surface, ray_origin, ray_vector);

  return py::array_t<double>(3, &intersection.x);
}

py::array_t<double> ray_multisurface_intersection(const MultiSurface &surface,
                                                  const py::array_t<double> &py_ray_origin,
                                                  const py::array_t<double> &py_ray_vector)
{

  Vector3D ray_origin(py_ray_origin.at(0), py_ray_origin.at(1), py_ray_origin.at(2));
  Vector3D ray_vector(py_ray_vector.at(0), py_ray_vector.at(1), py_ray_vector.at(2));

  auto intersection = Intersection::ray_multisurface_intersection(surface, ray_origin, ray_vector);

  return py::array_t<double>(3, &intersection.x);
}

py::array_t<size_t> statistical_outlier_finder(py::array_t<double> &points, size_t neighbors,
                                               double outlier_margin)
{
  auto points_r = points.unchecked<2>();
  size_t num_points = points_r.shape(0);
  std::vector<Vector3D> pc;
  for (size_t i = 0; i < num_points; i++)
  {
    pc.push_back(Vector3D(points_r(i, 0), points_r(i, 1), points_r(i, 2)));
  }

  auto outliers = PointCloudProcessor::statistical_outlier_finder(pc, neighbors, outlier_margin);
  return py::array_t<size_t>(outliers.size(), outliers.data());
}

} // namespace DTCC_BUILDER

PYBIND11_MODULE(_dtcc_builder, m)
{

#ifdef DTCC_HAVE_TRIANGLE
  constexpr bool have_triangle = true;
#else
  constexpr bool have_triangle = false;
#endif

  m.attr("HAVE_TRIANGLE") = py::bool_(have_triangle);

  m.def(
      "triangulation_backends",
      []()
      {
        std::vector<std::string> backends;
        if (have_triangle)
          backends.emplace_back("triangle");
        if (backends.empty())
          backends.emplace_back("earcut");
        return backends;
      },
      R"pbdoc(
            Report available triangulation backends compiled into the extension.
        )pbdoc");

  py::class_<DTCC_BUILDER::Vector2D>(m, "Vector2D")
      .def(py::init<>())
      .def(
          "__repr__", [](const DTCC_BUILDER::Vector3D &p)
          { return "<Vector3D (" + DTCC_BUILDER::str(p.x) + ", " + DTCC_BUILDER::str(p.y) + ")>"; })
      .def_readonly("x", &DTCC_BUILDER::Vector2D::x)
      .def_readonly("y", &DTCC_BUILDER::Vector2D::y);

  py::class_<DTCC_BUILDER::Vector3D>(m, "Vector3D")
      .def(py::init<>())
      .def("__repr__",
           [](const DTCC_BUILDER::Vector3D &p)
           {
             return "<Vector3D (" + DTCC_BUILDER::str(p.x) + ", " + DTCC_BUILDER::str(p.y) + ", " +
                    DTCC_BUILDER::str(p.z) + ")>";
           })
      .def_readonly("x", &DTCC_BUILDER::Vector3D::x)
      .def_readonly("y", &DTCC_BUILDER::Vector3D::y)
      .def_readonly("z", &DTCC_BUILDER::Vector3D::z);

  // py::class_<DTCC_BUILDER::Vector3D>(m, "Vector3D")
  //     .def(py::init<>())
  //     .def_readonly("x", &DTCC_BUILDER::Vector3D::x)
  //     .def_readonly("y", &DTCC_BUILDER::Vector3D::y)
  //     .def_readonly("z", &DTCC_BUILDER::Vector3D::z);

  py::class_<DTCC_BUILDER::BoundingBox2D>(m, "bounding_box")
      .def(py::init<>())
      .def_readonly("P", &DTCC_BUILDER::BoundingBox2D::P)
      .def_readonly("Q", &DTCC_BUILDER::BoundingBox2D::Q);

  py::class_<DTCC_BUILDER::Polygon>(m, "Polygon")
      .def(py::init<>())
      .def_readonly("vertices", &DTCC_BUILDER::Polygon::vertices)
      .def_readonly("holes", &DTCC_BUILDER::Polygon::holes);

  py::class_<DTCC_BUILDER::GridField>(m, "GridField")
      .def(py::init<>())
      .def_readonly("grid", &DTCC_BUILDER::GridField::grid)
      .def_readonly("values", &DTCC_BUILDER::GridField::values);

  py::class_<DTCC_BUILDER::Grid>(m, "Grid")
      .def(py::init<>())
      .def_readonly("xsize", &DTCC_BUILDER::Grid::xsize)
      .def_readonly("ysize", &DTCC_BUILDER::Grid::ysize)
      .def_readonly("xstep", &DTCC_BUILDER::Grid::xstep)
      .def_readonly("ystep", &DTCC_BUILDER::Grid::ystep);

  py::class_<DTCC_BUILDER::Simplex2D>(m, "Simplex2D")
      .def(py::init<>())
      .def_readonly("v0", &DTCC_BUILDER::Simplex2D::v0)
      .def_readonly("v1", &DTCC_BUILDER::Simplex2D::v1)
      .def_readonly("v2", &DTCC_BUILDER::Simplex2D::v2);

  py::class_<DTCC_BUILDER::Simplex3D>(m, "Simplex3D")
      .def(py::init<>())
      .def_readonly("v0", &DTCC_BUILDER::Simplex3D::v0)
      .def_readonly("v1", &DTCC_BUILDER::Simplex3D::v1)
      .def_readonly("v2", &DTCC_BUILDER::Simplex3D::v2)
      .def_readonly("v3", &DTCC_BUILDER::Simplex3D::v3);

  py::class_<DTCC_BUILDER::Mesh>(m, "Mesh")
      .def(py::init<>())
      .def_readonly("vertices", &DTCC_BUILDER::Mesh::vertices)
      .def_readonly("faces", &DTCC_BUILDER::Mesh::faces)
      .def_readonly("normals", &DTCC_BUILDER::Mesh::normals)
      .def_readonly("markers", &DTCC_BUILDER::Mesh::markers)
      .def(
          "from_cpp",
          [](const DTCC_BUILDER::Mesh &m)
          {
            py::object conv = py::module::import("dtcc_core.builder.model_conversion")
                                  .attr("builder_mesh_to_mesh");
            return conv(m);
          },
          R"pbdoc(
                 Convert this C++ Mesh back into a Python model.Mesh.
             )pbdoc");
  ;

  py::class_<DTCC_BUILDER::VolumeMesh>(m, "VolumeMesh")
      .def(py::init<>())
      .def_readonly("vertices", &DTCC_BUILDER::VolumeMesh::vertices)
      .def_readonly("cells", &DTCC_BUILDER::VolumeMesh::cells)
      .def_readonly("markers", &DTCC_BUILDER::VolumeMesh::markers)
      .def(
          "from_cpp",
          [](const DTCC_BUILDER::VolumeMesh &m)
          {
            py::object conv = py::module::import("dtcc_core.builder.model_conversion")
                                  .attr("builder_volume_mesh_to_volume_mesh");
            return conv(m);
          },
          R"pbdoc(
                 Convert this C++ Volume Mesh back into a Python model.VolumeMesh.
             )pbdoc");

  py::class_<DTCC_BUILDER::Surface>(m, "Surface")
      .def(py::init<>())
      .def_readonly("vertices", &DTCC_BUILDER::Surface::vertices)
      .def_readonly("holes", &DTCC_BUILDER::Surface::holes);

  py::class_<DTCC_BUILDER::MultiSurface>(m, "MultiSurface")
      .def(py::init<>())
      .def_readonly("surfaces", &DTCC_BUILDER::MultiSurface::surfaces);

  m.def("create_polygon", &DTCC_BUILDER::create_polygon, "Create C++ polygon");

  m.def("create_mesh", &DTCC_BUILDER::create_mesh, "Create C++ mesh",
        py::arg("vertices"), py::arg("faces"), py::arg("markers"),
        py::arg("normals") = py::array_t<double>(py::ssize_t{0}));

  m.def("create_volume_mesh", &DTCC_BUILDER::create_volume_mesh, "Create C++ volume mesh");

  m.def("mesh_as_arrays", &DTCC_BUILDER::mesh_as_arrays, "Create C++ mesh");
  m.def("volume_mesh_as_arrays", &DTCC_BUILDER::volume_mesh_as_arrays,
        "Return volume mesh as owning arrays");

  m.def("create_gridfield", &DTCC_BUILDER::create_gridfield, "Create C++ grid field");

  m.def("build_terrain_mesh_zemlya", &DTCC_BUILDER::build_terrain_mesh_zemlya,
        "Build a terrain mesh from a GridField using the Zemlya terrain mesher");

  m.def("extract_building_points", &DTCC_BUILDER::extract_building_points,
        "Compute building points from point cloud");

  m.def("points_in_polygons", &DTCC_BUILDER::points_in_polygons, "Find points inside polygons");

  m.def("build_city_flat_mesh", &DTCC_BUILDER::MeshBuilder::build_city_flat_mesh,
        "build city flat mesh");

  m.def("build_terrain_surface_mesh", &DTCC_BUILDER::MeshBuilder::build_terrain_surface_mesh,
        "build terrain surface mesh");
  m.def("build_terrain_surface_mesh_from_ground_mesh",
        &DTCC_BUILDER::MeshBuilder::build_terrain_surface_mesh_from_ground_mesh,
        "build terrain surface mesh from a prebuilt ground mesh");

  m.def("build_city_surface_mesh_from_terrain_mesh",
        &DTCC_BUILDER::MeshBuilder::build_city_surface_mesh_from_terrain_mesh,
        "build city surface mesh from a prebuilt terrain mesh");

  m.def("merge_meshes", &DTCC_BUILDER::MeshProcessor::merge_meshes,
        "Merge meshes into a single mesh");

  m.def("snap_mesh_vertices", &DTCC_BUILDER::MeshProcessor::snap_vertices, "Snap mesh vertices");

  m.def("create_surface", &DTCC_BUILDER::create_surface, "Create C++ surface");

  m.def("create_multisurface", &DTCC_BUILDER::create_multisurface, "Create C++ multisurface");

  m.def("mesh_surface", &DTCC_BUILDER::MeshBuilder::mesh_surface,
        "Create triangulated mesh from surface");

  m.def("mesh_multisurface", &DTCC_BUILDER::MeshBuilder::mesh_multisurface,
        "Create triangulated mesh from multisurface");
  m.def("mesh_multisurfaces", &DTCC_BUILDER::MeshBuilder::mesh_multisurfaces,
        "Create a lits of triangulated meshes from a list of multisurfaces");

  m.def("ray_surface_intersection", &DTCC_BUILDER::ray_surface_intersection,
        "Compute ray-surface intersection");

  m.def("ray_multisurface_intersection", &DTCC_BUILDER::ray_multisurface_intersection,
        "Compute ray-multisurface intersection");

  m.def("statistical_outlier_finder", &DTCC_BUILDER::statistical_outlier_finder,
        "Find statistical outliers in point cloud");
}
