// Copyright (C) 2023 Dag Wästberg
// Licensed under the MIT License
//
// Modified by Anders Logg 2023

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <pybind11/stl.h>

#include "BuildingProcessor.h"
#include "Intersection.h"
#include "MeshBuilder.h"
#include "VolumeMeshBuilder.h"
#include "MeshProcessor.h"
#include "Smoother.h"
#include "VertexSmoother.h"
#include "model/GridField.h"
#include "model/Mesh.h"
#include "model/Polygon.h"
#include "model/Simplices.h"
#include "model/Vector.h"

namespace py = pybind11;

namespace DTCC_BUILDER
{

  Polygon create_polygon(py::list vertices, py::list holes)
  {
    Polygon poly;
    for (size_t i = 0; i < vertices.size(); i++)
    {
      auto pt = vertices[i].cast<py::tuple>();
      poly.vertices.push_back(
          Vector2D(pt[0].cast<double>(), pt[1].cast<double>()));
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


  Mesh create_mesh(py::array_t<double> vertices,
                   py::array_t<size_t> faces,
                   py::array_t<int> markers)
  {
    Mesh mesh;
    auto verts_r = vertices.unchecked<2>();
    auto faces_r = faces.unchecked<2>();
    auto markers_r = markers.unchecked<1>();
    size_t num_vertices = verts_r.shape(0);
    size_t num_faces = faces_r.shape(0);
    size_t num_markers = markers_r.size();

    for (size_t i = 0; i < num_vertices; i++)
    {
      mesh.vertices.push_back(
          Vector3D(verts_r(i, 0), verts_r(i, 1), verts_r(i, 2)));
    }

    for (size_t i = 0; i < num_faces; i++)
    {
      mesh.faces.push_back(
          Simplex2D(faces_r(i, 0), faces_r(i, 1), faces_r(i, 2)));
    }

    for (size_t i = 0; i < num_markers; i++)
    {
      mesh.markers.push_back(markers_r(i));
    }

    return mesh;
  }

  py::tuple mesh_as_arrays(const Mesh &mesh)
  {
    py::array_t<double> py_vertices(mesh.vertices.size() * 3);
    py::array_t<size_t> py_faces(mesh.faces.size() * 3);
    py::array_t<int> py_markers(mesh.markers.size());
    for (size_t i = 0; i < mesh.vertices.size(); i++)
    {
      py_vertices.mutable_at(i * 3) = mesh.vertices[i].x;
      py_vertices.mutable_at(i * 3 + 1) = mesh.vertices[i].y;
      py_vertices.mutable_at(i * 3 + 2) = mesh.vertices[i].z;
    }
    for (size_t i = 0; i < mesh.faces.size(); i++)
    {
      py_faces.mutable_at(i * 3) = mesh.faces[i].v0;
      py_faces.mutable_at(i * 3 + 1) = mesh.faces[i].v1;
      py_faces.mutable_at(i * 3 + 2) = mesh.faces[i].v2;
    }
    for (size_t i = 0; i < mesh.markers.size(); i++)
    {
      py_markers.mutable_at(i) = mesh.markers[i];
    }
    return py::make_tuple(py_vertices, py_faces, py_markers);
  }

  Surface create_surface(py::array_t<double> vertices, py::list holes)
  {
    Surface surface;
    auto verts_r = vertices.unchecked<2>();
    size_t num_vertices = verts_r.shape(0);
    for (size_t i = 0; i < num_vertices; i++)
    {
      surface.vertices.push_back(
          Vector3D(verts_r(i, 0), verts_r(i, 1), verts_r(i, 2)));
    }
    for (size_t i = 0; i < holes.size(); i++)
    {
      auto hole = holes[i].cast<py::array_t<double>>();
      auto hole_r = hole.unchecked<2>();
      std::vector<Vector3D> hole_vertices;
      size_t num_hole_vertices = hole_r.shape(0);
      for (size_t j = 0; j < num_hole_vertices; j++)
      {
        hole_vertices.push_back(
            Vector3D(hole_r(j, 0), hole_r(j, 1), hole_r(j, 2)));
      }
      surface.holes.push_back(hole_vertices);
    }
    return surface;
  }

  py::list extract_building_points(std::vector<Polygon> &buildings,
                                   const py::array_t<double> &pts,
                                   bool statistical_outlier_remover,
                                   size_t neighbors,
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
        PointCloudProcessor::statistical_outlier_remover(rp, neighbors,
                                                         outlier_margin);
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


  GridField create_gridfield(py::array_t<double> data,
                             py::tuple bounds,
                             size_t xsize,
                             size_t ysize)
  {
    GridField grid_field;
    double px = bounds[0].cast<double>();
    double py = bounds[1].cast<double>();
    double qx = bounds[2].cast<double>();
    double qy = bounds[3].cast<double>();
    auto bbox = BoundingBox2D(Vector2D(px, py), Vector2D(qx, qy));

    grid_field.grid.bounding_box = bbox;
    grid_field.grid.xstep = (qx - px) / xsize;
    grid_field.grid.ystep = (qy - py) / ysize;

    grid_field.grid.xsize = xsize;
    grid_field.grid.ysize = ysize;

    auto data_r = data.unchecked<1>();
    size_t data_count = data_r.size();

    for (size_t i = 0; i < data_count; i++)
    {
      grid_field.values.push_back(data_r(i));
    }

    return grid_field;
  }

  py::array_t<double> ray_surface_intersection(const Surface &surface,
                                               const py::array_t<double> &py_ray_origin,
                                               const py::array_t<double> &py_ray_vector)
  {

    Vector3D ray_origin(py_ray_origin.at(0), py_ray_origin.at(1),
                        py_ray_origin.at(2));
    Vector3D ray_vector(py_ray_vector.at(0), py_ray_vector.at(1),
                        py_ray_vector.at(2));

    auto intersection =
        Intersection::ray_surface_intersection(surface, ray_origin, ray_vector);

    return py::array_t<double>(3, &intersection.x);
  }

  py::array_t<double> ray_multisurface_intersection(const MultiSurface &surface,
                                                    const py::array_t<double> &py_ray_origin,
                                                    const py::array_t<double> &py_ray_vector)
  {

    Vector3D ray_origin(py_ray_origin.at(0), py_ray_origin.at(1),
                        py_ray_origin.at(2));
    Vector3D ray_vector(py_ray_vector.at(0), py_ray_vector.at(1),
                        py_ray_vector.at(2));

    auto intersection =
        Intersection::ray_multisurface_intersection(surface, ray_origin, ray_vector);

    return py::array_t<double>(3, &intersection.x);
  }

py::array_t<size_t> statistical_outlier_finder(py::array_t<double> &points, size_t neighbors, double outlier_margin)
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


  py::class_<DTCC_BUILDER::Vector2D>(m, "Vector2D")
      .def(py::init<>())
      .def("__repr__",
           [](const DTCC_BUILDER::Vector3D &p)
           {
             return "<Vector3D (" + DTCC_BUILDER::str(p.x) + ", " +
                    DTCC_BUILDER::str(p.y) + ")>";
           })
      .def_readonly("x", &DTCC_BUILDER::Vector2D::x)
      .def_readonly("y", &DTCC_BUILDER::Vector2D::y);

  py::class_<DTCC_BUILDER::Vector3D>(m, "Vector3D")
      .def(py::init<>())
      .def("__repr__",
           [](const DTCC_BUILDER::Vector3D &p)
           {
             return "<Vector3D (" + DTCC_BUILDER::str(p.x) + ", " +
                    DTCC_BUILDER::str(p.y) + ", " + DTCC_BUILDER::str(p.z) +
                    ")>";
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
      .def_readonly("markers", &DTCC_BUILDER::Mesh::markers);

  py::class_<DTCC_BUILDER::VolumeMesh>(m, "VolumeMesh")
      .def(py::init<>())
      .def_readonly("num_layers", &DTCC_BUILDER::VolumeMesh::num_layers)
      .def_readonly("vertices", &DTCC_BUILDER::VolumeMesh::vertices)
      .def_readonly("cells", &DTCC_BUILDER::VolumeMesh::cells)
      .def_readonly("markers", &DTCC_BUILDER::VolumeMesh::markers);

  py::class_<DTCC_BUILDER::Surface>(m, "Surface")
      .def(py::init<>())
      .def_readonly("vertices", &DTCC_BUILDER::Surface::vertices)
      .def_readonly("holes", &DTCC_BUILDER::Surface::holes);

  py::class_<DTCC_BUILDER::MultiSurface>(m, "MultiSurface")
      .def(py::init<>())
      .def_readonly("surfaces", &DTCC_BUILDER::MultiSurface::surfaces);

  m.def("create_polygon", &DTCC_BUILDER::create_polygon, "Create C++ polygon");

  m.def("create_mesh", &DTCC_BUILDER::create_mesh, "Create C++ mesh");

  m.def("mesh_as_arrays", &DTCC_BUILDER::mesh_as_arrays, "Create C++ mesh");

  m.def("create_gridfield", &DTCC_BUILDER::create_gridfield,
        "Create C++ grid field");


  m.def("extract_building_points", &DTCC_BUILDER::extract_building_points,
        "Compute building points from point cloud");

  m.def("smooth_field", &DTCC_BUILDER::VertexSmoother::smooth_field,
        "Smooth grid field");


  // m.def("build_mesh", &DTCC_BUILDER::MeshBuilder::build_mesh,
  //       "build mesh for city, returning a list of meshes");

  m.def("build_ground_mesh", &DTCC_BUILDER::MeshBuilder::build_ground_mesh,
        "build ground mesh");

  m.def("build_terrain_mesh", &DTCC_BUILDER::MeshBuilder::build_terrain_mesh,
        "build terrain mesh");

  m.def("build_city_surface_mesh",
        &DTCC_BUILDER::MeshBuilder::build_city_surface_mesh,
        "build city surface mesh");

  m.def("layer_ground_mesh", &DTCC_BUILDER::MeshBuilder::layer_ground_mesh,
        "Layer ground mesh");

  m.def("smooth_volume_mesh", &DTCC_BUILDER::Smoother::smooth_volume_mesh,
        "Smooth volume mesh");

  m.def("trim_volume_mesh", &DTCC_BUILDER::MeshBuilder::trim_volume_mesh,
        "Trim volume mesh by removing cells inside buildings");

  // m.def("extrude_footprint", &DTCC_BUILDER::MeshBuilder::extrude_footprint,
  //       "Extrude footprint to a mesh");

  m.def("compute_boundary_mesh",
        &DTCC_BUILDER::MeshProcessor::compute_boundary_mesh,
        "Compute boundary mesh from volume mesh");

  m.def("compute_open_mesh", &DTCC_BUILDER::MeshProcessor::compute_open_mesh,
        "Compute open mesh from boundary, excluding top and sides");

  m.def("merge_meshes", &DTCC_BUILDER::MeshProcessor::merge_meshes,
        "Merge meshes into a single mesh");

  m.def("create_surface", &DTCC_BUILDER::create_surface, "Create C++ surface");

  m.def("create_multisurface", &DTCC_BUILDER::create_multisurface,
        "Create C++ multisurface");

  m.def("mesh_surface", &DTCC_BUILDER::MeshBuilder::mesh_surface,
        "Create triangulated mesh from surface");

  m.def("mesh_multisurface", &DTCC_BUILDER::MeshBuilder::mesh_multisurface,
        "Create triangulated mesh from multisurface");
  m.def("mesh_multisurfaces", &DTCC_BUILDER::MeshBuilder::mesh_multisurfaces,
        "Create a lits of triangulated meshes from a list of multisurfaces");

  m.def("ray_surface_intersection", &DTCC_BUILDER::ray_surface_intersection,
        "Compute ray-surface intersection");

  m.def("ray_multisurface_intersection", &DTCC_BUILDER::ray_multisurface_intersection, "Compute ray-multisurface intersection");

  m.def("statistical_outlier_finder", &DTCC_BUILDER::statistical_outlier_finder, "Find statistical outliers in point cloud");

  py::class_<DTCC_BUILDER::VolumeMeshBuilder>(m, "VolumeMeshBuilder")
      .def(py::init<const std::vector<DTCC_BUILDER::Surface> &, const DTCC_BUILDER::GridField &,
                    DTCC_BUILDER::Mesh &, double>(),
           py::arg("buildings"), py::arg("dem"), py::arg("ground_mesh"),
           py::arg("domain_height"),
           "Constructor for VolumeMeshBuilder taking city, dem, ground_mesh, "
           "and domain_height as arguments")
      .def("build", &DTCC_BUILDER::VolumeMeshBuilder::build,
           "Layers the ground mesh and returns a VolumeMesh")
      // Expose public variables directly
      .def_readwrite("domain_height",
                     &DTCC_BUILDER::VolumeMeshBuilder::domain_height)
      .def_readwrite("top_height", &DTCC_BUILDER::VolumeMeshBuilder::top_height)
      // If you need to expose std::vectors or similar, pybind11/stl.h header
      // takes care of this. For custom types like City, GridField, Mesh, ensure
      // you've also provided bindings for them.
      ;

}