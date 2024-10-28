// Copyright (C) 2019 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_MESH_H
#define DTCC_MESH_H

#include <vector>

#include "Logging.h"
#include "model/Simplices.h"
#include "model/Vector.h"

namespace DTCC_BUILDER
{

  /// Mesh represents a triangular mesh in 3D

  class Mesh : public Printable
  {
  public:
    /// Array of vertices
    std::vector<Vector3D> vertices{};

    /// Array of faces (triangles)
    std::vector<Simplex2D> faces{};

    /// Array of normals
    std::vector<Vector3D> normals{};

    /// Array of cell markers
    std::vector<int> markers{};

    Mesh() = default;
    virtual ~Mesh() {} // make the destructor virtual

    /// Compute midpoint of cell
    Vector3D mid_point(size_t cell_index) const
    {
      Vector3D c{};
      c += Vector3D(vertices[faces[cell_index].v0]);
      c += Vector3D(vertices[faces[cell_index].v1]);
      c += Vector3D(vertices[faces[cell_index].v2]);
      c /= 3.0;
      return c;
    }

    Vector3D centroid() const
    {
      Vector3D c{};
      for (auto &v : vertices)
      {
        c += Vector3D(v);
      }
      c /= vertices.size();
      return c;
    }

    std::vector<std::vector<Vector3D>> triangles() const
    {
      std::vector<std::vector<Vector3D>> triangles;
      for (auto &face : faces)
      {
        std::vector<Vector3D> triangle;
        triangle.push_back(vertices[face.v0]);
        triangle.push_back(vertices[face.v1]);
        triangle.push_back(vertices[face.v2]);
        triangles.push_back(triangle);
      }
      return triangles;
    }

    void calc_normals()
    {
      normals.resize(faces.size());
      for (size_t i = 0; i < faces.size(); i++)
      {
        auto v0 = vertices[faces[i].v0];
        auto v1 = vertices[faces[i].v1];
        auto v2 = vertices[faces[i].v2];
        normals[i] = (v1 - v0).cross(v2 - v0);
        normals[i].normalize();
      }
    }

    /// sort vertices of each face
    void sort_vertices(bool ccw = true)
    {
      for (auto &face : faces)
      {
        // Get the vertices of the face
        auto v0 = vertices[face.v0];
        auto v1 = vertices[face.v1];
        auto v2 = vertices[face.v2];
        auto normal = (v1 - v0).cross(v2 - v0);
        // Check if the normal vector is pointing in the desired direction
        bool swap = ccw ? normal.z < 0 : normal.z > 0;
        if (swap) // Assuming z-axis is up
        {
          // Swap v1 and v2 to change the orientation
          std::swap(face.v1, face.v2);
        }
      }
    }

    void normalize_normal_direction(bool outwards = true)
    {
      for (auto &face : faces)
      {
        auto centroid = this->centroid();
        auto v0 = vertices[face.v0];
        auto v1 = vertices[face.v1];
        auto v2 = vertices[face.v2];
        auto normal = (v1 - v0).cross(v2 - v0);
        // Calculate the vector from the centroid to the face
        auto to_face = (v0 + v1 + v2) / 3.0 - centroid;
        bool swap = outwards ? normal.dot(to_face) < 0 : normal.dot(to_face) > 0;
        if (swap)
        {
          // Swap v1 and v2 to change the orientation
          std::swap(face.v1, face.v2);
        }
      }
    }

    /// Pretty-print
    std::string __str__() const override
    {
      return "Mesh with " + str(vertices.size()) + " vertices and " +
             str(faces.size()) + " faces";
    }
  };

} // namespace DTCC_BUILDER

#endif
