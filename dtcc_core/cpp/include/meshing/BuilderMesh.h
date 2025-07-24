// include/meshing/LayerHeights.h
#pragma once
#ifndef DTCC_VOLUME_MESH_BUILDER_MESH_H
#define DTCC_VOLUME_MESH_BUILDER_MESH_H

#include <vector>
#include <unordered_set>
#include "model/Mesh.h"

namespace DTCC_BUILDER
{
class BuilderMesh
{
public:
  // References to the original mesh's attributes.

  /// Array of vertices
  std::vector<Vector3D> &vertices;

  /// Array of faces (triangles)
  std::vector<Simplex2D> &faces;

  /// Array of cell markers
  std::vector<int> &markers;

  // Stores vertex to face mapping for ground mesh.
  std::vector<std::unordered_set<int>> vf;

  // Stores vertex to face mapping for ground mesh.
  std::vector<std::unordered_set<int>> ff;

  std::vector<int> face_colors;

  std::vector<int> vertex_colors;

  std::vector<int> face_partitions;

  // BuilderMesh(){}

  BuilderMesh(Mesh &mesh)
      : faces(mesh.faces), vertices(mesh.vertices), markers(mesh.markers),
        face_colors(mesh.faces.size()), vertex_colors(mesh.vertices.size(), mesh.faces.size()),
        face_partitions(mesh.faces.size())
  {
  }

  ~BuilderMesh() {}

  void translate(const Vector3D &t)
    {
      for (auto &v : vertices)
        v += t;
    }
  
  // Build mapping from vertices to faces
  void build_vertex_to_face_mapping()
  {
    const size_t num_faces = faces.size();
    vf.resize(num_faces);
    for (size_t i = 0; i < num_faces; i++)
    {
      vf[faces[i].v0].insert(i);
      vf[faces[i].v1].insert(i);
      vf[faces[i].v2].insert(i);
    }
  }

  /// Build mapping from faces to faces
  void build_face_to_face_mapping()
  {
    const size_t num_faces = faces.size();
    ff.resize(num_faces);
    for (size_t i = 0; i < num_faces; i++)
    {
      for (size_t j = 0; j < 3; j++)
      {
        for (const int &v : vf[faces[i][j]])
        {
          ff[i].insert(v);
        }
      }
    }
  }

  /// Mapping face markers to Vertex markers for ground mesh.
  ///
  /// Each vertex adopts the max marker value from the faces it belongs to.
  ///
  /// Ground Mesh Markers:
  /// -2: Ground
  /// -1: Building halos
  ///  0: Building 0
  ///  1: Building 1
  ///  etc (non-negative integers mark faces inside buildings)
  std::vector<int> face_to_vertex_markers()
  {
    const size_t num_vertices = vertices.size();
    const size_t num_faces = faces.size();

    std::vector<int> vertex_markers(num_vertices, -2);

    if (!markers.size())
    {
      error("Ground mesh has no face Markers. Treating all "
            "faces as "
            "ground");
      return vertex_markers;
    }

    for (size_t f = 0; f < num_faces; f++)
    {
      if (markers[f] < -2)
      {
        error(" Not allowd value for marker:" + str(markers[f]));
      }

      const std::array<size_t, 3> I = {faces[f].v0, faces[f].v1, faces[f].v2};

      vertex_markers[I[0]] = std::max(vertex_markers[I[0]], markers[f]);
      vertex_markers[I[1]] = std::max(vertex_markers[I[1]], markers[f]);
      vertex_markers[I[2]] = std::max(vertex_markers[I[2]], markers[f]);
    }

    return vertex_markers;
  }




  // Reassign face colors to avoid big jumps
  void reassign_face_colors()
  {
    if (ff.empty()){
      build_face_to_face_mapping();
    }

    for (size_t i = 0; i < ff.size(); i++)
    {
      for (const auto &j : ff[i])
      {
        const auto diff = face_colors[i] - face_colors[j];
        if (diff > 1)
        {
          face_colors[i] -= diff - 1;
        }
      }
    }
  }

  // Check face colors (big jumps)
  size_t check_face_colors()
  { 
    if (ff.empty()){
      build_face_to_face_mapping();
    }

    size_t num_big_jumps = 0;
    for (size_t i = 0; i < ff.size(); i++)
    {
      size_t _num_big_jumps = 0;
      for (const auto &j : ff[i])
      {
        const auto jump = face_colors[i] - face_colors[j];
        if (jump > 1)
          _num_big_jumps++;
      }
      if (_num_big_jumps > 0)
        num_big_jumps++;
    }
    double percentage = 100.0 * num_big_jumps / ff.size();
    info("Big jumps: " + str(num_big_jumps) + " / " + str(ff.size()) + " (" +
         str(percentage, 2L) + "%)");

    return num_big_jumps;
  }

  /// Assign vertex colors based on minimum neighbor face color
  void assign_vertex_colors()
  {
    // const size_t num_vertices = mesh.vertices.size();
    const size_t num_faces = faces.size();
    // vertex_colors.resize(num_vertices, layer_heights.size());
    for (size_t i = 0; i < num_faces; i++)
    {
      const auto &face = faces[i];
      for (const auto &j : {face.v0, face.v1, face.v2})
        vertex_colors[j] = std::min(vertex_colors[j], face_colors[i]);
    }
  }

  // Assign face (prism) partitions based on face and vertex colors.
  //
  // Partition 0: 6 vertices, 3 tetrahedrons
  // Partition 1: 7 vertices, 4 tetrahedrons
  // Partition 2: 8 vertices, 5 tetrahedrons
  // Partition 3: 9 vertices, 6 tetrahedrons (reduntant)
  void assign_face_partitions()
  {
    const size_t num_faces = faces.size();
    // face_partitions.resize(num_faces);
    for (size_t i = 0; i < num_faces; i++)
    {
      const auto &face = faces[i];
      const size_t c0 = vertex_colors[face.v0];
      const size_t c1 = vertex_colors[face.v1];
      const size_t c2 = vertex_colors[face.v2];
      face_partitions[i] = 3 * face_colors[i] - (c0 + c1 + c2);
    }
  }
  // Eliminate type 3 partitions which are just two stacked prisms of type 0
  void eliminate_type_3_partitions()
  {
    const size_t num_faces = face_partitions.size();
    for (size_t i = 0; i < num_faces; i++)
    {
      if (face_partitions[i] == 3)
      {
        face_colors[i]--;
        face_partitions[i] = 0;
      }
    }
  }


  // Sort faces by vertex color and index
  void sort_faces_by_vertex_color_and_index()
  {
    for (auto &face : faces)
    {
      size_t c0 = vertex_colors[face.v0];
      size_t c1 = vertex_colors[face.v1];
      size_t c2 = vertex_colors[face.v2];
      if (c0 > c1)
      {
        std::swap(c0, c1);
        std::swap(face.v0, face.v1);
      }
      if (c0 > c2)
      {
        std::swap(c0, c2);
        std::swap(face.v0, face.v2);
      }
      if (c1 > c2)
      {
        std::swap(c1, c2);
        std::swap(face.v1, face.v2);
      }
      // Compare and swap first and second elements
      if (c0 > c1 || (c0 == c1 && face.v0 > face.v1))
      {
        std::swap(c0, c1);
        std::swap(face.v0, face.v1);
      }

      // Compare and swap first and third elements
      if (c0 > c2 || (c0 == c2 && face.v0 > face.v2))
      {
        std::swap(c0, c2);
        std::swap(face.v0, face.v2);
      }

      // Compare and swap second and third elements
      if (c1 > c2 || (c1 == c2 && face.v1 > face.v2))
      {
        std::swap(c1, c2);
        std::swap(face.v1, face.v2);
      }
    }
  }
};

}
#endif