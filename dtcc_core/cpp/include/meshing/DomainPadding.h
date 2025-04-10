#ifndef DTCC_VOLUME_MESH_DOMAIN_PADDING_H
#define DTCC_VOLUME_MESH_DOMAIN_PADDING_H

#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <utility>
#include <vector>

#include "../model/VolumeMesh.h"
#include "../model/Mesh.h"
#include "../Logging.h"

namespace DTCC_BUILDER
{

  class VolumeMeshDomainPadding
  {
    public:

    static Mesh extract_top_mesh(const VolumeMesh &volume_mesh,std::vector<size_t> &new_to_old_index)
    {
      info("Extracting surface from the top boundary of the domain");
      size_t num_vertices = 0;
      // std::vector<size_t> new_to_old_index;
      std::unordered_map<size_t, size_t> old_to_new_index;
      std::vector<Vector3D> top_vertices;
      std::vector<Simplex2D> top_faces;
      
  
      for (size_t i = 0; i < volume_mesh.vertices.size(); i++)
      {
        if (volume_mesh.markers[i] == -3)
        {
          top_vertices.push_back(volume_mesh.vertices[i]);
          old_to_new_index[i] = num_vertices;
          new_to_old_index.push_back(i);
          ++num_vertices;
        }
      }
      info("Top mesh vertices: "+ str(top_vertices.size()));
      std::vector<int> top_markers(top_vertices.size(), -3);
  
      for (size_t i = 0; i < volume_mesh.cells.size(); i++)
      {
        const Simplex3D cell = volume_mesh.cells[i];
  
        const std::array<size_t, 4> cell_vertices = {cell.v0, cell.v1, cell.v2, cell.v3};
  
        std::vector<size_t> face_vertices_top;
        for (const auto v : cell_vertices)
        {
          if (volume_mesh.markers[v] == -3)
            face_vertices_top.push_back(old_to_new_index[v]);
        }
  
        if (face_vertices_top.size() == 3)
        {
          top_faces.push_back(Simplex2D(face_vertices_top[0], face_vertices_top[1], face_vertices_top[2]));
        }
      }
      info("Top mesh faces: "+ str(top_faces.size()));
      Mesh top_mesh;
      top_mesh.vertices = top_vertices;
      top_mesh.faces = top_faces;

      info(top_mesh);
  
      return top_mesh;
    }



  }; 

}

#endif