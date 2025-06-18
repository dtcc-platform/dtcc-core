// Copyright (C) 2023 George Spaias
// Licensed under the MIT License

#ifndef DTCC_BOUNDARY_CONDITIONS_ELASTIC_H
#define DTCC_BOUNDARY_CONDITIONS_ELASTIC_H

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "Geometry.h"
#include "StiffnessMatrix.h"
#include "model/GridField.h"
#include "model/Surface.h"
#include "fem/SparseMatrix.h"

namespace DTCC_BUILDER
{

class BoundaryConditionsElastic
{
  typedef unsigned int uint;

  struct VertexBC {
    bool dofX;
    bool dofY;
    bool dofZ;
    };


public:
  // Vertex Boundary markers:
  // -4 : Neumann vertices
  // -3 : Top Boundary vertices
  // -2 : Ground Boundary vertices
  // -1 : Building Halos Boundary vertices
  const std::vector<int> &vertex_markers;

  const BoundingBox2D _domain_bounds;

  // Boundary values (flag if set and value)
  std::vector<std::pair<VertexBC, Vector3D>> values;

  // Building centroids
  std::vector<Vector3D> building_centroids;

  // Elevation for halo vertices based on min Cell elevation
  std::vector<double> halo_elevations;

  // Constructor
  BoundaryConditionsElastic(const VolumeMesh &volume_mesh,
                            const std::vector<Surface> &building_surfaces, const GridField &dtm,
                            double top_height, bool fix_buildings, bool fix_top,const BoundingBox2D &domain_bounds)
      : _volume_mesh(volume_mesh), _building_surfaces(building_surfaces), _dtm(dtm),
        top_height(top_height), vertex_markers(volume_mesh.markers), _domain_bounds(domain_bounds),
        values(volume_mesh.vertices.size(), std::make_pair(VertexBC{false,false,false}, Vector3D(0, 0, 0))),
        halo_elevations(volume_mesh.vertices.size(), std::numeric_limits<double>::max())
  {

    // Compute boundary values
    compute_boundary_values(fix_buildings, fix_top);
  }

  // Destructor
  ~BoundaryConditionsElastic() {}


  // Compute boundary values
  void compute_boundary_values(bool fix_buildings, bool fix_top)
  {
    info("Computing boundary values");

    // TODO: Check if Search tree has already been built
    //_city.build_search_tree(true,0.0);

    // Compute building centroids
    compute_building_centroids();

    // Compute halo elevations
    compute_halo_elevations();

    const int N = static_cast<int>(building_centroids.size());

    for (size_t i = 0; i < _volume_mesh.vertices.size(); ++i)
    {
      int m = vertex_markers[i];
      auto &bc = values[i];

      // reset to “free” by default
      bc.first = {false, false, false};
      bc.second = Vector3D(0.0, 0.0, 0.0);

      if (m >= 0 && m < N)
      {
        // wall of building k
        bc.first = {true, true, false}; // pin x,y
      }
      else if (m >= N && m < 2 * N && fix_buildings)
      {
        // roof of building k
        bc.first = {true, true, true}; // pin x,y,z
        double dz = building_centroids[m - N].z - _volume_mesh.vertices[i].z;
        bc.second = Vector3D(0.0, 0.0, dz); // snap up to roof plane
      }
      else if (m == -1)
      {
        // halo
        bc.first = {false, false, true};
        double dz = halo_elevations[i] - _volume_mesh.vertices[i].z;
        bc.second = Vector3D(0.0, 0.0, dz);
      }
      else if (m == -2)
      {
        // ground
        bc.first = {false, false, true};
        Vector2D p(_volume_mesh.vertices[i].x, _volume_mesh.vertices[i].y);
        double dz = _dtm(p) - _volume_mesh.vertices[i].z;
        bc.second = Vector3D(0.0, 0.0, dz);
      }
      else if (m == -3 && fix_top)
      {
        // top of domain
        bc.first = {false, false, true};
        double dz = top_height - _volume_mesh.vertices[i].z;
        bc.second = Vector3D(0.0, 0.0, dz);
      }
      else if (m == -4)
      {
        // building columns
        bc.first = {true, true, false};
      }
      else if (m == -5)
      {
        // fully free interior
        bc.first = {false, false, false};
      }
      else
      {
        // unexpected marker → maybe warn or assert
        assert(false && "Unknown vertex_marker");
      }
    }
  }

  // Apply boundary conditions on stiffness matrix for Elastic smoothing
  void apply(dtcc::SparseMatrix &A)
  {
    info("Applying boundary conditions to stiffness Sparse matrix for Elastic Smoothing... (test)");
    const double epsilon = 1e-6;
    const size_t num_vertices = _volume_mesh.vertices.size();
    for (size_t v = 0; v < num_vertices; v++)
    {
      if (values[v].first.dofX){
        A.set_boundary_condition(3*v + 0);
      }
      if (values[v].first.dofY)
      {
        A.set_boundary_condition(3*v + 1);
      }
      if (values[v].first.dofZ)
      {
        A.set_boundary_condition(3*v + 2);
      }

      if (_volume_mesh.vertices[v].x - _domain_bounds.P.x <= epsilon ||  _domain_bounds.Q.x - _volume_mesh.vertices[v].x  <= epsilon )
          A.set_boundary_condition(3*v + 0);
      if (_volume_mesh.vertices[v].y - _domain_bounds.P.y <= epsilon ||  _domain_bounds.Q.y - _volume_mesh.vertices[v].y  <= epsilon )
          A.set_boundary_condition(3*v + 1);
        
      
    }
  }

   // Apply boundary conditions on load vector
   void apply(std::vector<double> &b)
   {
     info("Applying boundary conditions to load vector");
     for (size_t i = 0; i < _volume_mesh.vertices.size(); i++)
     {
       if (values[i].first.dofX || values[i].first.dofX || values[i].first.dofZ){
        b[3*i+0] = values[i].second.x;
        b[3*i+1] = values[i].second.y;
        b[3*i+2] = values[i].second.z;
       }
     }
   }


private:
  

  const VolumeMesh &_volume_mesh;

  const std::vector<Surface> &_building_surfaces;

  const GridField &_dtm;

  const double top_height;

  // Compute building centroids
  void compute_building_centroids()
  {
    building_centroids.resize(_building_surfaces.size());

    for (size_t i = 0; i < _building_surfaces.size(); i++)
    {
      auto centroid = Geometry::surface_centroid(_building_surfaces[i]);
      building_centroids[i] = centroid;
    }
  }

  // Compute halos elevation based on min elevation in containing cells
  void compute_halo_elevations()
  {
    std::array<uint, 4> I = {0};

    for (size_t c = 0; c < _volume_mesh.cells.size(); c++)
    {
      I[0] = _volume_mesh.cells[c].v0;
      I[1] = _volume_mesh.cells[c].v1;
      I[2] = _volume_mesh.cells[c].v2;
      I[3] = _volume_mesh.cells[c].v3;

      double z_min = std::numeric_limits<double>::max();

      for (size_t i = 0; i < 4; i++)
      {
        const Vector2D p(_volume_mesh.vertices[I[i]].x, _volume_mesh.vertices[I[i]].y);
        const double z = _dtm(p);

        z_min = std::min(z_min, z);
      }

      halo_elevations[I[0]] = std::min(halo_elevations[I[0]], z_min);
      halo_elevations[I[1]] = std::min(halo_elevations[I[1]], z_min);
      halo_elevations[I[2]] = std::min(halo_elevations[I[2]], z_min);
      halo_elevations[I[3]] = std::min(halo_elevations[I[3]], z_min);
    }
  }
};

} // namespace DTCC_BUILDER

#endif
