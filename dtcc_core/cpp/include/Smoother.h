// Copyright (C) 2023 George Spaias
// Licensed under the MIT License

#ifndef DTCC_SMOOTHER_H
#define DTCC_SMOOTHER_H

// Testing include of new FEM classes
#include "fem/Assembler.h"
#include "fem/Elasticity.h"
#include "fem/LinearSolver.h"
#include "fem/Poisson.h"
#include "fem/SparseMatrix.h"

#include "BoundaryConditions.h"
#include "StiffnessMatrix.h"
#include "Timer.h"
#include "model/Mesh.h"
#include "model/Surface.h"

namespace DTCC_BUILDER
{

class Smoother
{
  typedef unsigned int uint;

public:
  // Smooth mesh using Laplacian smoothing
  static VolumeMesh smooth_volume_mesh(const VolumeMesh &volume_mesh,
                                       const std::vector<Surface> &building_surfaces,
                                       const GridField &dem, double top_height, bool fix_buildings,
                                       bool fix_top, size_t max_iterations,
                                       double relative_tolerance)

  {
    info("Smoothing volume mesh...");
    info(volume_mesh.__str__());

    // Compute (local) stiffness matrices
    StiffnessMatrix AK(volume_mesh);

    // Create solution vector and load vector
    std::vector<double> u(volume_mesh.vertices.size(), 0);
    std::vector<double> b(volume_mesh.vertices.size(), 0);

    // Apply boundary conditions
    BoundaryConditions bc(volume_mesh, building_surfaces, dem, top_height, fix_buildings, fix_top);
    bc.apply(AK);
    bc.apply(b);

    // Set initial guess
    //    if (!fix_buildings)
    //      set_initial_guess(u, volume_mesh, dem, top_height, bc);
    //   else
    u = b;

    // Solve linear system
    solve_unassembled_gauss_seidel(volume_mesh, AK, b, u, max_iterations, relative_tolerance);

    // Update mesh coordinates
    VolumeMesh _volume_mesh{volume_mesh};
    for (std::size_t i = 0; i < volume_mesh.vertices.size(); i++)
      _volume_mesh.vertices[i].z += u[i];

    return _volume_mesh;
  }

private:
  // Solve linear system using unassembled Gauss-Seidel iterations
  static void solve_unassembled_gauss_seidel(const VolumeMesh &volume_mesh, StiffnessMatrix &AK,
                                             std::vector<double> &b, std::vector<double> &u,
                                             const size_t max_iterations,
                                             const double relative_tolerance)
  {
    info("Solving linear system using unassembled Gauss-Seidel");

    // Sum of non-diagonal elements
    std::vector<double> C(volume_mesh.vertices.size());

    // Vertex indices of current cell
    std::array<uint, 4> I = {0};

    // Compute the number of cells that each vertex belongs
    std::vector<uint> vertex_degrees(volume_mesh.vertices.size());
    std::vector<uint> _vertex_degrees(volume_mesh.vertices.size());
    compute_vertex_degrees(vertex_degrees, volume_mesh);

    // Gauss-Seidel iterations
    size_t iterations;
    double residual;
    for (iterations = 0; iterations < max_iterations; iterations++)
    {
      // Initialize right-hand side and residual
      C = b;
      residual = 0;

      // Initialize vertex degrees
      _vertex_degrees = vertex_degrees;

      // Iterate over cells
      for (size_t c = 0; c < volume_mesh.cells.size(); c++)
      {
        // Get vertex indices
        I[0] = volume_mesh.cells[c].v0;
        I[1] = volume_mesh.cells[c].v1;
        I[2] = volume_mesh.cells[c].v2;
        I[3] = volume_mesh.cells[c].v3;

        // Gauss-Seidel update
        for (uint8_t i = 0; i < 4; i++)
        {
          // Update right-hand side
          C[I[i]] -= AK._data[c * 16 + i * 4 + (i + 1) % 4] * u[I[(i + 1) % 4]] +
                     AK._data[c * 16 + i * 4 + (i + 2) % 4] * u[I[(i + 2) % 4]] +
                     AK._data[c * 16 + i * 4 + (i + 3) % 4] * u[I[(i + 3) % 4]];

          // Divide by diagonal when fully updated
          _vertex_degrees[I[i]]--;
          if (_vertex_degrees[I[i]] == 0)
          {
            double res = u[I[i]];
            u[I[i]] = C[I[i]] / AK.diagonal[I[i]];
            res = std::abs(res - u[I[i]]);
            residual = std::max(residual, res);
          }
        }
      }

      // Check convergence
      if (residual < relative_tolerance)
        break;
    }

    // Check convergence
    if (iterations == max_iterations)
    {
      error("Failed to converge in " + str(max_iterations) + " iterations with residual " +
            str(residual));
    }
    else
    {
      info("Converged in " + str(iterations) + "/" + str(max_iterations) +
           " iterations with residual " + str(residual));
    }
  }

  // Set initial guess for solution vector
  static void set_initial_guess(std::vector<double> &u, const VolumeMesh &volume_mesh,
                                const GridField &dem, double top_height, BoundaryConditions &bc)
  {
    info("Setting initial guess for solution vector");

    for (size_t i = 0; i < volume_mesh.vertices.size(); i++)
    {
      if (bc.vertex_markers[i] == -4)
      {
        const Vector2D p(volume_mesh.vertices[i].x, volume_mesh.vertices[i].y);
        u[i] = dem(p) * (1 - volume_mesh.vertices[i].z / top_height);
      }
      else
        u[i] = 0.0;
    }
  }

  static std::vector<double>
  get_adjusted_building_heights(const VolumeMesh &volume_mesh,
                                const std::vector<Surface> &building_surfaces)
  {
    std::vector<double> adj_heights(building_surfaces.size(), 0.0);

    for (size_t i = 0; i < volume_mesh.vertices.size(); i++)
    {
      int marker = volume_mesh.markers[i];
      if (marker >= 0)
      {
        adj_heights[marker] = volume_mesh.vertices[i].z;
      }
    }
    for (size_t i = 0; i < adj_heights.size(); i++)
    {
      std::cout << i << ") Building max height: " << building_surfaces[i].max_height()
                << " adj: " << adj_heights[i] << std::endl;
    }

    return adj_heights;
  }

  // Compute the number of cells to which each vertex belongs
  static void compute_vertex_degrees(std::vector<uint> &vertex_degrees,
                                     const VolumeMesh &volume_mesh)
  {
    for (size_t c = 0; c < volume_mesh.cells.size(); c++)
    {
      vertex_degrees[volume_mesh.cells[c].v0]++;
      vertex_degrees[volume_mesh.cells[c].v1]++;
      vertex_degrees[volume_mesh.cells[c].v2]++;
      vertex_degrees[volume_mesh.cells[c].v3]++;
    }
  }
};

} // namespace DTCC_BUILDER

#endif // DTCC_LAPLACIAN_SMOOTHER_NEW_H
