// Copyright (C) 2025 Anders Logg
// Licensed under the MIT License

#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include "BilinearForm.h"
#include "SparseMatrix.h"
#include "model/VolumeMesh.h"

// FIXME: Should be namespace dtcc
using namespace DTCC_BUILDER;

namespace dtcc
{

class Assembler
{
public:
  // Assemble a sparse matrix from a bilinear form and a mesh
  static SparseMatrix assemble(const BilinearForm &a, const VolumeMesh &mesh)
  {
    // FIXME: Incorrect size for elasticity

    // Initialize sparse matrix
    SparseMatrix A(mesh.vertices.size(), mesh.vertices.size());

    // Initialize element dofs and element matrix
    std::vector<size_t> element_dofs(a.dim(), 0);
    std::vector<double> element_matrix(a.dim() * a.dim(), 0.0);

    // Iterate over elements and assemble
    for (const auto &cell : mesh.cells)
    {
      // Get vertex coordinates
      const Vector3D &v0 = mesh.vertices[cell.v0];
      const Vector3D &v1 = mesh.vertices[cell.v1];
      const Vector3D &v2 = mesh.vertices[cell.v2];
      const Vector3D &v3 = mesh.vertices[cell.v3];

      // Compute element dofs and element matrix
      a.compute_element_dofs(element_dofs, cell);
      a.compute_element_matrix(element_matrix, v0, v1, v2, v3);

      // Insert element matrix into sparse matrix
      A.insert_element_matrix(element_matrix, element_dofs);
    }

    return A;
  }
};

} // namespace dtcc

#endif // ASSEMBLER_H
