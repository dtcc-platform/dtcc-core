// Copyright (C) 2025 Anders Logg
// Licensed under the MIT License

#ifndef ELASTICITY_H
#define ELASTICITY_H

#include "BilinearForm.h"
#include "forms/_Elasticity.h"
#include <memory>

namespace dtcc
{

// Wrapper for UFC form Elasticity
class Elasticity : public BilinearForm
{
public:
  // Constructor
  Elasticity()
  {
    _elasticity_form_0 _ufc_form;
    _cell_integral = std::unique_ptr<ufc::cell_integral>(_ufc_form.create_default_cell_integral());
  }

  // Return local dimension
  std::size_t local_dimension() const { return 12; }

  // Return global dimension
  std::size_t global_dimension(const VolumeMesh &mesh) const { return 3 * mesh.vertices.size(); }

  // Compute element dofs
  void compute_element_dofs(std::vector<size_t> &dofs, const Simplex3D &cell) const
  {
    dofs[0] = 3 * cell.v0;
    dofs[1] = 3 * cell.v1;
    dofs[2] = 3 * cell.v2;
    dofs[3] = 3 * cell.v3;
    dofs[4] = 3 * cell.v0 + 1;
    dofs[5] = 3 * cell.v1 + 1;
    dofs[6] = 3 * cell.v2 + 1;
    dofs[7] = 3 * cell.v3 + 1;
    dofs[8] = 3 * cell.v0 + 2;
    dofs[9] = 3 * cell.v1 + 2;
    dofs[10] = 3 * cell.v2 + 2;
    dofs[11] = 3 * cell.v3 + 2;
  }

  // Compute element matrix
  void compute_element_matrix(std::vector<double> &element_matrix, const Vector3D &v0,
                              const Vector3D &v1, const Vector3D &v2, const Vector3D &v3) const
  {
    update_vertex_coordinates(v0, v1, v2, v3);
    _cell_integral->tabulate_tensor(element_matrix.data(), 0, vertex_coordinates.data(), 0);
  }

private:
  std::unique_ptr<ufc::cell_integral> _cell_integral;
};

} // namespace dtcc

#endif // ELASTICITY_H
