// Copyright (C) 2025 Anders Logg
// Licensed under the MIT License

#ifndef POISSON_H
#define POISSON_H

#include "BilinearForm.h"
#include "forms/_Poisson.h"
#include <memory>

namespace dtcc
{

// Wrapper for UFC form Poisson
class Poisson : public BilinearForm
{
public:
  // Constructor
  Poisson()
  {
    _poisson_form_0 _ufc_form;
    _cell_integral = std::unique_ptr<ufc::cell_integral>(_ufc_form.create_default_cell_integral());
  }

  // Return local dimension
  std::size_t local_dimension() const { return 4; }

  // Return global dimension
  std::size_t global_dimension(const VolumeMesh &mesh) const { return mesh.vertices.size(); }

  // Compute element dofs
  void compute_element_dofs(std::vector<size_t> &dofs, const Simplex3D &cell) const
  {
    dofs[0] = cell.v0;
    dofs[1] = cell.v1;
    dofs[2] = cell.v2;
    dofs[3] = cell.v3;
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

#endif // POISSON_H
