// Copyright (C) 2025 Anders Logg
// Licensed under the MIT License

#ifndef BILINEAR_FORM_H
#define BILINEAR_FORM_H

#include "model/Simplices.h"
#include "model/Vector.h"
#include <vector>

// FIXME
using namespace DTCC_BUILDER;

namespace dtcc
{

class BilinearForm
{
public:
  // Constructor (initialize vertex coordinates)
  BilinearForm() : vertex_coordinates(12, 0.0) {}

  // Return local dimension (overridden in derived classes)
  virtual std::size_t local_dimension() const = 0;

  // Return global dimension (overridden in derived classes)
  virtual std::size_t global_dimension(const VolumeMesh &mesh) const = 0;

  // Compute element dofs (overridden in derived classes)
  virtual void compute_element_dofs(std::vector<size_t> &dofs, const Simplex3D &cell) const = 0;

  // Compute element matrix (overridden in derived classes)
  virtual void compute_element_matrix(std::vector<double> &element_matrix, const Vector3D &v0,
                                      const Vector3D &v1, const Vector3D &v2,
                                      const Vector3D &v3) const = 0;

protected:
  // Update vertex coordinates
  void update_vertex_coordinates(const Vector3D &v0, const Vector3D &v1, const Vector3D &v2,
                                 const Vector3D &v3) const
  {
    vertex_coordinates[0] = v0.x;
    vertex_coordinates[1] = v0.y;
    vertex_coordinates[2] = v0.z;
    vertex_coordinates[3] = v1.x;
    vertex_coordinates[4] = v1.y;
    vertex_coordinates[5] = v1.z;
    vertex_coordinates[6] = v2.x;
    vertex_coordinates[7] = v2.y;
    vertex_coordinates[8] = v2.z;
    vertex_coordinates[9] = v3.x;
    vertex_coordinates[10] = v3.y;
    vertex_coordinates[11] = v3.z;
  }

  // Vector of flattened vertex coordinates
  mutable std::vector<double> vertex_coordinates;
};

} // namespace dtcc

#endif // BILINEAR_FORM_H
