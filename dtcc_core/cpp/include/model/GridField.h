// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_GRID_FIELD_H
#define DTCC_GRID_FIELD_H

#include <cassert>

#include "Geometry.h"
#include "Grid.h"
#include "Logging.h"

namespace DTCC_BUILDER
{

/// GridField represents a scalar field on a uniform grid.
/// The field can be efficiently evaluated at arbitrary points inside
/// the grid domain. The value is computed by bilinear interpolation.
class GridField : public Printable
{
public:
  /// The grid
  Grid grid{};

  /// Array of values (vertex values)
  std::vector<double> values{};

  /// Create empty field
  GridField() = default;
  virtual ~GridField() {} // make the destructor virtual

  /// Create zero field on given grid
  ///
  /// @param grid The grid
  explicit GridField(const Grid &_grid) : grid(_grid)
  {
    // Initialize values to zero
    values.resize(_grid.num_vertices());
    std::fill(values.begin(), values.end(), 0.0);
  }

  /// Evaluate field at given point
  ///
  /// @param p The point
  /// @return Value at point
  double operator()(const Vector2D &p) const
  {
    // Map point to cell
    size_t i{};
    double x{}, y{};
    grid.point_to_cell(p, i, x, y);

    // Get value at lower left corner of cell
    const double v00 = values[i];

    // Get indices for other corners of cell
    const size_t i10 = i + 1;
    const size_t i01 = i + grid.xsize;
    const size_t i11 = i + grid.xsize + 1;
    double v10{}, v01{}, v11{};

    // Get value at lower right corner of cell
    if (i10 % grid.xsize != 0 && i10 < values.size())
      v10 = values[i10];
    else
      v10 = v00;

    // Get value at upper left corner of cell
    if (i01 < values.size())
      v01 = values[i01];
    else
      v01 = v00;

    // Get value at upper right corner of cell
    if (i11 % grid.xsize != 0 && i11 < values.size())
      v11 = values[i11];
    else
      v11 = v10;

    // Compute value by bilinear interpolation
    return DTCC_BUILDER::Grid::interpolate(x, y, v00, v10, v01, v11);
  }

  /// Evaluate field at given 3D point (using only x and y)
  double operator()(const Vector3D &p) const
  {
    Vector2D _p(p.x, p.y);
    return (*this)(_p);
  }

  double nearest(const Vector2D &p) const
  {
    size_t i{};
    double x{}, y{};
    grid.point_to_cell(p, i, x, y);
    return values[i];
  }

  /// interpolate given field at vertices.
  ///
  /// @param field The field to be interpolated
  void interpolate(const GridField &field)
  {
    // Iterate over vertices and evaluate field
    for (size_t i = 0; i < grid.num_vertices(); i++)
      values[i] = field(grid.index_to_point(i));
  }

  /// Compute minimal vertex value.
  ///
  /// @return Minimal vertex value
  double min() const { return *std::min_element(values.begin(), values.end()); }

  /// Compute maximal of vertex value.
  ///
  /// @return Maximal vertex value
  double max() const { return *std::max_element(values.begin(), values.end()); }

  /// Compute mean vertex value
  ///
  /// @return mean vertex value
  double mean() const
  {
    double mean = 0.0;
    for (const auto &value : values)
      mean += value;
    return mean / static_cast<double>(values.size());
  }

  /// Pretty-print
  std::string __str__() const { return "2D field on " + str(grid); }
};

} // namespace DTCC_BUILDER

#endif
