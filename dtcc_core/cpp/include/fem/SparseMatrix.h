// Copyright (C) 2025 Anders Logg
// Licensed under the MIT License

#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <unordered_map>
#include <vector>

namespace dtcc
{

class SparseMatrix
{
public:
  SparseMatrix(const size_t rows, const size_t cols) : M(rows), N(cols), matrix(rows) {}

  // Get number of rows
  size_t num_rows() const { return M; }

  // Get number of columns
  size_t num_cols() const { return N; }

  // Insert element matrix
  void insert_element_matrix(const std::vector<double> &element_matrix,
                             const std::vector<size_t> &element_dofs)
  {
    const size_t n = element_dofs.size();
    for (size_t i = 0; i < n; ++i)
    {
      const size_t I = element_dofs[i];
      for (size_t j = 0; j < n; ++j)
      {
        const size_t J = element_dofs[j];
        matrix[I][J] += element_matrix[i * n + j];
      }
    }
  }

  // Set boundary conditions by replacing a rows with the identity vector
  void set_boundary_condition(const size_t row)
  {
    matrix[row].clear();
    matrix[row][row] = 1.0;
  }

  // Set boundary conditions by replacing rows with identity vectors
  void set_boundary_conditions(const std::vector<size_t> &rows)
  {
    for (const size_t row : rows)
    {
      matrix[row].clear();
      matrix[row][row] = 1.0;
    }
  }

  // Convert the sparse matrix to CSR format
  void to_csr(std::vector<int> &ptr, std::vector<int> &col, std::vector<double> &val) const
  {
    size_t row_idx = 0;
    for (const auto &row : matrix)
    {
      ptr[row_idx] = col.size();
      for (const auto &entry : row)
      {
        col.push_back(entry.first);
        val.push_back(entry.second);
      }
      ++row_idx;
    }
    ptr[row_idx] = col.size();
  }

private:
  // Matrix dimensions
  const size_t M, N;

  // Sparse matrix data structure
  std::vector<std::unordered_map<size_t, double>> matrix;
};

} // namespace dtcc

#endif // SPARSE_MATRIX_H
