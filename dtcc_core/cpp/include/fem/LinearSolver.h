// Copyright (C) 2025 Anders Logg
// Licensed under the MIT License

#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include "SparseMatrix.h"
// #include <amgcl/backend/builtin.hpp>
// #include <amgcl/make_solver.hpp>
// #include <amgcl/preconditioner/smoothed_aggregation.hpp>
// #include <amgcl/solver/bicgstab.hpp>
#include <vector>

namespace dtcc
{

class LinearSolver
{
public:
  static std::vector<double> solve(const SparseMatrix &A, const std::vector<double> &b)
  {
    // Get matrix size (assume square matrix)
    const size_t N = A.num_rows();

    // Set initial guess
    std::vector<double> x(N, 0.0);

    // Convert SparseMatrix to AMGCL-compatible format (CSR)
    std::vector<int> ptr(N + 1, 0);
    std::vector<int> col;
    std::vector<double> val;
    A.to_csr(ptr, col, val);

    // Define the AMGCL solver (BiCGStab + Smoothed Aggregation)
    // using Solver = amgcl::make_solver<
    //    amgcl::preconditioner::smoothed_aggregation<amgcl::backend::builtin<double>>,
    //    amgcl::solver::bicgstab<amgcl::backend::builtin<double>>>;
    // Solver solver(std::tie(N, ptr, col, val));

    // Solve the system
    // size_t iters;
    // double residual;
    // std::tie(iters, residual) = solver(b, x);

    return x;
  }
};

} // namespace dtcc

#endif // LINEAR_SOLVER_H
