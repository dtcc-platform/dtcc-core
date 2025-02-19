// Copyright (C) 2025 Anders Logg
// Licensed under the MIT License

#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include "SparseMatrix.h"
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/gmres.hpp>

#include <amgcl/profiler.hpp>

#include <vector>

namespace dtcc
{

  // Checks if the stiffness matrix in CSR form is "solvable"
  // by verifying that each row has at least one nonzero entry and a non-nearly-zero diagonal.
  bool isSolvableCSR(const std::vector<ptrdiff_t>& ptr,
                    const std::vector<ptrdiff_t>& col,
                    const std::vector<double>& val,
                    double tol = 1e-12)
  {
      // Ensure that the CSR pointer vector is not empty.
      if (ptr.empty()) {
          std::cerr << "CSR pointer vector is empty." << std::endl;
          return false;
      }
      
      // The number of rows is ptr.size() - 1.
      size_t n_rows = ptr.size() - 1;
      
      // Check each row.
      for (size_t i = 0; i < n_rows; ++i)
      {
          bool nonZeroRow = false;  // Does the row contain any nonzero entry?
          bool foundDiag = false;   // Is the diagonal entry present?
          double diagValue = 0.0;
          
          // Iterate over the entries in row i.
          for (ptrdiff_t j = ptr[i]; j < ptr[i+1]; ++j)
          {
              double currentValue = val[j];
              // Check for nonzero entry in the row.
              if (std::fabs(currentValue) > tol) {
                  nonZeroRow = true;
              }
              
              // Check if this entry is on the diagonal.
              if (col[j] == static_cast<ptrdiff_t>(i)) {
                  foundDiag = true;
                  diagValue = currentValue;
              }
          }
          
          // If the row has no nonzero entries, it's a sign of singularity.
          if (!nonZeroRow)
          {
              std::cerr << "Row " << i << " is completely zero." << std::endl;
              return false;
          }
          
          // If the diagonal entry is missing, we cannot guarantee stability.
          if (!foundDiag)
          {
              std::cerr << "Row " << i << " is missing its diagonal entry." << std::endl;
              return false;
          }
          
          // If the diagonal is nearly zero, the row may be ill-conditioned.
          if (std::fabs(diagValue) < tol)
          {
              std::cerr << "Row " << i << " has a nearly zero diagonal value (" 
                        << diagValue << ")." << std::endl;
              return false;
          }
      }
      
      // If all rows pass the checks, the matrix is (likely) solvable.
      return true;
  }

class LinearSolver
{
public:
  static std::vector<double> solve(const SparseMatrix &A, 
                                  const std::vector<double> &b,
                                  size_t max_iterations,
                                  double relative_tolerance)
  {

    // The profiler:
    amgcl::profiler<> prof("Smoothing");
    // Get matrix size (assume square matrix)
   
    ptrdiff_t rows = A.num_rows();
    // ptrdiff_t cols = A.num_cols();
    // std::vector<ptrdiff_t> ptr, col;
    std::vector<double> val, rhs;
    // Set initial guess
    std::vector<double> x(rows, 0.0);

    // Convert SparseMatrix to AMGCL-compatible format (CSR)
    std::vector<ptrdiff_t> ptr(rows + 1, 0);
    std::vector<ptrdiff_t> col;

    prof.tic("To CSR");
    A.to_csr(ptr, col, val);
    prof.toc("To CSR");
    prof.tic("Check CSR");
    if (!dtcc::isSolvableCSR(ptr, col, val)){
      std::cout << "Stiffness Matrix in CSR format is not solvable!" << std::endl;
    }
    prof.toc("Check CSR");

    // Compose the solver type
    //   the solver backend:
    typedef amgcl::backend::builtin<double> Backend;
    // //   the preconditioner backend:
    // #ifdef MIXED_PRECISION
    //     typedef amgcl::backend::builtin<float> PBackend;
    // #else
    //     typedef amgcl::backend::builtin<double> PBackend;
    // #endif
        
    typedef amgcl::make_solver<
            amgcl::amg<
                Backend,
                amgcl::coarsening::smoothed_aggregation,
                amgcl::relaxation::spai0
                >,
            amgcl::solver::gmres<Backend>
            > Solver;

    Solver::params prm; 
    prm.solver.M = max_iterations;
    prm.solver.tol = relative_tolerance;
    
    // Define the AMGCL solver (BiCGStab + Smoothed Aggregation)
    // using Solver = amgcl::make_solver<
    //    amgcl::preconditioner::dummy<amgcl::backend::builtin<double>>,
    //    amgcl::solver::bicgstab<amgcl::backend::builtin<double>>>;

    prof.tic("setup");
    Solver solver(std::tie(rows, ptr, col, val),prm);
    prof.toc("setup");
    // // Solve the system
    size_t iters;
    double residual;
    prof.tic("solving");
    std::tie(iters, residual) = solver(b, x);
    prof.toc("solving");
    // Output the number of iterations, the relative error,
    // and the profiling data:
    std::cout << "Smoothing Completed in:" << std::endl;
    std::cout << "Iters: " << iters << std::endl
              << "Residual: " << residual << std::endl
              << prof << std::endl;
    

    return x;
  }

  static void solve(const SparseMatrix &A, 
                  const std::vector<double> &b, 
                  std::vector<double> &u,
                  size_t max_iterations,
                  double relative_tolerance)
  {

    // The profiler:
    amgcl::profiler<> prof("Smoothing");
    // Get matrix size (assume square matrix)
   
    ptrdiff_t rows = A.num_rows();
    // ptrdiff_t cols = A.num_cols();
    // std::vector<ptrdiff_t> ptr, col;
    std::vector<double> val, rhs;
    
    // Convert SparseMatrix to AMGCL-compatible format (CSR)
    std::vector<ptrdiff_t> ptr(rows + 1, 0);
    std::vector<ptrdiff_t> col;

    prof.tic("To CSR");
    A.to_csr(ptr, col, val);
    prof.toc("To CSR");
    prof.tic("Check CSR");
    if (!dtcc::isSolvableCSR(ptr, col, val)){
      std::cout << "Stiffness Matrix in CSR format is Solvable!" <<std::endl;
    }
    prof.toc("Check CSR");

    // Compose the solver type
    //   the solver backend:
    typedef amgcl::backend::builtin<double> Backend;
    // //   the preconditioner backend:
    // #ifdef MIXED_PRECISION
    //     typedef amgcl::backend::builtin<float> PBackend;
    // #else
    //     typedef amgcl::backend::builtin<double> PBackend;
    // #endif
        
    typedef amgcl::make_solver<
            amgcl::amg<
                Backend,
                amgcl::coarsening::smoothed_aggregation,
                amgcl::relaxation::spai0
                >,
            amgcl::solver::gmres<Backend>
            > Solver;

    Solver::params prm; 
    prm.solver.M = max_iterations;
    prm.solver.tol = relative_tolerance;

    prof.tic("setup");
    Solver solver(std::tie(rows, ptr, col, val),prm);
    prof.toc("setup");
    // // Solve the system
    size_t iters;
    double residual;
    prof.tic("solving");
    std::tie(iters, residual) = solver(b, u);
    prof.toc("solving");
    // Output the number of iterations, the relative error,
    // and the profiling data:
    std::cout << "Smoothing Completed in:" << std::endl;
    std::cout << "Iters: " << iters << std::endl
              << "Residual: " << residual << std::endl
              << prof << std::endl;
    
  }
};

} // namespace dtcc

#endif // LINEAR_SOLVER_H
