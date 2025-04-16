#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstddef> // for ptrdiff_t


// Function to write a matrix in CSR format (val, ptr, col)
// to a Matrix Market file in coordinate format.
// 'numRows' is the number of rows and we assume the matrix is square.
void writeMatrixMarketMatrix(const std::string &filename,
  const std::vector<double> &val,
  const std::vector<ptrdiff_t> &ptr,
  const std::vector<ptrdiff_t> &col,
  std::size_t numRows,
  std::size_t numCols)
{
std::ofstream file(filename);
if (!file) {
std::cerr << "Error opening file: " << filename << std::endl;
return;
}

// Write the Matrix Market header for a coordinate real general matrix.
file << "%%MatrixMarket matrix coordinate real general\n";
file << "% Matrix stored in CSR format by writeMatrixMarketMatrix\n";

// Count the number of nonzero entries.
std::size_t nnz = val.size();

// Write matrix dimensions and the number of nonzeros.
file << numRows << " " << numCols << " " << nnz << "\n";

// Loop over each row.
for (std::size_t i = 0; i < numRows; i++) {
// CSR: entries for row i are stored from ptr[i] to ptr[i+1]-1.
for (ptrdiff_t j = ptr[i]; j < ptr[i+1]; j++) {
// Matrix Market format requires 1-indexed rows and columns.
file << (i + 1) << " " << (col[j] + 1) << " " << val[j] << "\n";
}
}

file.close();
}


// Function to write a right hand side vector to a Matrix Market file in array format.
void writeMatrixMarketVector(const std::string &filename,
  const std::vector<double> &rhs)
{
std::ofstream file(filename);
if (!file) {
std::cerr << "Error opening file: " << filename << std::endl;
return;
}

// Write the Matrix Market header for an array format (a dense vector).
file << "%%MatrixMarket matrix array real general\n";
file << "% Vector stored in array format by writeMatrixMarketVector\n";

// The vector is treated as an n-by-1 matrix.
file << rhs.size() << " " << 1 << "\n";

// Write each component of the vector.
for (std::size_t i = 0; i < rhs.size(); i++) {
file << rhs[i] << "\n";
}

file.close();
}

