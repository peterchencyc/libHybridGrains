// SparseMatrixIO.h
//
// Breannan Smith
// Last updated: 10/14/2015

#ifndef SPARSE_MATRIX_IO_H
#define SPARSE_MATRIX_IO_H

#include "scisim/Math/MathDefines.h"

namespace SparseMatrixIO {
// Saves a sparse matrix in a triplet format that can be imported into Matlab
// via:
//   load sparse_matrix.dat
//   H = spconvert( sparse_matrix )
bool writeToMatlabTripletText(const SparseMatrixsc &matrix,
                              const std::string &file_name);
} // namespace SparseMatrixIO

#endif
