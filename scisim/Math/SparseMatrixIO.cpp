// SparseMatrixIO.cpp
//
// Breannan Smith
// Last updated: 10/15/2015

#include "SparseMatrixIO.h"

#include <fstream>

bool SparseMatrixIO::writeToMatlabTripletText(const SparseMatrixsc &matrix,
                                              const std::string &file_name) {
  std::ofstream output_file{file_name};
  if (!output_file.is_open()) {
    return false;
  }
  for (int outer_idx = 0; outer_idx < matrix.outerSize(); ++outer_idx) {
    for (SparseMatrixsc::InnerIterator it{matrix, outer_idx}; it; ++it) {
      // Matlab is 1 indexed
      output_file << it.row() + 1 << "\t" << it.col() + 1 << "\t" << it.value()
                  << std::endl;
    }
  }
  return true;
}
