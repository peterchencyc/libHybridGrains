// MathUtilities.cpp
//
// Breannan Smith
// Last updated: 09/28/2015

#include "scisim/Math/MathUtilities.h"

#include <fstream>
#include <iostream>

bool MathUtilities::isRightHandedOrthoNormal(const Vector2s &a,
                                             const Vector2s &b,
                                             const scalar &tol) {
  // All basis vectors should be unit
  if (fabs(a.norm() - 1.0) > tol) {
    return false;
  }
  if (fabs(b.norm() - 1.0) > tol) {
    return false;
  }
  // All basis vectors should be mutually orthogonal
  if (fabs(a.dot(b)) > tol) {
    return false;
  }
  // Coordinate system should be right handed
  assert(fabs(cross(a, b) - 1.0) <= 1.0e-6 ||
         fabs(cross(a, b) + 1.0) <= 1.0e-6);
  if (cross(a, b) <= 0.0) {
    return false;
  }
  return true;
}

bool MathUtilities::isRightHandedOrthoNormal(const Vector3s &a,
                                             const Vector3s &b,
                                             const Vector3s &c,
                                             const scalar &tol) {
  // All basis vectors should be unit
  if (fabs(a.norm() - 1.0) > tol) {
    return false;
  }
  if (fabs(b.norm() - 1.0) > tol) {
    return false;
  }
  if (fabs(c.norm() - 1.0) > tol) {
    return false;
  }
  // All basis vectors should be mutually orthogonal
  if (fabs(a.dot(b)) > tol) {
    return false;
  }
  if (fabs(a.dot(c)) > tol) {
    return false;
  }
  if (fabs(b.dot(c)) > tol) {
    return false;
  }
  // Coordinate system should be right handed
  if ((a.cross(b) - c).lpNorm<Eigen::Infinity>() > tol) {
    return false;
  }
  return true;
}

bool MathUtilities::isSquare(const SparseMatrixsc &matrix) {
  return matrix.rows() == matrix.cols();
}

bool MathUtilities::isIdentity(const SparseMatrixsc &A, const scalar &tol) {
  if (!isSquare(A)) {
    return false;
  }
  for (int outer_idx = 0; outer_idx < A.outerSize(); ++outer_idx) {
    for (SparseMatrixsc::InnerIterator it(A, outer_idx); it; ++it) {
      if (it.row() == it.col()) {
        if (fabs(it.value() - 1.0) > tol) {
          return false;
        }
      } else {
        if (fabs(it.value()) > tol) {
          return false;
        }
      }
    }
  }
  return true;
}

bool MathUtilities::isSymmetric(const SparseMatrixsc &A, const scalar &tol) {
  const SparseMatrixsc &B{A - SparseMatrixsc{A.transpose()}};
  for (int outer_idx = 0; outer_idx < B.outerSize(); ++outer_idx) {
    for (SparseMatrixsc::InnerIterator it(B, outer_idx); it; ++it) {
      if (fabs(it.value()) > tol) {
        return false;
      }
    }
  }
  return true;
}

unsigned MathUtilities::computeNumDigits(unsigned n) {
  if (n == 0) {
    return 1;
  }
  unsigned num_digits{0};
  while (n != 0) {
    n /= 10;
    ++num_digits;
  }
  return num_digits;
}

int MathUtilities::nzLowerTriangular(const SparseMatrixsc &A) {
  int num{0};
  for (int col = 0; col < A.outerSize(); ++col) {
    for (SparseMatrixsc::InnerIterator it(A, col); it; ++it) {
      // Skip entries above the diagonal
      if (col > it.row()) {
        continue;
      }
      ++num;
    }
  }
  return num;
}

int MathUtilities::sparsityPatternLowerTriangular(const SparseMatrixsc &A,
                                                  int *rows, int *cols) {
  assert(rows != nullptr);
  assert(cols != nullptr);

  int curel{0};
  for (int col = 0; col < A.outerSize(); ++col) {
    for (SparseMatrixsc::InnerIterator it(A, col); it; ++it) {
      if (col > it.row())
        continue;
      rows[curel] = it.row();
      cols[curel] = col;
      ++curel;
    }
  }

  return curel;
}

int MathUtilities::values(const SparseMatrixsc &A, scalar *vals) {
  assert(vals != nullptr);

  int curel{0};
  for (int col = 0; col < A.outerSize(); ++col) {
    for (SparseMatrixsc::InnerIterator it(A, col); it; ++it) {
      vals[curel] = it.value();
      ++curel;
    }
  }

  assert(curel == A.nonZeros());
  return curel;
}

int MathUtilities::valuesLowerTriangular(const SparseMatrixsc &A,
                                         scalar *vals) {
  assert(vals != nullptr);

  int curel{0};
  for (int col = 0; col < A.outerSize(); ++col) {
    for (SparseMatrixsc::InnerIterator it(A, col); it; ++it) {
      if (col > it.row())
        continue;
      vals[curel] = it.value();
      ++curel;
    }
  }

  return curel;
}

void MathUtilities::extractDataCCS(const SparseMatrixsc &A, VectorXi &col_ptr,
                                   VectorXi &row_ind, VectorXs &val) {
  col_ptr.resize(A.cols() + 1);
  row_ind.resize(A.nonZeros());
  val.resize(A.nonZeros());

  col_ptr(0) = 0;
  for (int col = 0; col < A.outerSize(); ++col) {
    col_ptr(col + 1) = col_ptr(col);
    for (SparseMatrixsc::InnerIterator it(A, col); it; ++it) {
      const int row{int(it.row())};

      val(col_ptr(col + 1)) = it.value();
      row_ind(col_ptr(col + 1)) = row;
      ++col_ptr(col + 1);
    }
  }

  assert(col_ptr(col_ptr.size() - 1) == row_ind.size());
}

// TODO: Pull the outerIndexPtr arithmetic into a helper function
void MathUtilities::extractColumns(const SparseMatrixsc &A0,
                                   const std::vector<unsigned> &cols,
                                   SparseMatrixsc &A1) {
  const unsigned ncols_to_extract{static_cast<unsigned>(cols.size())};

  assert(ncols_to_extract <= static_cast<unsigned>(A0.cols()));
#ifndef NDEBUG
  for (unsigned i = 0; i < ncols_to_extract; ++i) {
    assert(cols[i] < unsigned(A0.cols()));
  }
#endif

  // Compute the number of nonzeros in each column of the new matrix
  VectorXi column_nonzeros{ncols_to_extract};
  for (unsigned i = 0; i < ncols_to_extract; ++i) {
    column_nonzeros(i) =
        A0.outerIndexPtr()[cols[i] + 1] - A0.outerIndexPtr()[cols[i]];
  }

  // Resize A1 and reserve space
  A1.resize(A0.rows(), ncols_to_extract);
  A1.reserve(column_nonzeros);
  // Copy the data over, column by column
  for (unsigned cur_col = 0; cur_col < ncols_to_extract; ++cur_col) {
    for (SparseMatrixsc::InnerIterator it(A0, cols[cur_col]); it; ++it) {
      A1.insert(it.row(), cur_col) = it.value();
    }
  }

  A1.makeCompressed();

#ifndef NDEBUG
  for (int i = 0; i < A1.cols(); ++i) {
    assert((A1.outerIndexPtr()[i + 1] - A1.outerIndexPtr()[i]) ==
           column_nonzeros(i));
  }
#endif
}

void MathUtilities::serialize(const SparseMatrixsc &A, std::ostream &stm) {
  assert(stm.good());
  assert(A.isCompressed());

  Utilities::serializeBuiltInType(A.rows(), stm);
  Utilities::serializeBuiltInType(A.cols(), stm);
  Utilities::serializeBuiltInType(A.nonZeros(), stm);
  // TODO: Write a utility class to save out pointers to arrays of data
  stm.write(reinterpret_cast<const char *>(A.innerIndexPtr()),
            A.nonZeros() * sizeof(SparseMatrixsc::StorageIndex));
  stm.write(reinterpret_cast<const char *>(A.outerIndexPtr()),
            (A.outerSize() + 1) * sizeof(SparseMatrixsc::StorageIndex));
  stm.write(reinterpret_cast<const char *>(A.valuePtr()),
            A.nonZeros() * sizeof(SparseMatrixsc::Scalar));
}

void MathUtilities::deserialize(SparseMatrixsc &A, std::istream &stm) {
  assert(stm.good());

  const Eigen::Index rows{Utilities::deserialize<Eigen::Index>(stm)};
  const Eigen::Index cols{Utilities::deserialize<Eigen::Index>(stm)};
  const Eigen::Index nnz{Utilities::deserialize<Eigen::Index>(stm)};

  Eigen::Matrix<SparseMatrixsc::StorageIndex, Eigen::Dynamic, 1> inner_indices{
      nnz};
  stm.read(reinterpret_cast<char *>(inner_indices.data()),
           inner_indices.size() * sizeof(SparseMatrixsc::StorageIndex));

  Eigen::Matrix<SparseMatrixsc::StorageIndex, Eigen::Dynamic, 1> outter_indices{
      cols + 1};
  stm.read(reinterpret_cast<char *>(outter_indices.data()),
           outter_indices.size() * sizeof(SparseMatrixsc::StorageIndex));

  Eigen::Matrix<SparseMatrixsc::Scalar, Eigen::Dynamic, 1> values{nnz};
  stm.read(reinterpret_cast<char *>(values.data()),
           values.size() * sizeof(SparseMatrixsc::Scalar));

  const Eigen::Map<const SparseMatrixsc> matrix_map{
      rows,         cols, nnz, outter_indices.data(), inner_indices.data(),
      values.data()};

  A = matrix_map;
  assert(A.isCompressed());
}

bool MathUtilities::writeToMatlabTripletText(const SparseMatrixsc &matrix,
                                             const std::string &file_name) {
  std::ofstream output_file(file_name);
  if (!output_file.is_open()) {
    return false;
  }
  for (int outer_idx = 0; outer_idx < matrix.outerSize(); ++outer_idx) {
    for (SparseMatrixsc::InnerIterator it(matrix, outer_idx); it; ++it) {
      // Matlab is 1 indexed
      output_file << it.row() + 1 << "\t" << it.col() + 1 << "\t" << it.value()
                  << std::endl;
    }
  }
  return true;
}

bool MathUtilities::writeToMatlabVectorText(const VectorXs &vector,
                                            const std::string &file_name) {
  std::ofstream output_file(file_name);
  if (!output_file.is_open()) {
    return false;
  }
  output_file << vector << std::endl;
  return true;
}
