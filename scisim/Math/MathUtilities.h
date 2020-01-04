// MathUtilities.h
//
// Breannan Smith
// Last updated: 09/28/2015

// TODO: Break these utilities into separate dense and sparse utility
// collections

#ifndef MATH_UTILITIES_H
#define MATH_UTILITIES_H

#include "scisim/Math/MathDefines.h"
#include <Eigen/Core>
#include <Eigen/LU>

#include "scisim/Utilities.h"

namespace MathUtilities {

// A 2D cross product
inline scalar cross(const Vector2s &a, const Vector2s &b) {
  return a.x() * b.y() - a.y() * b.x();
}

unsigned computeNumDigits(unsigned n);

// Checks if the the two vectors when stacked as [ a b ] form an orthonormal
// matrix with positive determinant
bool isRightHandedOrthoNormal(const Vector2s &a, const Vector2s &b,
                              const scalar &tol);
// Checks if the the three vectors when stacked as [ a b c ] form an orthonormal
// matrix with positive determinant
bool isRightHandedOrthoNormal(const Vector3s &a, const Vector3s &b,
                              const Vector3s &c, const scalar &tol);

// SPARSE MATRIX

bool isSquare(const SparseMatrixsc &matrix);

bool isIdentity(const SparseMatrixsc &A, const scalar &tol);

bool isSymmetric(const SparseMatrixsc &A, const scalar &tol);

// The number of non-zeros on or below the diagonal
int nzLowerTriangular(const SparseMatrixsc &A);

// Determine which elements on or below the diagonal are non-zero
int sparsityPatternLowerTriangular(const SparseMatrixsc &A, int *rows,
                                   int *cols);

// Extract elements
int values(const SparseMatrixsc &A, scalar *vals);

// Extract elements on or below the diagonal
int valuesLowerTriangular(const SparseMatrixsc &A, scalar *vals);

// Generates a Compressed Column Sparse array representation of a sparse matrix
void extractDataCCS(const SparseMatrixsc &A, VectorXi &col_ptr,
                    VectorXi &row_ind, VectorXs &val);

// Extracts columns in cols from A0, in order, and places them in A1
void extractColumns(const SparseMatrixsc &A0, const std::vector<unsigned> &cols,
                    SparseMatrixsc &A1);

void serialize(const SparseMatrixsc &A, std::ostream &stm);
void deserialize(SparseMatrixsc &A, std::istream &stm);

template <typename ScalarType>
inline void serialize(const Eigen::Quaternion<ScalarType> &a,
                      std::ostream &stm) {
  assert(stm.good());
  Utilities::serializeBuiltInType(a.x(), stm);
  Utilities::serializeBuiltInType(a.y(), stm);
  Utilities::serializeBuiltInType(a.z(), stm);
  Utilities::serializeBuiltInType(a.w(), stm);
}

template <typename ScalarType>
inline void deserialize(Eigen::Quaternion<ScalarType> &a, std::istream &stm) {
  assert(stm.good());
  a.x() = Utilities::deserialize<ScalarType>(stm);
  a.y() = Utilities::deserialize<ScalarType>(stm);
  a.z() = Utilities::deserialize<ScalarType>(stm);
  a.w() = Utilities::deserialize<ScalarType>(stm);
}

template <typename ScalarType>
inline Eigen::Quaternion<ScalarType> deserializeQuaternion(std::istream &stm) {
  assert(stm.good());
  Eigen::Quaternion<ScalarType> a;
  deserialize(a, stm);
  return a;
}

template <typename Derived>
void serialize(const Eigen::DenseBase<Derived> &eigen_variable,
               std::ostream &stm) {
  assert(stm.good());
  // If the size of either dimension is dynamic, it must be serialized
  if (Derived::RowsAtCompileTime == Eigen::Dynamic) {
    typename Derived::Index nrows = eigen_variable.rows();
    stm.write(reinterpret_cast<char *>(&nrows),
              sizeof(typename Eigen::DenseBase<Derived>::Index));
  }
  if (Derived::ColsAtCompileTime == Eigen::Dynamic) {
    typename Derived::Index ncols = eigen_variable.cols();
    stm.write(reinterpret_cast<char *>(&ncols),
              sizeof(typename Eigen::DenseBase<Derived>::Index));
  }
  // Write the data
  stm.write(const_cast<char *>(reinterpret_cast<const char *>(
                eigen_variable.derived().data())),
            eigen_variable.rows() * eigen_variable.cols() *
                sizeof(typename Derived::Scalar));
  assert(stm.good());
}

// Deserialization for dense Eigen types of fixed row count, fixed column count
template <typename Derived>
typename std::enable_if<Derived::RowsAtCompileTime != Eigen::Dynamic &&
                            Derived::ColsAtCompileTime != Eigen::Dynamic,
                        Derived>::type
deserialize(std::istream &stm) {
  assert(stm.good());
  Derived output_matrix;
  assert(output_matrix.rows() == Derived::RowsAtCompileTime);
  assert(output_matrix.cols() == Derived::ColsAtCompileTime);
  stm.read(reinterpret_cast<char *>(output_matrix.data()),
           Derived::RowsAtCompileTime * Derived::ColsAtCompileTime *
               sizeof(typename Derived::Scalar));
  assert(stm.good());
  return output_matrix;
}

// Deserialization for dense Eigen types of dynamic col count, fixed row count
template <typename Derived>
typename std::enable_if<Derived::RowsAtCompileTime != Eigen::Dynamic &&
                            Derived::ColsAtCompileTime == Eigen::Dynamic,
                        Derived>::type
deserialize(std::istream &stm) {
  assert(stm.good());
  Derived output_matrix;
  assert(output_matrix.rows() == Derived::RowsAtCompileTime);
  {
    typename Derived::Index ncols;
    stm.read(reinterpret_cast<char *>(&ncols), sizeof(typename Derived::Index));
    output_matrix.resize(Derived::RowsAtCompileTime, ncols);
  }
  stm.read(reinterpret_cast<char *>(output_matrix.data()),
           Derived::RowsAtCompileTime * output_matrix.cols() *
               sizeof(typename Derived::Scalar));
  assert(stm.good());
  return output_matrix;
}

// Deserialization for dense Eigen types of fixed col count, dynamic row count
template <typename Derived>
typename std::enable_if<Derived::RowsAtCompileTime == Eigen::Dynamic &&
                            Derived::ColsAtCompileTime != Eigen::Dynamic,
                        Derived>::type
deserialize(std::istream &stm) {
  assert(stm.good());
  Derived output_matrix;
  assert(output_matrix.cols() == Derived::ColsAtCompileTime);
  {
    typename Derived::Index nrows;
    stm.read(reinterpret_cast<char *>(&nrows), sizeof(typename Derived::Index));
    output_matrix.resize(nrows, Derived::ColsAtCompileTime);
  }
  stm.read(reinterpret_cast<char *>(output_matrix.data()),
           output_matrix.rows() * Derived::ColsAtCompileTime *
               sizeof(typename Derived::Scalar));
  assert(stm.good());
  return output_matrix;
}

template <typename T>
void serialize(const std::vector<T> &vector, std::ostream &stm) {
  assert(stm.good());
  // Write out the length of the vector
  Utilities::serializeBuiltInType(vector.size(), stm);
  // Output each element of the vector
  for (typename std::vector<T>::size_type idx = 0; idx < vector.size(); ++idx) {
    MathUtilities::serialize(vector[idx], stm);
  }
}

template <typename T>
void serialize(const std::vector<T, Eigen::aligned_allocator<T>> &vector,
               std::ostream &stm) {
  assert(stm.good());
  // Write out the length of the vector
  Utilities::serializeBuiltInType(vector.size(), stm);
  // Output each element of the vector
  for (typename std::vector<T>::size_type idx = 0; idx < vector.size(); ++idx) {
    MathUtilities::serialize(vector[idx], stm);
  }
}

// Deserialize a vector of Eigen types
template <typename T>
typename std::enable_if<Utilities::is_specialization<T, std::vector>::value,
                        T>::type
deserialize(std::istream &input_stream) {
  assert(input_stream.good());
  T vector;
  typename T::size_type length;
  Utilities::deserializeBuiltInType<typename T::size_type>(length,
                                                           input_stream);
  vector.reserve(length);
  for (typename T::size_type idx = 0; idx < length; ++idx) {
    vector.emplace_back(
        MathUtilities::deserialize<typename T::value_type>(input_stream));
    assert(input_stream.good());
  }
  assert(vector.size() == length);
  return vector;
}

// Saves a sparse matrix in a triplet format that can be imported into Matlab
// via:
//   load sparse_matrix.dat
//   H = spconvert( sparse_matrix )
bool writeToMatlabTripletText(const SparseMatrixsc &matrix,
                              const std::string &file_name);
bool writeToMatlabVectorText(const VectorXs &vector,
                             const std::string &file_name);

// TODO: Deserialization for dense Eigen types of dynamic row count, dynamic col
// count

} // namespace MathUtilities

#endif
