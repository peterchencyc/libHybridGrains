// HDF5File.h
//
// Breannan Smith
// Last updated: 09/28/2015

// TODO: Store sparse matrix components in a struct to prevent polution of
// namespace
// TODO: Support routines for users to create structs
// TODO: Don't make explicit copies if conversion needed, just use hyperslabs
// TODO: Throw error if file not open when attempting to write
// TODO: Check close return types

#ifndef HDF5_FILE_H
#define HDF5_FILE_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <string>

#ifdef USE_HDF5
#include "hdf5.h"

// Wrapper for HDF5 group id. Functionality for HDF5 types differs only in the
// required close operaiton.
template <herr_t H5CloseOperation(hid_t id)> class HDFID final {

public:
  HDFID() : m_hid_t(-1) {}

  explicit HDFID(const hid_t value) : m_hid_t(value) {}

  ~HDFID() {
    if (m_hid_t >= 0) {
      H5CloseOperation(m_hid_t);
    }
  }

  HDFID(HDFID &&other) : m_hid_t(other.m_hid_t) { other.m_hid_t = -1; }

  HDFID &operator=(HDFID &&other) {
    // NB: Backup hid_t to guard against self assignment
    const hid_t others_hid_t{other.m_hid_t};
    other.m_hid_t = -1;
    m_hid_t = others_hid_t;
    return *this;
  }

  operator hid_t() const { return m_hid_t; }

private:
  HDFID(const HDFID &) = delete;
  HDFID &operator=(const HDFID &) = delete;

  hid_t m_hid_t;
};
#endif

namespace HDF5SupportedTypes {
template <typename T> constexpr bool isSupportedEigenType() { return false; }

template <> constexpr bool isSupportedEigenType<int>() { return true; }

template <> constexpr bool isSupportedEigenType<long>() { return true; }

template <> constexpr bool isSupportedEigenType<unsigned>() { return true; }

template <> constexpr bool isSupportedEigenType<unsigned long>() {
  return true;
}

template <> constexpr bool isSupportedEigenType<float>() { return true; }

template <> constexpr bool isSupportedEigenType<double>() { return true; }
} // namespace HDF5SupportedTypes

enum class HDF5AccessType : std::uint8_t { READ_ONLY, READ_WRITE };

class HDF5File final {

public:
  HDF5File();
#ifndef USE_HDF5
  [[noreturn]]
#endif
  HDF5File(const std::string &file_name, const HDF5AccessType &access_type);
  ~HDF5File();

#ifdef USE_HDF5
  hid_t fileID();
#endif

#ifndef USE_HDF5
  [[noreturn]]
#endif
  void
  open(const std::string &file_name, const HDF5AccessType &access_type);

  bool is_open() const;

#ifdef USE_HDF5
  HDFID<H5Gclose> getGroup(const std::string &group_name) const;
#endif

#ifdef USE_HDF5
  HDFID<H5Gclose> findGroup(const std::string &group_name) const;
#endif

#ifndef USE_HDF5
  [[noreturn]]
#endif
  void
  writeString(const std::string &group, const std::string &variable_name,
              const std::string &string_variable) const;

  void readString(const std::string &group, const std::string &variable_name,
                  std::string &string_variable) const;

  template <typename Scalar>
  void writeScalar(const std::string &group, const std::string &variable_name,
                   const Scalar &variable) const {
    Eigen::Matrix<Scalar, 1, 1> output_mat;
    output_mat << variable;
    writeMatrix(group, variable_name, output_mat);
  }

  template <typename Scalar>
  void readScalar(const std::string &group, const std::string &variable_name,
                  Scalar &variable) const {
    Eigen::Matrix<Scalar, 1, 1> input_mat;
    readMatrix(group, variable_name, input_mat);
    variable = input_mat(0, 0);
  }

  template <typename Derived>
  void writeMatrix(const std::string &group, const std::string &variable_name,
                   const Eigen::DenseBase<Derived> &eigen_variable) const {
#ifdef USE_HDF5
    using HDFSID = HDFID<H5Sclose>;
    using HDFGID = HDFID<H5Gclose>;
    using HDFDID = HDFID<H5Dclose>;

    using Scalar = typename Derived::Scalar;
    static_assert(HDF5SupportedTypes::isSupportedEigenType<Scalar>(),
                  "Error, scalar type of Eigen variable must be float, double, "
                  "unsigned or integer");

    assert(eigen_variable.rows() >= 0);
    assert(eigen_variable.cols() >= 0);
    const hsize_t dims[2] = {hsize_t(eigen_variable.rows()),
                             hsize_t(eigen_variable.cols())};
    const HDFSID dataspace_id{H5Screate_simple(2, dims, nullptr)};
    if (dataspace_id < 0) {
      throw std::string{"Failed to create HDF data space"};
    }

    // Open the requested group
    const HDFGID grp_id{getGroup(group)};

    const HDFDID dataset_id{H5Dcreate2(
        grp_id, variable_name.c_str(), computeHDFType(eigen_variable),
        dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
    if (dataset_id < 0) {
      throw std::string{"Failed to create HDF data set"};
    }

    // Convert to row major format, if needed
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        col_major_output_data;
    if (isColumnMajor(eigen_variable)) {
      col_major_output_data.resize(eigen_variable.rows(),
                                   eigen_variable.cols());
      col_major_output_data = eigen_variable.derived().matrix();
    }

    const herr_t status_write{H5Dwrite(
        dataset_id, computeHDFType(eigen_variable), H5S_ALL, H5S_ALL,
        H5P_DEFAULT,
        isColumnMajor(eigen_variable) ? col_major_output_data.data()
                                      : eigen_variable.derived().data())};
    if (status_write < 0) {
      throw std::string{"Failed to write HDF data"};
    }
#else
    throw std::string{"HDF5File::writeMatrix not compiled with HDF5 support"};
#endif
  }

  template <typename Derived>
  void readMatrix(const std::string &group, const std::string &variable_name,
                  Eigen::DenseBase<Derived> &eigen_variable) const {
#ifdef USE_HDF5
    using HDFDID = HDFID<H5Dclose>;
    using HDFGID = HDFID<H5Gclose>;

    using Scalar = typename Derived::Scalar;
    static_assert(HDF5SupportedTypes::isSupportedEigenType<Scalar>(),
                  "Error, scalar type of Eigen variable must be float, double, "
                  "unsigned or integer");

    // Open the requested group
    const HDFGID grp_id{findGroup(group)};

    const HDFDID dataset_id{
        H5Dopen2(grp_id, variable_name.c_str(), H5P_DEFAULT)};
    if (dataset_id < 0) {
      throw std::string{"Failed to open HDF data set"};
    }
    if (getNativeType(dataset_id) != computeHDFType(eigen_variable)) {
      throw std::string{
          "Requested HDF data set is not of given type from Eigen variable"};
    }

    Eigen::ArrayXi dimensions;
    getDimensions(dataset_id, dimensions);
    if (dimensions.size() != 2) {
      throw std::string{"Invalid dimensions for Eigen matrix type in file"};
    }
    if ((dimensions < 0).any()) {
      throw std::string{"Negative dimensions for Eigen matrix type in file"};
    }
    // If the Eigen type has a fixed dimension, ensure it is correct
    if (rowsFixed(eigen_variable) && eigen_variable.rows() != dimensions(0)) {
      throw std::string{
          "Eigen type of fixed row size does not have correct number of rows"};
    }
    if (colsFixed(eigen_variable) && eigen_variable.cols() != dimensions(1)) {
      throw std::string{
          "Eigen type of fixed cols size does not have correct number of cols"};
    }
    // Resize the Eigen type
    eigen_variable.derived().resize(dimensions(0), dimensions(1));

    // Convert to row major format, if needed
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        col_major_input_data;
    if (isColumnMajor(eigen_variable)) {
      col_major_input_data.resize(eigen_variable.rows(), eigen_variable.cols());
    }

    const herr_t read_status{H5Dread(
        dataset_id, getNativeType(dataset_id), H5S_ALL, H5S_ALL, H5P_DEFAULT,
        isColumnMajor(eigen_variable) ? col_major_input_data.data()
                                      : eigen_variable.derived().data())};
    if (read_status < 0) {
      throw std::string{"Failed to read data from HDF file"};
    }
    if (isColumnMajor(eigen_variable)) {
      eigen_variable.derived().matrix() = col_major_input_data;
    }
#else
    throw std::string{"HDF5File::readMatrix not compiled with HDF5 support"};
#endif
  }

  template <typename Derived>
  void writeSparseMatrix(
      const std::string &group, const std::string &variable_name,
      const Eigen::SparseMatrixBase<Derived> &sparse_matrix) const {
    static_assert(!Derived::IsRowMajor, "Error, HDF5 sparse matrix output only "
                                        "supported for column major matrices");

    // Grab the sparse data
    Eigen::Matrix<typename Derived::StorageIndex, Eigen::Dynamic, 1> col_ptr;
    Eigen::Matrix<typename Derived::StorageIndex, Eigen::Dynamic, 1> row_ind;
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> val;
    extractDataCCS(sparse_matrix, col_ptr, row_ind, val);

    // Write out the sparse data
    const Eigen::Index nrows{sparse_matrix.rows()};
    writeScalar(group, variable_name + "_nrows", nrows);
    const Eigen::Index ncols{sparse_matrix.cols()};
    writeScalar(group, variable_name + "_ncols", ncols);
    writeMatrix(group, variable_name + "_col_ptr", col_ptr);
    writeMatrix(group, variable_name + "_row_ind", row_ind);
    writeMatrix(group, variable_name + "_val", val);
  }

  template <typename Derived>
  void readSparseMatrix(const std::string &group,
                        const std::string &variable_name,
                        Eigen::SparseMatrixBase<Derived> &sparse_matrix) const {
    static_assert(!Derived::IsRowMajor, "Error, HDF5 sparse matrix input only "
                                        "supported for column major matrices");

    // Read in the sparse data
    std::vector<Eigen::Triplet<typename Derived::Scalar>> triplet_list;
    {
      typename Derived::Index nrows;
      readScalar(group, variable_name + "_nrows", nrows);
      typename Derived::Index ncols;
      readScalar(group, variable_name + "_ncols", ncols);
      Eigen::Matrix<typename Derived::Index, Eigen::Dynamic, 1> col_ptr;
      readMatrix(group, variable_name + "_col_ptr", col_ptr);
      if (col_ptr.size() != ncols + 1) {
        throw std::string{
            "Column count and size of col_ptr do not agree for sparse matrix"};
      }
      Eigen::Matrix<typename Derived::Index, Eigen::Dynamic, 1> row_ind;
      readMatrix(group, variable_name + "_row_ind", row_ind);
      Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> val;
      readMatrix(group, variable_name + "_val", val);
      if (row_ind.size() != val.size()) {
        throw std::string{
            "Row indices and value size do not agree for sparse matrix"};
      }
      sparse_matrix.derived().resize(nrows, ncols);
      sparse_matrix.derived().reserve(int(val.size()));
      for (typename Derived::Index col_num = 0; col_num < ncols; ++col_num) {
        assert(col_ptr(col_num) <= col_ptr(col_num + 1));
        for (typename Derived::Index idx = col_ptr(col_num);
             idx < col_ptr(col_num + 1); ++idx) {
          typename Derived::Index row_num = row_ind(idx);
          triplet_list.push_back(Eigen::Triplet<typename Derived::Scalar>(
              row_num, col_num, val(idx)));
        }
      }
    }
    sparse_matrix.derived().setFromTriplets(triplet_list.cbegin(),
                                            triplet_list.cend());
    sparse_matrix.derived().makeCompressed();
  }

  // Doesn't work because HOFFSET is a macro... find a workaround
  // template < typename StructType, typename ScalarType, unsigned N >
  // void createStructType( const hid_t struct_tid, const ScalarType(&
  // struct_entry)[N], const std::string& struct_entry_name )
  //{
  //  #ifdef USE_HDF5
  //  using HDFTID = HDFID<H5Tclose>;
  //
  //  const hsize_t array_dim[]{ N };
  //  const HDFTID array_tid = H5Tarray_create2( computeHDFType<ScalarType>(),
  //  1, array_dim ); if( array_tid < 0 )
  //  {
  //    throw std::string{ "Failed to create HDF " + struct_entry_name + " type
  //    for static planes" };
  //  }
  //  if( H5Tinsert( struct_tid, struct_entry_name.c_str(),
  //  HOFFSET(StructType,struct_entry), array_tid ) < 0 )
  //  {
  //    throw std::string{ "Failed to insert " + struct_entry_name + " in HDF
  //    struct for static planes" };
  //  }
  //  #else
  //  throw std::string{ "HDF5File::createStructType not compiled with HDF5
  //  support" }; #endif
  //}

private:
  // TODO: Don't extract data first, just directly save out to disk
  template <typename Derived>
  static void extractDataCCS(
      const Eigen::SparseMatrixBase<Derived> &A,
      Eigen::Matrix<typename Derived::StorageIndex, Eigen::Dynamic, 1>
          &outer_ptr,
      Eigen::Matrix<typename Derived::StorageIndex, Eigen::Dynamic, 1>
          &inner_ptr,
      Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> &val) {
    outer_ptr = Eigen::Map<
        const Eigen::Matrix<typename Derived::StorageIndex, Eigen::Dynamic, 1>>(
        A.derived().outerIndexPtr(), A.outerSize() + 1);
    inner_ptr = Eigen::Map<
        const Eigen::Matrix<typename Derived::StorageIndex, Eigen::Dynamic, 1>>(
        A.derived().innerIndexPtr(), A.derived().nonZeros());
    val = Eigen::Map<
        const Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>>(
        A.derived().valuePtr(), A.derived().nonZeros());
  }

#ifdef USE_HDF5
  static hid_t getNativeType(const hid_t dataset_id) {
    using HDFTID = HDFID<H5Tclose>;

    const HDFTID dtype_id{H5Dget_type(dataset_id)};
    if (dtype_id < 0) {
      throw std::string{"Failed to get HDF5 datatype from dataset"};
    }
    if (H5Tequal(dtype_id, H5T_NATIVE_DOUBLE) > 0) {
      return H5T_NATIVE_DOUBLE;
    } else if (H5Tequal(dtype_id, H5T_NATIVE_FLOAT) > 0) {
      return H5T_NATIVE_FLOAT;
    } else if (H5Tequal(dtype_id, H5T_NATIVE_INT) > 0) {
      return H5T_NATIVE_INT;
    } else if (H5Tequal(dtype_id, H5T_NATIVE_UINT) > 0) {
      return H5T_NATIVE_UINT;
    } else
      return -1;
  }

  static void getDimensions(const hid_t dataset_id,
                            Eigen::ArrayXi &dimensions) {
    using HDFSID = HDFID<H5Sclose>;

    const HDFSID space_id{H5Dget_space(dataset_id)};
    if (space_id < 0) {
      throw std::string{"Failed to open data space"};
    }
    const int rank{H5Sget_simple_extent_ndims(space_id)};
    if (rank < 0) {
      throw std::string{"Failed to get rank"};
    }
    dimensions.resize(rank);
    std::vector<hsize_t> dims(
        static_cast<std::vector<hsize_t>::size_type>(rank));
    assert(int(dims.size()) == rank);
    const herr_t status_get_simple_extent_dims{
        H5Sget_simple_extent_dims(space_id, dims.data(), nullptr)};
    if (status_get_simple_extent_dims < 0) {
      throw std::string{"Failed to get extents"};
    }
    for (int i = 0; i < rank; ++i) {
      dimensions(i) = int(dims[i]);
    }
  }

  template <typename ScalarType> static constexpr hid_t computeHDFType() {
    using std::is_same;
    return is_same<ScalarType, double>::value
               ? H5T_NATIVE_DOUBLE
               : is_same<ScalarType, float>::value
                     ? H5T_NATIVE_FLOAT
                     : is_same<ScalarType, int>::value
                           ? H5T_NATIVE_INT
                           : is_same<ScalarType, unsigned>::value
                                 ? H5T_NATIVE_UINT
                                 : is_same<ScalarType, unsigned long>::value
                                       ? H5T_NATIVE_ULONG
                                       : is_same<ScalarType, long>::value
                                             ? H5T_NATIVE_LONG
                                             : -1;
  }

  template <typename Derived>
  static constexpr hid_t computeHDFType(const Eigen::EigenBase<Derived> &) {
    return computeHDFType<typename Derived::Scalar>();
  }

  template <typename Derived>
  static constexpr bool isColumnMajor(const Eigen::EigenBase<Derived> &) {
    return !Derived::IsRowMajor;
  }

  template <typename Derived>
  static constexpr bool rowsFixed(const Eigen::EigenBase<Derived> &) {
    return Derived::RowsAtCompileTime != Eigen::Dynamic;
  }

  template <typename Derived>
  static constexpr bool colsFixed(const Eigen::EigenBase<Derived> &) {
    return Derived::ColsAtCompileTime != Eigen::Dynamic;
  }
#endif

#ifdef USE_HDF5
  hid_t m_hdf_file_id;
#endif
  bool m_file_opened;
};

#endif
