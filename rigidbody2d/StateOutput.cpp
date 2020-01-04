// StateOutput.cpp
//
// Breannan Smith
// Last updated: 01/08/2016

#include "StateOutput.h"

#include "BoxGeometry.h"
#include "CircleGeometry.h"
#include "PlanarPortal.h"
#include "RigidBody2DGeometry.h"
#include "RigidBody2DStaticDrum.h"
#include "RigidBody2DStaticPlane.h"
#include "scisim/HDF5File.h"
#include "scisim/Utilities.h"

#include <iostream>

void RigidBody2DStateOutput::writeGeometryIndices(
    const std::vector<std::unique_ptr<RigidBody2DGeometry>> &geometry,
    const VectorXu &indices, const std::string &group, HDF5File &output_file) {
#ifdef USE_HDF5
  using enum_type =
      std::uint8_t; // std::underlying_type<RigidBody2DGeometryType>::type;
  using HDFSID = HDFID<H5Sclose>;
  using HDFGID = HDFID<H5Gclose>;
  using HDFDID = HDFID<H5Dclose>;

  // Map global indices to 'local' indices; that is, given a global index, gives
  // the index into the particular type (e.g. global index 10 could be sphere
  // number 3)
  VectorXu global_local_geo_mapping(geometry.size());
  {
    unsigned current_geo_idx{0};
    Vector3u geo_type_local_indices{Vector3u::Zero()};
    for (const std::unique_ptr<RigidBody2DGeometry> &current_geo : geometry) {
      assert(enum_type(current_geo->type()) <=
             geo_type_local_indices.size() - 1);
      global_local_geo_mapping(current_geo_idx++) =
          geo_type_local_indices(enum_type(current_geo->type()))++;
    }
  }

  // Create an HDF5 dataspace
  const hsize_t dim[]{2, hsize_t(indices.size())};
  const HDFSID data_space{H5Screate_simple(2, dim, nullptr)};
  if (data_space < 0) {
    throw std::string{"Failed to create HDF dataspace for geometry indices"};
  }

  const HDFGID grp_id{output_file.getGroup(group)};

  // Create an HDF5 dataset
  const HDFDID data_set{H5Dcreate2(grp_id, "geometry_indices", H5T_NATIVE_UINT,
                                   data_space, H5P_DEFAULT, H5P_DEFAULT,
                                   H5P_DEFAULT)};
  if (data_set < 0) {
    throw std::string{"Failed to create HDF dataset for geometry indices"};
  }

  // Create an HDF5 memspace to allow us to insert elements one by one
  const HDFSID mem_space{H5Screate_simple(1, dim, nullptr)};
  if (mem_space < 0) {
    throw std::string{"Failed to create HDF memspace for geometry indices"};
  }

  // Insert the geometry indices one by one
  unsigned current_geo{0};
  for (unsigned idx = 0; idx < indices.size(); ++idx) {
    const unsigned geo_idx{indices(idx)};
    const hsize_t count[]{2, 1};
    const hsize_t offset[]{0, current_geo++};
    const hsize_t mem_offset[]{0, 0};
    H5Sselect_hyperslab(data_space, H5S_SELECT_SET, offset, nullptr, count,
                        nullptr);
    H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, mem_offset, nullptr, count,
                        nullptr);
    const unsigned data[]{unsigned(geometry[geo_idx]->type()),
                          global_local_geo_mapping(geo_idx)};
    if (H5Dwrite(data_set, H5T_NATIVE_UINT, mem_space, data_space, H5P_DEFAULT,
                 &data) < 0) {
      throw std::string{"Failed to write geometry indices to HDF"};
    }
  }
  assert(current_geo == indices.size());
#else
  throw std::string{"writeGeometryIndices not compiled with HDF5 support"};
#endif
}

// TODO: Combine the write circle geometry and write box geometry code into one
// function
#ifndef USE_HDF5
[[noreturn]]
#endif
static void
writeCircleGeometry(
    const std::vector<std::unique_ptr<RigidBody2DGeometry>> &geometry,
    const unsigned circle_count, const std::string &group,
    HDF5File &output_file) {
#ifdef USE_HDF5
  using HDFTID = HDFID<H5Tclose>;
  using HDFSID = HDFID<H5Sclose>;
  using HDFGID = HDFID<H5Gclose>;
  using HDFDID = HDFID<H5Dclose>;

  struct CircleData {
    scalar r;
  };

  // Create an HDF5 dataspace
  const hsize_t dim[]{circle_count};
  const HDFSID data_space{H5Screate_simple(1, dim, nullptr)};
  if (data_space < 0) {
    throw std::string{"Failed to create HDF dataspace for circle geometry"};
  }

  // Create an HDF5 struct for the data
  const HDFTID struct_tid{H5Tcreate(H5T_COMPOUND, sizeof(CircleData))};
  if (struct_tid < 0) {
    throw std::string{"Failed to create HDF struct"};
  }
  // Insert the r type in the struct
  if (H5Tinsert(struct_tid, "r", HOFFSET(CircleData, r), H5T_NATIVE_DOUBLE) <
      0) {
    throw std::string{"Failed to create HDF r type for circle geometry"};
  }

  const HDFGID grp_id{output_file.getGroup(group)};

  // Create an HDF5 dataset
  const HDFDID data_set{H5Dcreate2(grp_id, "circles", struct_tid, data_space,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
  if (data_set < 0) {
    throw std::string{"Failed to create HDF dataset for circle geometry"};
  }

  // Create an HDF5 memspace to allow us to insert elements one by one
  const HDFSID mem_space{H5Screate_simple(1, dim, nullptr)};
  if (mem_space < 0) {
    throw std::string{"Failed to create HDF memspace for circle geometry"};
  }

  // Insert the spheres one by one
  unsigned current_circle{0};
  CircleData data;
  for (const std::unique_ptr<RigidBody2DGeometry> &geometry_instance :
       geometry) {
    switch (geometry_instance->type()) {
    case RigidBody2DGeometryType::CIRCLE: {
      const CircleGeometry &circle{
          static_cast<const CircleGeometry &>(*geometry_instance)};
      hsize_t count[]{1};
      hsize_t offset[]{current_circle++};
      hsize_t mem_offset[]{0};
      H5Sselect_hyperslab(data_space, H5S_SELECT_SET, offset, nullptr, count,
                          nullptr);
      H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, mem_offset, nullptr, count,
                          nullptr);
      data.r = circle.r();
      if (H5Dwrite(data_set, struct_tid, mem_space, data_space, H5P_DEFAULT,
                   &data) < 0) {
        throw std::string{"Failed to write circle geometry struct to HDF"};
      }
      break;
    }
    case RigidBody2DGeometryType::BOX: {
      // Nothing to do for boxes here
      break;
    }
    case RigidBody2DGeometryType::ANNULUS: {
      // Nothing to do for annuli here
      break;
    }
    }
  }
  assert(current_circle == circle_count);
#else
  throw std::string{"writeCircleGeometry not compiled with HDF5 support"};
#endif
}

#ifndef USE_HDF5
[[noreturn]]
#endif
static void
writeBoxGeometry(
    const std::vector<std::unique_ptr<RigidBody2DGeometry>> &geometry,
    const unsigned box_count, const std::string &group, HDF5File &output_file) {
#ifdef USE_HDF5
  using HDFTID = HDFID<H5Tclose>;
  using HDFSID = HDFID<H5Sclose>;
  using HDFGID = HDFID<H5Gclose>;
  using HDFDID = HDFID<H5Dclose>;

  struct BoxData {
    scalar r[2];
  };

  // Create an HDF5 dataspace
  const hsize_t dim[]{box_count};
  const HDFSID data_space{H5Screate_simple(1, dim, nullptr)};
  if (data_space < 0) {
    throw std::string{"Failed to create HDF dataspace for box geometry"};
  }

  // Create an HDF5 struct for the data
  const HDFTID struct_tid{H5Tcreate(H5T_COMPOUND, sizeof(BoxData))};
  if (struct_tid < 0) {
    throw std::string{"Failed to create HDF struct"};
  }
  // Insert the r type in the struct
  {
    const hsize_t array_dim[]{2};
    const HDFTID array_tid{H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, array_dim)};
    if (array_tid < 0) {
      throw std::string{"Failed to create HDF r type for box geometry"};
    }
    if (H5Tinsert(struct_tid, "r", HOFFSET(BoxData, r), array_tid) < 0) {
      throw std::string{"Failed to insert r in HDF struct for box geometry"};
    }
  }

  const HDFGID grp_id{output_file.getGroup(group)};

  // Create an HDF5 dataset
  const HDFDID data_set{H5Dcreate2(grp_id, "boxes", struct_tid, data_space,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)};
  if (data_set < 0) {
    throw std::string{"Failed to create HDF dataset for box geometry"};
  }

  // Create an HDF5 memspace to allow us to insert elements one by one
  const HDFSID mem_space{H5Screate_simple(1, dim, nullptr)};
  if (mem_space < 0) {
    throw std::string{"Failed to create HDF memspace for box geometry"};
  }

  // Insert the boxes one by one
  unsigned current_box{0};
  BoxData data;
  for (const std::unique_ptr<RigidBody2DGeometry> &geometry_instance :
       geometry) {
    switch (geometry_instance->type()) {
    case RigidBody2DGeometryType::CIRCLE: {
      // Nothing to do for circles here
      break;
    }
    case RigidBody2DGeometryType::BOX: {
      const BoxGeometry &box{
          static_cast<const BoxGeometry &>(*geometry_instance)};
      Eigen::Map<Vector2s>{data.r} = box.r();
      const hsize_t count[]{1};
      const hsize_t offset[]{current_box++};
      const hsize_t mem_offset[]{0};
      H5Sselect_hyperslab(data_space, H5S_SELECT_SET, offset, nullptr, count,
                          nullptr);
      H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, mem_offset, nullptr, count,
                          nullptr);
      if (H5Dwrite(data_set, struct_tid, mem_space, data_space, H5P_DEFAULT,
                   &data) < 0) {
        throw std::string{"Failed to write box geometry struct to HDF"};
      }
      break;
    }
    case RigidBody2DGeometryType::ANNULUS: {
      // Nothing to do for annuli here
      break;
    }
    }
  }
  assert(current_box == box_count);
#else
  throw std::string{"writeBoxGeometry not compiled with HDF5 support"};
#endif
}

void RigidBody2DStateOutput::writeGeometry(
    const std::vector<std::unique_ptr<RigidBody2DGeometry>> &geometry,
    const std::string &group, HDF5File &output_file) {
  using enum_type =
      std::uint8_t; // std::underlying_type<RigidBody2DGeometryType>::type;

  Vector3u body_count{Vector3u::Zero()};
  for (const std::unique_ptr<RigidBody2DGeometry> &geometry_instance :
       geometry) {
    assert(enum_type(geometry_instance->type()) < body_count.size());
    ++body_count(enum_type(geometry_instance->type()));
  }
  assert(body_count.sum() == geometry.size());

  if (body_count(0) != 0) {
    std::cerr << "Error, annulus output not yet coded up for "
                 "RigidBody2DStateOutput::writeGeometry!"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (body_count(1) != 0) {
    writeCircleGeometry(geometry, body_count(1), group, output_file);
  }
  if (body_count(2) != 0) {
    writeBoxGeometry(geometry, body_count(2), group, output_file);
  }
}

// TODO: Abstract out the struct type creation into own routines
void RigidBody2DStateOutput::writeStaticPlanes(
    const std::vector<RigidBody2DStaticPlane> &static_planes,
    const std::string &group, HDF5File &output_file) {
#ifdef USE_HDF5
  using HDFTID = HDFID<H5Tclose>;
  using HDFSID = HDFID<H5Sclose>;
  using HDFGID = HDFID<H5Gclose>;
  using HDFDID = HDFID<H5Dclose>;

  struct LocalStaticPlaneData {
    scalar x[2];
    scalar n[2];
    // scalar t[2]; // Tangent can be recovered from n by user
  };

  // Create an HDF5 dataspace
  const hsize_t dim[]{static_planes.size()};
  const HDFSID data_space{H5Screate_simple(1, dim, nullptr)};
  if (data_space < 0) {
    throw std::string{"Failed to create HDF dataspace for static planes"};
  }

  // Create an HDF5 struct for the data
  const HDFTID struct_tid{
      H5Tcreate(H5T_COMPOUND, sizeof(LocalStaticPlaneData))};
  if (struct_tid < 0) {
    throw std::string{"Failed to create HDF struct for static planes"};
  }
  // Insert the x type in the struct
  {
    const hsize_t array_dim[]{2};
    const HDFTID array_tid{H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, array_dim)};
    if (array_tid < 0) {
      throw std::string{"Failed to create HDF x type for static planes"};
    }
    if (H5Tinsert(struct_tid, "x", HOFFSET(LocalStaticPlaneData, x),
                  array_tid) < 0) {
      throw std::string{"Failed to insert x in HDF struct for static planes"};
    }
  }
  // Insert the R type in the struct
  {
    const hsize_t array_dim[]{2};
    const HDFTID array_tid{H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, array_dim)};
    if (array_tid < 0) {
      throw std::string{"Failed to create HDF n type for static planes"};
    }
    if (H5Tinsert(struct_tid, "n", HOFFSET(LocalStaticPlaneData, n),
                  array_tid) < 0) {
      throw std::string{"Failed to insert n in HDF struct for static planes"};
    }
  }

  const HDFGID grp_id{output_file.getGroup(group)};

  // Create an HDF5 dataset
  const HDFDID data_set{H5Dcreate2(grp_id, "static_planes", struct_tid,
                                   data_space, H5P_DEFAULT, H5P_DEFAULT,
                                   H5P_DEFAULT)};
  if (data_set < 0) {
    throw std::string{"Failed to create HDF dataset for static planes"};
  }

  // Create an HDF5 memspace to allow us to insert elements one by one
  const HDFSID mem_space{H5Screate_simple(1, dim, nullptr)};
  if (mem_space < 0) {
    throw std::string{"Failed to create HDF memspace for static planes"};
  }

  // Insert the static planes one by one
  unsigned current_plane{0};
  LocalStaticPlaneData local_data;
  for (const RigidBody2DStaticPlane &plane : static_planes) {
    Eigen::Map<Vector2s>{local_data.x} = plane.x();
    Eigen::Map<Vector2s>{local_data.n} = plane.n();
    const hsize_t count[]{1};
    const hsize_t offset[]{current_plane++};
    const hsize_t mem_offset[]{0};
    H5Sselect_hyperslab(data_space, H5S_SELECT_SET, offset, nullptr, count,
                        nullptr);
    H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, mem_offset, nullptr, count,
                        nullptr);
    if (H5Dwrite(data_set, struct_tid, mem_space, data_space, H5P_DEFAULT,
                 &local_data) < 0) {
      throw std::string{"Failed to write static plane struct to HDF"};
    }
  }
  assert(current_plane == static_planes.size());
#else
  throw std::string{"writeStaticPlanes not compiled with HDF5 support"};
#endif
}

void RigidBody2DStateOutput::writeStaticDrums(
    const std::vector<RigidBody2DStaticDrum> &static_drums,
    const std::string &group, HDF5File &output_file) {
#ifdef USE_HDF5
  using HDFTID = HDFID<H5Tclose>;
  using HDFSID = HDFID<H5Sclose>;
  using HDFGID = HDFID<H5Gclose>;
  using HDFDID = HDFID<H5Dclose>;

  struct LocalStaticDrumData final {
    scalar x[2];
    scalar r;
    scalar theta;
  };

  // Create an HDF5 dataspace
  const hsize_t dim[]{static_drums.size()};
  const HDFSID data_space{H5Screate_simple(1, dim, nullptr)};
  if (data_space < 0) {
    throw std::string{"Failed to create HDF dataspace for static drums"};
  }

  // Create an HDF5 struct for the data
  const HDFTID struct_tid{H5Tcreate(H5T_COMPOUND, sizeof(LocalStaticDrumData))};
  if (struct_tid < 0) {
    throw std::string{"Failed to create HDF struct for static drums"};
  }
  // Insert the x type in the struct
  {
    const hsize_t array_dim[]{2};
    const HDFTID array_tid{H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, array_dim)};
    if (array_tid < 0) {
      throw std::string{"Failed to create HDF x type for static drums"};
    }
    if (H5Tinsert(struct_tid, "x", HOFFSET(LocalStaticDrumData, x), array_tid) <
        0) {
      throw std::string{"Failed to insert x in HDF struct for static drums"};
    }
  }
  // Insert the r type in the struct
  if (H5Tinsert(struct_tid, "r", HOFFSET(LocalStaticDrumData, r),
                H5T_NATIVE_DOUBLE) < 0) {
    throw std::string{"Failed to create HDF r type for static drums"};
  }
  // Insert the theta type in the struct
  if (H5Tinsert(struct_tid, "theta", HOFFSET(LocalStaticDrumData, theta),
                H5T_NATIVE_DOUBLE) < 0) {
    throw std::string{"Failed to create HDF theta type for static drums"};
  }

  const HDFGID grp_id{output_file.getGroup(group)};

  // Create an HDF5 dataset
  const HDFDID data_set{H5Dcreate2(grp_id, "static_drums", struct_tid,
                                   data_space, H5P_DEFAULT, H5P_DEFAULT,
                                   H5P_DEFAULT)};
  if (data_set < 0) {
    throw std::string{"Failed to create HDF dataset for static drums"};
  }

  // Create an HDF5 memspace to allow us to insert elements one by one
  const HDFSID mem_space{H5Screate_simple(1, dim, nullptr)};
  if (mem_space < 0) {
    throw std::string{"Failed to create HDF memspace for static drums"};
  }

  // Insert the static planes one by one
  unsigned current_drum{0};
  LocalStaticDrumData local_drum;
  for (const RigidBody2DStaticDrum &drum : static_drums) {
    Eigen::Map<Vector2s>{local_drum.x} = drum.x();
    local_drum.r = drum.r();
    local_drum.theta = drum.theta();
    const hsize_t count[]{1};
    const hsize_t offset[]{current_drum++};
    const hsize_t mem_offset[]{0};
    H5Sselect_hyperslab(data_space, H5S_SELECT_SET, offset, nullptr, count,
                        nullptr);
    H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, mem_offset, nullptr, count,
                        nullptr);
    if (H5Dwrite(data_set, struct_tid, mem_space, data_space, H5P_DEFAULT,
                 &local_drum) < 0) {
      throw std::string{"Failed to write static drum struct to HDF"};
    }
  }
  assert(current_drum == static_drums.size());
#else
  throw std::string{"writeStaticDrums not compiled with HDF5 support"};
#endif
}

void RigidBody2DStateOutput::writePlanarPortals(
    const std::vector<PlanarPortal> &planar_portals, const std::string &group,
    HDF5File &output_file) {
#ifdef USE_HDF5
  using HDFTID = HDFID<H5Tclose>;
  using HDFSID = HDFID<H5Sclose>;
  using HDFGID = HDFID<H5Gclose>;
  using HDFDID = HDFID<H5Dclose>;

  struct LocalPlanarPortalData final {
    // Data for the first plane
    scalar x0[2];
    scalar n0[2];
    // Data for the second plane
    scalar x1[2];
    scalar n1[2];
    // Lees-Edwards data
    scalar v;
    scalar width;
    scalar delta;
  };

  // Create an HDF5 dataspace
  const hsize_t dim[]{planar_portals.size()};
  const HDFSID data_space{H5Screate_simple(1, dim, nullptr)};
  if (data_space < 0) {
    throw std::string{"Failed to create HDF dataspace for planar portals"};
  }

  // Create an HDF5 struct for the data
  const HDFTID struct_tid{
      H5Tcreate(H5T_COMPOUND, sizeof(LocalPlanarPortalData))};
  if (struct_tid < 0) {
    throw std::string{"Failed to create HDF struct for planar portals"};
  }
  // Insert the x0 type in the struct
  {
    const hsize_t array_dim[]{2};
    const HDFTID array_tid{H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, array_dim)};
    if (array_tid < 0) {
      throw std::string{"Failed to create HDF x0 type for planar portals"};
    }
    if (H5Tinsert(struct_tid, "x0", HOFFSET(LocalPlanarPortalData, x0),
                  array_tid) < 0) {
      throw std::string{"Failed to insert x0 in HDF struct for planar portals"};
    }
  }
  // Insert the n0 type in the struct
  {
    const hsize_t array_dim[]{2};
    const HDFTID array_tid{H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, array_dim)};
    if (array_tid < 0) {
      throw std::string{"Failed to create HDF n0 type for planar portals"};
    }
    if (H5Tinsert(struct_tid, "n0", HOFFSET(LocalPlanarPortalData, n0),
                  array_tid) < 0) {
      throw std::string{"Failed to insert n0 in HDF struct for planar portals"};
    }
  }
  // Insert the x1 type in the struct
  {
    const hsize_t array_dim[]{2};
    const HDFTID array_tid{H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, array_dim)};
    if (array_tid < 0) {
      throw std::string{"Failed to create HDF x1 type for planar portals"};
    }
    if (H5Tinsert(struct_tid, "x1", HOFFSET(LocalPlanarPortalData, x1),
                  array_tid) < 0) {
      throw std::string{"Failed to insert x1 in HDF struct for planar portals"};
    }
  }
  // Insert the n1 type in the struct
  {
    const hsize_t array_dim[]{2};
    const HDFTID array_tid{H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, array_dim)};
    if (array_tid < 0) {
      throw std::string{"Failed to create HDF n1 type for planar portals"};
    }
    if (H5Tinsert(struct_tid, "n1", HOFFSET(LocalPlanarPortalData, n1),
                  array_tid) < 0) {
      throw std::string{"Failed to insert n1 in HDF struct for planar portals"};
    }
  }
  // Insert the v type in the struct
  if (H5Tinsert(struct_tid, "v", HOFFSET(LocalPlanarPortalData, v),
                H5T_NATIVE_DOUBLE) < 0) {
    throw std::string{"Failed to create HDF v type for planar portals"};
  }
  // Insert the width type in the struct
  if (H5Tinsert(struct_tid, "width", HOFFSET(LocalPlanarPortalData, width),
                H5T_NATIVE_DOUBLE) < 0) {
    throw std::string{"Failed to create HDF width type for planar portals"};
  }
  // Insert the delta type in the struct
  if (H5Tinsert(struct_tid, "delta", HOFFSET(LocalPlanarPortalData, delta),
                H5T_NATIVE_DOUBLE) < 0) {
    throw std::string{"Failed to create HDF delta type for planar portals"};
  }

  const HDFGID grp_id{output_file.getGroup(group)};

  // Create an HDF5 dataset
  const HDFDID data_set{H5Dcreate2(grp_id, "planar_portals", struct_tid,
                                   data_space, H5P_DEFAULT, H5P_DEFAULT,
                                   H5P_DEFAULT)};
  if (data_set < 0) {
    throw std::string{"Failed to create HDF dataset for planar portals"};
  }

  // Create an HDF5 memspace to allow us to insert elements one by one
  const HDFSID mem_space{H5Screate_simple(1, dim, nullptr)};
  if (mem_space < 0) {
    throw std::string{"Failed to create HDF memspace for planar portals"};
  }

  // Insert the planar portals one by one
  unsigned current_portal{0};
  LocalPlanarPortalData local_data;
  for (const PlanarPortal &planar_portal : planar_portals) {
    Eigen::Map<Vector2s>{local_data.x0} = planar_portal.planeA().x();
    Eigen::Map<Vector2s>{local_data.n0} = planar_portal.planeA().n();
    Eigen::Map<Vector2s>{local_data.x1} = planar_portal.planeB().x();
    Eigen::Map<Vector2s>{local_data.n1} = planar_portal.planeB().n();
    local_data.v = planar_portal.v();
    local_data.width = planar_portal.bounds();
    local_data.delta = planar_portal.delta();
    const hsize_t count[]{1};
    const hsize_t offset[]{current_portal++};
    const hsize_t mem_offset[]{0};
    H5Sselect_hyperslab(data_space, H5S_SELECT_SET, offset, nullptr, count,
                        nullptr);
    H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, mem_offset, nullptr, count,
                        nullptr);
    if (H5Dwrite(data_set, struct_tid, mem_space, data_space, H5P_DEFAULT,
                 &local_data) < 0) {
      throw std::string{"Failed to write planar_portal struct to HDF"};
    }
  }
  assert(current_portal == planar_portals.size());
#else
  throw std::string{"writePlanarPortals not compiled with HDF5 support"};
#endif
}
