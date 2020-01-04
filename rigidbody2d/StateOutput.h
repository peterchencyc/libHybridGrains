// StateOutput.h
//
// Breannan Smith
// Last updated: 10/01/2015

#ifndef RIGID_BODY_2D_STATE_OUTPUT_H
#define RIGID_BODY_2D_STATE_OUTPUT_H

// TODO: Template on float/double type, add cmake option to dump contents as
// float

#include <iosfwd>
#include <memory>
#include <vector>

#include "scisim/Math/MathDefines.h"

class RigidBody2DGeometry;
class HDF5File;
class RigidBody2DStaticDrum;
class RigidBody2DStaticPlane;
class PlanarPortal;

namespace RigidBody2DStateOutput {

#ifndef USE_HDF5
[[noreturn]]
#endif
  void writeGeometryIndices( const std::vector<std::unique_ptr<RigidBody2DGeometry>>& geometry, const VectorXu& indices, const std::string& group, HDF5File& output_file );

void writeGeometry(
    const std::vector<std::unique_ptr<RigidBody2DGeometry>> &geometry,
    const std::string &group, HDF5File &output_file);

#ifndef USE_HDF5
[[noreturn]]
#endif
  void writeStaticPlanes( const std::vector<RigidBody2DStaticPlane>& static_planes, const std::string& group, HDF5File& output_file );

#ifndef USE_HDF5
[[noreturn]]
#endif
  void writeStaticDrums( const std::vector<RigidBody2DStaticDrum>& static_drums, const std::string& group, HDF5File& output_file );

#ifndef USE_HDF5
[[noreturn]]
#endif
  void writePlanarPortals( const std::vector<PlanarPortal>& planar_portals, const std::string& group, HDF5File& output_file );

} // namespace RigidBody2DStateOutput

#endif
