// RigidBody2DUtilities.cpp
//
// Breannan Smith
// Last updated: 12/08/2015

#include "RigidBody2DUtilities.h"

#include <cassert>
#include <iostream>

#include "scisim/StringUtilities.h"
#include "scisim/UnconstrainedMaps/UnconstrainedMap.h"

void RigidBody2DUtilities::serialize(
    const std::unique_ptr<UnconstrainedMap> &unconstrained_map,
    std::ostream &output_stream) {
  assert(output_stream.good());

  if (unconstrained_map != nullptr) {
    StringUtilities::serializeString(unconstrained_map->name(), output_stream);
    unconstrained_map->serialize(output_stream);
  } else {
    StringUtilities::serializeString("NULL", output_stream);
  }
}
