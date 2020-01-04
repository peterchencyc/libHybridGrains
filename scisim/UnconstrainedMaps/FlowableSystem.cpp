#include "FlowableSystem.h"

#include <iostream>

FlowableSystem::~FlowableSystem() = default;

unsigned FlowableSystem::numBodies() const {
  return nvdofs() / numVelDoFsPerBody();
}

void FlowableSystem::zeroOutForcesOnFixedBodies(VectorXs &F) const {
  std::cerr << "FlowableSystem::zeroOutForcesOnFixedBodies not coded up for "
               "this simulation type."
            << std::endl;
  std::exit(EXIT_FAILURE);
}

const std::vector<bool> &FlowableSystem::fixed() const {
  std::cerr << "FlowableSystem::fixed not coded up for this simulation type."
            << std::endl;
  std::exit(EXIT_FAILURE);
}
