#include "ConstrainedSystem.h"

#include <iostream>

ConstrainedSystem::~ConstrainedSystem() = default;

void ConstrainedSystem::computeActiveSetNew(
    const VectorXs &q, const VectorXs &v,
    std::vector<std::vector<std::unique_ptr<Constraint>>> &con_table) {
  std::cerr << "Error, computeActiveSetNew not implemented for this "
               "ConstrainedSystem implementation. Exiting."
            << std::endl;
  std::exit(EXIT_FAILURE);
}
