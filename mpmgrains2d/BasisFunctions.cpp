#include "BasisFunctions.h"

#include "scisim/Utilities.h"

#include "CubicBasisFunctions.h"
#include "LinearBasisFunctions.h"
#include "uGIMPLinearBasisFunctions.h"

#ifndef CMAKE_DETECTED_CLANG_COMPILER
#include <iostream>
#endif

BasisFunctions::~BasisFunctions() = default;

void BasisFunctionsTools::serialize(
    const std::unique_ptr<BasisFunctions> &basis_functions,
    std::ostream &output_stream) {
  Utilities::serializeBuiltInType(basis_functions->type(), output_stream);
}

std::unique_ptr<BasisFunctions>
BasisFunctionsTools::deserialize(std::istream &input_stream) {
  const BasisFunctionType basis_type{
      Utilities::deserialize<BasisFunctionType>(input_stream)};
  switch (basis_type) {
  case BasisFunctionType::Linear:
    return std::make_unique<LinearBasisFunctions>();
  case BasisFunctionType::ThirdOrder:
    return std::make_unique<ThirdOrderBasisFunctions>();
  case BasisFunctionType::uGIMPLinear:
    return std::make_unique<uGIMPLinearBasisFunctions>();
#ifndef CMAKE_DETECTED_CLANG_COMPILER
  default: {
    std::cerr << "Invalid basis type in BasisFunctionsTools::deserialize, this "
                 "is a bug."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif
  }
}
