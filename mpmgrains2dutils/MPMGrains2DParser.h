#ifndef MPM_GRAINS_2D_PARSER_H
#define MPM_GRAINS_2D_PARSER_H

#include <string>

struct InitialSimulationState;

namespace MPMGrains2DParser {

bool readXMLFile(const std::string &file_name,
                 InitialSimulationState &initial_state);

}

#endif
