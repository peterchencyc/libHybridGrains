#ifndef HYBRID_GRAINS_2D_SCENE_PARSER_H
#define HYBRID_GRAINS_2D_SCENE_PARSER_H

#include <string>
template <typename T> class Rational;
struct HybridIntegratorSettings;

namespace HybridGrains2DSceneParser {

bool parseXMLSceneFile(const std::string &file_name,
                       HybridIntegratorSettings &hybrid_settings);

}

#endif
