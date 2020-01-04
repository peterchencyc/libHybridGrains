#include "HybridGrains2DSceneParser.h"

#include <fstream>
#include <iostream>
#include <vector>

#include "rapidxml.hpp"

#include "hybridgrains2dnew/HybridDefinitions.h"
#include "hybridgrains2dnew/HybridGrains2DSim.h"
#include "scisim/Math/Rational.h"
#include "scisim/StringUtilities.h"

static bool loadTextFileIntoVector(const std::string &filename,
                                   std::vector<char> &xmlchars) {
  assert(xmlchars.empty());

  // Attempt to open the text file for reading
  std::ifstream textfile(filename);
  if (!textfile.is_open()) {
    return false;
  }

  // Read the entire file into a single string
  std::string line;
  while (getline(textfile, line)) {
    std::copy(line.begin(), line.end(), back_inserter(xmlchars));
  }
  xmlchars.emplace_back('\0');

  return true;
}

static bool loadXMLFile(const std::string &filename,
                        std::vector<char> &xmlchars,
                        rapidxml::xml_document<> &doc) {
  assert(xmlchars.empty());

  // Attempt to read the text from the user-specified xml file
  if (!loadTextFileIntoVector(filename, xmlchars)) {
    std::cerr << "Failed to read scene file: " << filename << std::endl;
    return false;
  }

  // Initialize the xml parser with the character vector
  try {
    doc.parse<0>(xmlchars.data());
  } catch (const rapidxml::parse_error &e) {
    std::cerr << "Failed to parse scene file: " << filename << std::endl;
    std::cerr << "Error message: " << e.what() << std::endl;
    return false;
  }

  return true;
}

static bool loadDiscreteFileName(const rapidxml::xml_node<> &node,
                                 std::string &discrete_file_name) {
  const rapidxml::xml_node<> *const discrete_sim_node{
      node.first_node("discrete_sim")};
  if (discrete_sim_node == nullptr) {
    std::cerr << "Failed to locate discrete_sim node." << std::endl;
    return false;
  }

  const rapidxml::xml_attribute<> *const file_name_attrib{
      discrete_sim_node->first_attribute("file_name")};
  if (file_name_attrib == nullptr) {
    std::cerr << "Failed to locate file_name attribute for discrete_sim node."
              << std::endl;
    return false;
  }
  discrete_file_name = file_name_attrib->value();

  return true;
}

static bool loadContinuumFileName(const rapidxml::xml_node<> &node,
                                  std::string &continuum_file_name) {
  const rapidxml::xml_node<> *const continuum_sim_node{
      node.first_node("continuum_sim")};
  if (continuum_sim_node == nullptr) {
    std::cerr << "Failed to locate continuum_sim node." << std::endl;
    return false;
  }

  const rapidxml::xml_attribute<> *const file_name_attrib{
      continuum_sim_node->first_attribute("file_name")};
  if (file_name_attrib == nullptr) {
    std::cerr << "Failed to locate file_name attribute for continuum_sim node."
              << std::endl;
    return false;
  }
  continuum_file_name = file_name_attrib->value();

  return true;
}

static bool loadTimeStep(const rapidxml::xml_node<> &node,
                         Rational<std::intmax_t> &time_step,
                         std::string &time_step_string) {
  const rapidxml::xml_node<> *const time_step_node{
      node.first_node("time_step")};
  if (time_step_node == nullptr) {
    std::cerr << "Failed to locate time_step node." << std::endl;
    return false;
  }
  const rapidxml::xml_attribute<> *dt_attrib{
      time_step_node->first_attribute("dt")};
  if (!dt_attrib) {
    std::cerr << "Failed to locate dt attribute for time_step node."
              << std::endl;
    return false;
  }
  if (!RationalTools::extractFromString(std::string{dt_attrib->value()},
                                        time_step) ||
      !time_step.positive()) {
    std::cerr << "Failed to load dt attribute for time_step. Must provide a "
                 "positive number."
              << std::endl;
    return false;
  }
  time_step_string = dt_attrib->value();
  return true;
}

static bool loadEndTime(const rapidxml::xml_node<> &node,
                        Rational<std::intmax_t> &end_time) {
  const rapidxml::xml_node<> *const end_time_node{node.first_node("end_time")};
  if (end_time_node == nullptr) {
    std::cerr << "Failed to locate end_time node." << std::endl;
    return false;
  }
  const rapidxml::xml_attribute<> *end_time_attrib{
      end_time_node->first_attribute("t")};
  if (!end_time_attrib) {
    std::cerr << "Failed to locate t attribute for end_time node." << std::endl;
    return false;
  }
  if (!RationalTools::extractFromString(std::string{end_time_attrib->value()},
                                        end_time) ||
      !end_time.positive()) {
    std::cerr << "Failed to load t attribute for end_time. Must provide a "
                 "positive number."
              << std::endl;
    return false;
  }
  return true;
}

static bool
loadHybridIntegratorSettings(const rapidxml::xml_node<> &node,
                             HybridIntegratorSettings &hybrid_settings) {
  const rapidxml::xml_node<> *const coupling_force_node{
      node.first_node("coupling_force")};

  if (coupling_force_node == nullptr) {
    hybrid_settings.poorman_settings.enabled = false;
    return true;
  }

  hybrid_settings.poorman_settings.enabled = true;

  {
    hybrid_settings.poorman_settings.rzone_half_thickness = 0.0;
    const rapidxml::xml_attribute<> *rzone_half_thickness_attrib{
        coupling_force_node->first_attribute("rzoneHalfThickness")};
    if (rzone_half_thickness_attrib) {
      StringUtilities::extractFromString(
          rzone_half_thickness_attrib->value(),
          hybrid_settings.poorman_settings.rzone_half_thickness);
    }

    hybrid_settings.poorman_settings.dem_r_mean = -1.0;
    const rapidxml::xml_attribute<> *dem_r_mean_attrib{
        coupling_force_node->first_attribute("demRMean")};
    if (dem_r_mean_attrib) {
      StringUtilities::extractFromString(
          dem_r_mean_attrib->value(),
          hybrid_settings.poorman_settings.dem_r_mean);
    }

    hybrid_settings.poorman_settings.dem_r_std = 0.0;
    const rapidxml::xml_attribute<> *dem_r_std_attrib{
        coupling_force_node->first_attribute("demRStd")};
    if (dem_r_std_attrib) {
      StringUtilities::extractFromString(
          dem_r_std_attrib->value(),
          hybrid_settings.poorman_settings.dem_r_std);
    }

    const rapidxml::xml_attribute<> *phi_window_size_attrib{
        coupling_force_node->first_attribute("phi_window_size")};
    if (phi_window_size_attrib) {
      StringUtilities::extractFromString(
          std::string{phi_window_size_attrib->value()},
          hybrid_settings.poorman_settings.phi_window_size);
    }

    const rapidxml::xml_attribute<> *level_set_cell_width_attrib{
        coupling_force_node->first_attribute("level_set_cell_width")};
    if (level_set_cell_width_attrib) {
      StringUtilities::extractFromString(
          std::string{level_set_cell_width_attrib->value()},
          hybrid_settings.poorman_settings.level_set_cell_width);
    }

    const rapidxml::xml_attribute<> *phi_samples_per_cell_side_attrib{
        coupling_force_node->first_attribute("phi_samples_per_cell_side")};
    if (phi_samples_per_cell_side_attrib) {
      StringUtilities::extractFromString(
          std::string{phi_samples_per_cell_side_attrib->value()},
          hybrid_settings.poorman_settings.phi_samples_per_cell_side);
    }

    const rapidxml::xml_attribute<> *target_dense_packing_fraction_attrib{
        coupling_force_node->first_attribute("target_dense_packing_fraction")};
    if (target_dense_packing_fraction_attrib) {
      StringUtilities::extractFromString(
          std::string{target_dense_packing_fraction_attrib->value()},
          hybrid_settings.poorman_settings.target_dense_packing_fraction);
    }

    hybrid_settings.poorman_settings.avoidAVoidFreq = 20;
    const rapidxml::xml_attribute<> *avoid_a_void_freq_attrib{
        coupling_force_node->first_attribute("avoidAVoidFreq")};
    if (avoid_a_void_freq_attrib) {
      StringUtilities::extractFromString(
          avoid_a_void_freq_attrib->value(),
          hybrid_settings.poorman_settings.avoidAVoidFreq);
    }

    hybrid_settings.poorman_settings.phi_threshold = 0.5;
    const rapidxml::xml_attribute<> *phi_threshold_attrib{
        coupling_force_node->first_attribute("phiThreshold")};
    if (phi_threshold_attrib) {
      StringUtilities::extractFromString(
          phi_threshold_attrib->value(),
          hybrid_settings.poorman_settings.phi_threshold);
    }

    hybrid_settings.poorman_settings.rzone_level_set = 0.0;
    const rapidxml::xml_attribute<> *rzone_level_set_attrib{
        coupling_force_node->first_attribute("rzoneLevelSet")};
    if (rzone_level_set_attrib) {
      StringUtilities::extractFromString(
          rzone_level_set_attrib->value(),
          hybrid_settings.poorman_settings.rzone_level_set);
    }

    hybrid_settings.poorman_settings.rho_dem = 0.96;
    const rapidxml::xml_attribute<> *rho_dem_attrib{
        coupling_force_node->first_attribute("rho_dem")};
    if (rho_dem_attrib) {
      StringUtilities::extractFromString(
          rho_dem_attrib->value(), hybrid_settings.poorman_settings.rho_dem);
    }

    hybrid_settings.poorman_settings
        .allow_direct_transitions_between_discrete_and_continuum = false;
    const rapidxml::xml_attribute<> *allow_direct_transitions_attrib{
        coupling_force_node->first_attribute(
            "allow_direct_transitions_between_discrete_and_continuum")};
    if (allow_direct_transitions_attrib) {
      int val;
      StringUtilities::extractFromString(
          std::string{allow_direct_transitions_attrib->value()}, val);
      if (val != 0)
        hybrid_settings.poorman_settings
            .allow_direct_transitions_between_discrete_and_continuum = true;
    }

    hybrid_settings.poorman_settings.newly_converted_continuum_stress_free =
        false;
    const rapidxml::xml_attribute<>
        *newly_converted_continuum_stress_free_attrib{
            coupling_force_node->first_attribute(
                "newly_converted_continuum_stress_free")};
    if (newly_converted_continuum_stress_free_attrib) {
      int val;
      StringUtilities::extractFromString(
          std::string{newly_converted_continuum_stress_free_attrib->value()},
          val);
      if (val != 0) {
        hybrid_settings.poorman_settings.newly_converted_continuum_stress_free =
            true;
        std::cout << "debugging: newly_converted_continuum_stress_free enabled"
                  << std::endl;
      }
    }

    hybrid_settings.poorman_settings.homogenize_stress = false;
    const rapidxml::xml_attribute<> *homogenize_stress_attrib{
        coupling_force_node->first_attribute("homogenize_stress")};
    if (homogenize_stress_attrib) {
      int val;
      StringUtilities::extractFromString(
          std::string{homogenize_stress_attrib->value()}, val);
      if (val != 0) {
        hybrid_settings.poorman_settings.homogenize_stress = true;
        std::cout << "debugging: homogenize_stress enabled" << std::endl;
      }
    }

    hybrid_settings.poorman_settings.grid_smoothing_homogenized_stress = false;
    const rapidxml::xml_attribute<> *grid_smoothing_homogenized_stress_attrib{
        coupling_force_node->first_attribute(
            "grid_smoothing_homogenized_stress")};
    if (grid_smoothing_homogenized_stress_attrib) {
      int val;
      StringUtilities::extractFromString(
          std::string{grid_smoothing_homogenized_stress_attrib->value()}, val);
      if (val != 0) {
        hybrid_settings.poorman_settings.grid_smoothing_homogenized_stress =
            true;
        std::cout << "debugging: grid_smoothing_homogenized_stress enabled"
                  << std::endl;
      }
    }

    hybrid_settings.poorman_settings.homogenize_velocity = false;
    const rapidxml::xml_attribute<> *homogenize_velocity_attrib{
        coupling_force_node->first_attribute("homogenize_velocity")};
    if (homogenize_velocity_attrib) {
      int val;
      StringUtilities::extractFromString(
          std::string{homogenize_velocity_attrib->value()}, val);
      if (val != 0) {
        hybrid_settings.poorman_settings.homogenize_velocity = true;
        std::cout << "debugging: homogenize_velocity enabled" << std::endl;
      }
    }

    hybrid_settings.poorman_settings.grid_smoothing_homogenized_velocity =
        false;
    const rapidxml::xml_attribute<> *grid_smoothing_homogenized_velocity_attrib{
        coupling_force_node->first_attribute(
            "grid_smoothing_homogenized_velocity")};
    if (grid_smoothing_homogenized_velocity_attrib) {
      int val;
      StringUtilities::extractFromString(
          std::string{grid_smoothing_homogenized_velocity_attrib->value()},
          val);
      if (val != 0) {
        hybrid_settings.poorman_settings.grid_smoothing_homogenized_velocity =
            true;
        std::cout << "debugging: grid_smoothing_homogenized_velocity enabled"
                  << std::endl;
      }
    }
  }

  return true;
}

static bool
loadIntegratorStyle(const rapidxml::xml_node<> &node,
                    HybridIntegratorState::IntegratorStyle &integrator_style) {
  const rapidxml::xml_node<> *const style_node{
      node.first_node("integrator_style")};
  if (style_node == nullptr) {
    std::cerr << "Failed to locate integrator_style node." << std::endl;
    return false;
  }

  const rapidxml::xml_attribute<> *const style_attrib{
      style_node->first_attribute("style")};
  if (style_attrib == nullptr) {
    std::cerr << "Failed to locate style attribute for integrator_style node."
              << std::endl;
    return false;
  }

  const std::string style_string{style_attrib->value()};

  if (style_string == "old_version") {
    integrator_style = HybridIntegratorState::IntegratorStyle::OLD_VERSION;
  } else if (style_string == "old_version_node_node") {
    integrator_style =
        HybridIntegratorState::IntegratorStyle::OLD_VERSION_NODE_NODE;
  } else if (style_string == "prediction_correction_advection") {
    integrator_style =
        HybridIntegratorState::IntegratorStyle::PREDICTION_CORRECTION_ADVECTION;
  } else if (style_string == "prediction_correction_advection_node_node") {
    integrator_style = HybridIntegratorState::IntegratorStyle::
        PREDICTION_CORRECTION_ADVECTION_NODE_NODE;
  } else if (style_string == "iterative") {
    integrator_style = HybridIntegratorState::IntegratorStyle::ITERATIVE;
  } else {
    std::cerr
        << "Invalid style value provided for integrator_style node. Valid "
           "values are: old_version, prediction_correction_advection, iterative"
        << std::endl;
    return false;
  }

  return true;
}

static bool loadHybridKinematicScript(
    const rapidxml::xml_node<> &node,
    HybridIntegratorState::HybridKinematicScript &hybrid_kinematic_script) {
  const rapidxml::xml_node<> *const hybrid_kinematic_script_node{
      node.first_node("hybrid_kinematic_script")};
  if (hybrid_kinematic_script_node == nullptr) {
    std::cerr << "Failed to locate hybrid_kinematic_script node." << std::endl;
    return false;
  }

  const rapidxml::xml_attribute<> *const hybrid_kinematic_script_attrib{
      hybrid_kinematic_script_node->first_attribute("type")};
  if (hybrid_kinematic_script_attrib == nullptr) {
    std::cerr
        << "Failed to locate type attribute for hybrid_kinematic_script node."
        << std::endl;
    return false;
  }

  const std::string type_string{hybrid_kinematic_script_attrib->value()};

  if (type_string == "none") {
    hybrid_kinematic_script =
        HybridIntegratorState::HybridKinematicScript::NONE;
  } else if (type_string == "script_hybrid_fronts") {
    hybrid_kinematic_script =
        HybridIntegratorState::HybridKinematicScript::SCRIPT_HYBRID_FRONTS;
  } else {
    std::cerr << "Invalid type value provided for hybrid_kinematic_script "
                 "node. Valid values are: none, script_hybrid_fronts"
              << std::endl;
    return false;
  }

  return true;
}

bool HybridGrains2DSceneParser::parseXMLSceneFile(
    const std::string &file_name, HybridIntegratorSettings &hybrid_settings) {
  // Attempt to load the xml document
  std::vector<char> xmlchars;
  rapidxml::xml_document<> doc;
  if (!loadXMLFile(file_name, xmlchars, doc)) {
    return false;
  }

  // Attempt to locate the root node
  if (doc.first_node("hybridgrains2d_scene") == nullptr) {
    std::cerr
        << "Failed to locate root node hybridgrains2d_scene in xml scene file: "
        << file_name << std::endl;
    return false;
  }
  const rapidxml::xml_node<> &root_node{
      *doc.first_node("hybridgrains2d_scene")};

  HybridIntegratorSettings new_integrator_settings;

  // Attempt to load the discrete scene name
  if (!loadDiscreteFileName(root_node,
                            new_integrator_settings.discrete_file_name)) {
    return false;
  }

  // Attempt to load the continuum scene name
  if (!loadContinuumFileName(root_node,
                             new_integrator_settings.continuum_file_name)) {
    return false;
  }

  // Attempt to load the timestep
  if (!loadTimeStep(root_node, new_integrator_settings.overall_dt,
                    new_integrator_settings.time_step_string)) {
    return false;
  }

  // Attempt to load the end time
  if (!loadEndTime(root_node, new_integrator_settings.end_time)) {
    return false;
  }

  if (!loadHybridIntegratorSettings(root_node, new_integrator_settings)) {
    return false;
  }

  // Attempt to load the integration type
  // if( !loadIntegratorType( root_node, new_integrator_settings.integrator_type
  // ) )
  // {
  //   return false;
  // }

  // Attempt to load the integration style
  if (!loadIntegratorStyle(root_node,
                           new_integrator_settings.integrator_style)) {
    return false;
  }

  // Attempt to load the hybrid kinematic script settings
  if (!loadHybridKinematicScript(
          root_node,
          new_integrator_settings.kinematically_scripted_hybrid_fronts)) {
    return false;
  }

  // Append the scene names to the input path
  std::string input_path;
  std::string input_name;
  StringUtilities::splitAtLastCharacterOccurence(file_name, input_path,
                                                 input_name, '/');
  new_integrator_settings.discrete_file_name =
      !input_name.empty()
          ? input_path + '/' + new_integrator_settings.discrete_file_name
          : new_integrator_settings.discrete_file_name;
  new_integrator_settings.continuum_file_name =
      !input_name.empty()
          ? input_path + '/' + new_integrator_settings.continuum_file_name
          : new_integrator_settings.continuum_file_name;

  hybrid_settings = std::move(new_integrator_settings);

  return true;
}
