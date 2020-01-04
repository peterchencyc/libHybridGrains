#include "MPMGrains2DParser.h"

#include "rapidxml.hpp"

#include <fstream>
#include <iostream>

#include "mpmgrains2d/InitialSimulationState.h"
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
    std::copy(line.cbegin(), line.cend(), back_inserter(xmlchars));
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

static bool loadGravityForce(const rapidxml::xml_node<> &node,
                             InitialSimulationState &state) {
  state.near_earth_gravity.setZero();

  for (rapidxml::xml_node<> *nd = node.first_node("near_earth_gravity"); nd;
       nd = nd->next_sibling("near_earth_gravity")) {
    const rapidxml::xml_attribute<> *const f_attrib{nd->first_attribute("f")};
    if (f_attrib == nullptr) {
      std::cerr << "Failed to locate f attribute for near_earth_gravity node."
                << std::endl;
      return false;
    }
    Vector2s f;
    if (!StringUtilities::readScalarList(f_attrib->value(), 2, ' ', f)) {
      std::cerr << "Failed to load f attribute for near_earth_gravity node, "
                   "must provide two scalars."
                << std::endl;
      return false;
    }
    state.near_earth_gravity = f;
  }

  return true;
}

static bool loadIntegrator(const rapidxml::xml_node<> &node,
                           InitialSimulationState &state) {
  // Attempt to locate the integrator node
  const rapidxml::xml_node<> *nd{node.first_node("integrator")};
  if (nd == nullptr) {
    std::cerr << "Failed to locate integrator node." << std::endl;
    return false;
  }

  // Attempt to load the timestep
  std::string dt_string;
  Rational<std::intmax_t> dt;
  {
    const rapidxml::xml_attribute<> *dtnd{nd->first_attribute("dt")};
    if (dtnd == nullptr) {
      std::cerr << "Failed to locate dt attribute for integrator node."
                << std::endl;
      return false;
    }
    if (!RationalTools::extractFromString(std::string{dtnd->value()}, dt) ||
        !dt.positive()) {
      std::cerr << "Failed to load dt attribute for integrator. Must provide a "
                   "positive number."
                << std::endl;
      return false;
    }
    dt_string = dtnd->value();
  }

  // Attempt to load the constant that controls the PIC/FLIP interpolation
  scalar alpha;
  {
    const rapidxml::xml_attribute<> *attrib{nd->first_attribute("alpha")};
    if (attrib == nullptr) {
      std::cerr << "Failed to locate alpha attribute for integrator node."
                << std::endl;
      return false;
    }
    std::stringstream ss;
    ss << attrib->value();
    if (!(ss >> alpha)) {
      std::cerr << "Failed to load alpha attribute for integrator. Must "
                   "provide a scalar."
                << std::endl;
      return false;
    }
    if (alpha < 0.0 || alpha > 1.0) {
      std::cerr << "Failed to load alpha attribute for integrator. Must "
                   "provide a scalar in [0, 1]."
                << std::endl;
      return false;
    }
  }

  // Attempt to load the order of the basis functions
  unsigned order;
  {
    const rapidxml::xml_attribute<> *attrib{
        nd->first_attribute("shape_function_order")};
    if (attrib == nullptr) {
      std::cerr << "Failed to locate shape_function_order attribute for "
                   "integrator node."
                << std::endl;
      return false;
    }
    std::stringstream ss;
    ss << attrib->value();
    int bf_order;
    if (!(ss >> bf_order)) {
      std::cerr << "Failed to load shape_function_order attribute for "
                   "integrator. Must provide a positive integer."
                << std::endl;
      return false;
    }
    if (bf_order <= 0) {
      std::cerr << "Failed to load shape_function_order attribute for "
                   "integrator. Must provide a positive integer."
                << std::endl;
      return false;
    }
    if (bf_order != 1 && bf_order != 3) {
      std::cerr << "Failed to load shape_function_order attribute for "
                   "integrator. Supported values are 1 (linear) and 3 (cubic)."
                << std::endl;
      return false;
    }
    order = unsigned(bf_order);
  }

  BasisFunctionCategory category;
  {
    const rapidxml::xml_attribute<> *attrib{
        nd->first_attribute("shape_function_category")};
    if (attrib == nullptr) {
      // set the default basis category;
      category = BasisFunctionCategory::Standard;
    } else {
      std::stringstream ss;
      ss << attrib->value();

      if (strcmp(ss.str().c_str(), "uGIMP") == 0) {
        category = BasisFunctionCategory::uGIMP;
      } else if (strcmp(ss.str().c_str(), "Standard") == 0) {
        category = BasisFunctionCategory::Standard;
      } else {
        std::cerr << "Failed to load shape_function_category attribute. "
                     "Supported values are 'uGIMP' and 'Standard'."
                  << std::endl;
        return false;
      }
    }
  }

  state.dt_string = dt_string;
  state.dt = dt;
  state.alpha = alpha;
  state.basis_function_order = order;
  state.basis_function_category = category;

  return true;
}

static bool loadEndTime(const rapidxml::xml_node<> &node, scalar &end_time) {
  // Attempt to parse the time setting
  const rapidxml::xml_attribute<> *t_attrib{node.first_attribute("t")};
  if (!t_attrib) {
    std::cerr << "Failed to locate t attribute for end_time node." << std::endl;
    return false;
  }
  if (!StringUtilities::extractFromString(t_attrib->value(), end_time)) {
    std::cerr << "Failed to parse t attribute for end_time. Value must be a "
                 "positive scalar."
              << std::endl;
    return false;
  }
  if (end_time <= 0.0) {
    std::cerr << "Failed to parse t attribute for end_time. Value must be a "
                 "positive scalar."
              << std::endl;
    return false;
  }

  return true;
}

static bool
loadMaterialPointRectangles(const rapidxml::xml_node<> &node,
                            std::vector<RectangleRegion> &rectangle_regions) {
  for (rapidxml::xml_node<> *nd = node.first_node("rectangle"); nd;
       nd = nd->next_sibling("rectangle")) {
    // Attempt to read the rectangle's minimum
    Array2s min;
    {
      const rapidxml::xml_attribute<> *attrib{nd->first_attribute("min")};
      if (attrib == nullptr) {
        std::cerr << "Failed to locate min attribute for rectangle node."
                  << std::endl;
        return false;
      }
      if (!StringUtilities::readScalarList(attrib->value(), 2, ' ', min)) {
        std::cerr << "Failed to load min attribute for rectangle node, must "
                     "provide 2 scalars."
                  << std::endl;
        return false;
      }
    }

    // Attempt to read the rectangle's maximum
    Array2s max;
    {
      const rapidxml::xml_attribute<> *attrib{nd->first_attribute("max")};
      if (attrib == nullptr) {
        std::cerr << "Failed to locate min attribute for rectangle node."
                  << std::endl;
        return false;
      }
      if (!StringUtilities::readScalarList(attrib->value(), 2, ' ', max)) {
        std::cerr << "Failed to load min attribute for rectangle node, must "
                     "provide 2 scalars."
                  << std::endl;
        return false;
      }
    }

    if ((min >= max).any()) {
      std::cerr << "Failed to load rectangle region. Region max must be "
                   "greater than region min."
                << std::endl;
      return false;
    }

    // Attempt to read this bit of material's density
    scalar density;
    {
      const rapidxml::xml_attribute<> *attrib{nd->first_attribute("density")};
      if (attrib == nullptr) {
        std::cerr << "Failed to locate density attribute for rectangle node."
                  << std::endl;
        return false;
      }
      if (!StringUtilities::extractFromString(std::string{attrib->value()},
                                              density) ||
          density <= 0.0) {
        std::cerr << "Failed to load density attribute for rectangle node. "
                     "Must provide a positive scalar."
                  << std::endl;
        return false;
      }
    }

    // Attempt to read the number of samples along each dimension
    int samples_per_dim;
    {
      const rapidxml::xml_attribute<> *attrib{
          nd->first_attribute("cell_samples_per_dim")};
      if (attrib == nullptr) {
        std::cerr << "Failed to locate cell_samples_per_dim attribute for "
                     "rectangle node."
                  << std::endl;
        return false;
      }
      if (!StringUtilities::extractFromString(std::string{attrib->value()},
                                              samples_per_dim) ||
          samples_per_dim < 0) {
        std::cerr << "Failed to load cell_samples_per_dim attribute for "
                     "rectangle node. Must provide a positive integer."
                  << std::endl;
        return false;
      }
    }

    // Attempt to read the rectangle's initial velocity
    Vector2s vel;
    {
      const rapidxml::xml_attribute<> *attrib{nd->first_attribute("vel")};
      if (attrib == nullptr) {
        std::cerr << "Failed to locate vel attribute for rectangle node."
                  << std::endl;
        return false;
      }
      if (!StringUtilities::readScalarList(attrib->value(), 2, ' ', vel)) {
        std::cerr << "Failed to load vel attribute for rectangle node, must "
                     "provide 2 scalars."
                  << std::endl;
        return false;
      }
    }

    // Attempt to read the rectangle's initial angular velocity
    scalar omega;
    {
      const rapidxml::xml_attribute<> *attrib{nd->first_attribute("omega")};
      if (attrib == nullptr) {
        std::cerr << "Failed to locate omega attribute for rectangle node."
                  << std::endl;
        return false;
      }
      if (!StringUtilities::extractFromString(std::string{attrib->value()},
                                              omega)) {
        std::cerr << "Failed to load omega attribute for rectangle node, must "
                     "provide a scalar."
                  << std::endl;
        return false;
      }
    }

    rectangle_regions.emplace_back(min, max, density, unsigned(samples_per_dim),
                                   vel, omega);
  }

  return true;
}

static bool loadStaticPlanes(const rapidxml::xml_node<> &node,
                             std::vector<MPMStaticPlane> &planes) {
  for (rapidxml::xml_node<> *nd = node.first_node("static_plane"); nd;
       nd = nd->next_sibling("static_plane")) {
    // Read a point on the plane
    VectorXs x;
    {
      const rapidxml::xml_attribute<> *const attrib{nd->first_attribute("x")};
      if (!attrib) {
        std::cerr << "Failed to locate x attribute for static_plane."
                  << std::endl;
        return false;
      }
      if (!StringUtilities::readScalarList(attrib->value(), 2, ' ', x)) {
        std::cerr << "Failed to load x attribute for static_plane, must "
                     "provide 2 scalars."
                  << std::endl;
        return false;
      }
      assert(x.size() == 2);
    }
    // Read the normal
    VectorXs n;
    {
      const rapidxml::xml_attribute<> *const attrib{nd->first_attribute("n")};
      if (!attrib) {
        std::cerr << "Failed to locate n attribute for static_plane."
                  << std::endl;
        return false;
      }
      if (!StringUtilities::readScalarList(attrib->value(), 2, ' ', n)) {
        std::cerr << "Failed to load n attribute for static_plane, must "
                     "provide 2 scalars."
                  << std::endl;
        return false;
      }
      assert(n.size() == 2);
    }
    // Read the boundary treatment
    MPMStaticPlane::BoundaryBehavior boundary_behavior;
    {
      const rapidxml::xml_attribute<> *const attrib{
          nd->first_attribute("boundary_behavior")};
      if (!attrib) {
        std::cerr
            << "Failed to locate boundary_behavior attribute for static_plane."
            << std::endl;
        return false;
      }
      const std::string boundary_behavior_string{attrib->value()};
      if (boundary_behavior_string == "sticking") {
        boundary_behavior = MPMStaticPlane::BoundaryBehavior::STICKING;
      } else if (boundary_behavior_string == "sliding") {
        boundary_behavior = MPMStaticPlane::BoundaryBehavior::SLIDING;
      } else {
        std::cerr << "Invalid boundary_behavior attribute specified for "
                     "static_plane. Valid values are: sticking, sliding"
                  << std::endl;
        return false;
      }
    }
    //    read lower bound, useful for funnel
    scalar lower_bound;
    {
      const rapidxml::xml_attribute<> *const attrib{
          nd->first_attribute("lower_bound")};
      if (!attrib) {
        lower_bound = -SCALAR_INFINITY;
      } else {
        std::stringstream ss;
        ss << attrib->value();
        if (!(ss >> lower_bound)) {
          std::cerr << "Failed to load shear_modulus attribute for material. "
                       "Must provide a scalar."
                    << std::endl;
          return false;
        }
      }
    }
    planes.emplace_back(x, n, boundary_behavior, lower_bound);
  }

  return true;
}

static bool loadMaterialSettings(const rapidxml::xml_node<> &root_node,
                                 InitialSimulationState &state) {
  const rapidxml::xml_node<> *node{root_node.first_node("material")};
  if (node == nullptr) {
    std::cerr << "Failed to locate material node." << std::endl;
    return false;
  }

  // Attempt to load the shear modulus
  scalar shear_modulus;
  {
    const rapidxml::xml_attribute<> *const attrib{
        node->first_attribute("shear_modulus")};
    if (!attrib) {
      std::cerr << "Failed to locate shear_modulus attribute for material."
                << std::endl;
      return false;
    }
    std::stringstream ss;
    ss << attrib->value();
    if (!(ss >> shear_modulus)) {
      std::cerr << "Failed to load shear_modulus attribute for material. Must "
                   "provide a scalar."
                << std::endl;
      return false;
    }
    if (shear_modulus < 0.0) {
      std::cerr << "Failed to load shear_modulus attribute for material. Must "
                   "provide a non-negative scalar."
                << std::endl;
      return false;
    }
  }

  // Attempt to load the bulk modulus
  scalar bulk_modulus;
  {
    const rapidxml::xml_attribute<> *const attrib{
        node->first_attribute("bulk_modulus")};
    if (!attrib) {
      std::cerr << "Failed to locate bulk_modulus attribute for material."
                << std::endl;
      return false;
    }
    std::stringstream ss;
    ss << attrib->value();
    if (!(ss >> bulk_modulus)) {
      std::cerr << "Failed to load bulk_modulus attribute for material. Must "
                   "provide a scalar."
                << std::endl;
      return false;
    }
    if (bulk_modulus < 0.0) {
      std::cerr << "Failed to load bulk_modulus attribute for material. Must "
                   "provide a non-negative scalar."
                << std::endl;
      return false;
    }
  }

  // Attempt to load the alpha parameter of the Drucker-Prager yield condition
  scalar Drucker_Prager_alpha;
  {
    const rapidxml::xml_attribute<> *const attrib{
        node->first_attribute("alpha")};
    if (!attrib) {
      std::cerr << "Failed to locate alpha attribute for material."
                << std::endl;
      return false;
    }
    std::stringstream ss;
    ss << attrib->value();
    if (!(ss >> Drucker_Prager_alpha)) {
      std::cerr << "Failed to load alpha attribute for material. Must provide "
                   "a scalar."
                << std::endl;
      return false;
    }
    if (Drucker_Prager_alpha < 0.0) {
      std::cerr << "Failed to load alpha attribute for material. Must provide "
                   "a non-negative scalar."
                << std::endl;
      return false;
    }
  }

  state.bulk_modulus = bulk_modulus;
  state.shear_modulus = shear_modulus;
  state.Drucker_Prager_alpha = Drucker_Prager_alpha;

  return true;
}

static bool loadGridSettings(const rapidxml::xml_node<> &root_node,
                             GridSettings &grid) {
  const rapidxml::xml_node<> *node{root_node.first_node("grid")};
  if (node == nullptr) {
    std::cerr << "Failed to locate grid node." << std::endl;
    return false;
  }

  // Attempt to load the grid minimum
  VectorXs min;
  {
    const rapidxml::xml_attribute<> *const attrib{node->first_attribute("min")};
    if (!attrib) {
      std::cerr << "Failed to locate min attribute for grid node." << std::endl;
      return false;
    }
    if (!StringUtilities::readScalarList(attrib->value(), 2, ' ', min)) {
      std::cerr << "Failed to load min attribute for grid node. Must provide 2 "
                   "scalars."
                << std::endl;
      return false;
    }
  }

  // Attempt to load the grid maximum
  VectorXs max;
  {
    const rapidxml::xml_attribute<> *const attrib{node->first_attribute("max")};
    if (!attrib) {
      std::cerr << "Failed to locate max attribute for grid node." << std::endl;
      return false;
    }
    if (!StringUtilities::readScalarList(attrib->value(), 2, ' ', max)) {
      std::cerr << "Failed to load max attribute for grid node. Must provide 2 "
                   "scalars."
                << std::endl;
      return false;
    }
  }

  // Attempt to load the cell width
  scalar cell_width;
  {
    const rapidxml::xml_attribute<> *const attrib{
        node->first_attribute("cell_width")};
    if (!attrib) {
      std::cerr << "Failed to locate cell_width attribute for grid node."
                << std::endl;
      return false;
    }
    std::stringstream ss;
    ss << attrib->value();
    if (!(ss >> cell_width)) {
      std::cerr << "Failed to load cell_width attribute for grid node. Must "
                   "provide a scalar."
                << std::endl;
      return false;
    }
    if (cell_width <= 0.0) {
      std::cerr << "Failed to load cell_width attribute for grid node. Must "
                   "provide a positive scalar."
                << std::endl;
      return false;
    }
  }

  if ((max.array() <= min.array()).any()) {
    std::cerr << "Failed to load grid node. Grid minimum must be less than "
                 "grid maximum."
              << std::endl;
    return false;
  }

  grid.min = min;
  grid.max = max;
  grid.cell_width = cell_width;

  return true;
}

static bool readSimulationState(const rapidxml::xml_node<> &root_node,
                                InitialSimulationState &sim_state) {
  // Attempt to load any initial rectangles of material
  if (!loadMaterialPointRectangles(root_node, sim_state.rectangle_regions)) {
    return false;
  }

  // Attempt to load the material settings
  if (!loadMaterialSettings(root_node, sim_state)) {
    return false;
  }

  // Attempt to load the grid settings
  if (!loadGridSettings(root_node, sim_state.grid_settings)) {
    return false;
  }

  // Attempt to load the integrator
  if (!loadIntegrator(root_node, sim_state)) {
    std::cerr << "Failed to load integrator in xml scene file." << std::endl;
    return false;
  }

  // Attempt to load static planes
  if (!loadStaticPlanes(root_node, sim_state.static_planes)) {
    std::cerr << "Failed to load static_plane in xml scene file." << std::endl;
    return false;
  }

  // Attempt to load forces
  // Attempt to load a gravity force
  if (!loadGravityForce(root_node, sim_state)) {
    return false;
  }

  // Attempt to load the end time, if present
  sim_state.end_time = SCALAR_INFINITY;
  if (root_node.first_node("end_time") != nullptr) {
    if (!loadEndTime(*root_node.first_node("end_time"), sim_state.end_time)) {
      std::cerr << "Failed to parse end_time node." << std::endl;
      return false;
    }
  }

  return true;
}

bool MPMGrains2DParser::readXMLFile(const std::string &file_name,
                                    InitialSimulationState &initial_state) {
  // Attempt to load the xml document
  std::vector<char> xmlchars;
  rapidxml::xml_document<> doc;
  if (!loadXMLFile(file_name, xmlchars, doc)) {
    return false;
  }

  // Attempt to locate the root node
  if (doc.first_node("mpmgrains2d") == nullptr) {
    std::cerr << "Failed to locate root node mpmgrains2d in xml scene file: "
              << file_name << std::endl;
    return false;
  }
  const rapidxml::xml_node<> &root_node{*doc.first_node("mpmgrains2d")};

  InitialSimulationState sim_state;
  if (!readSimulationState(root_node, sim_state)) {
    return false;
  }

  using std::swap;
  swap(sim_state, initial_state);

  return true;
}