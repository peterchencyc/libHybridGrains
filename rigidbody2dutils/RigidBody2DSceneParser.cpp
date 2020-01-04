#include "RigidBody2DSceneParser.h"

#include <fstream>
#include <iostream>

#include "rapidxml.hpp"

#include "rigidbody2d/AnnulusGeometry.h"
#include "rigidbody2d/BoxGeometry.h"
#include "rigidbody2d/CircleGeometry.h"
#include "rigidbody2d/NearEarthGravityForce.h"
#include "rigidbody2d/PenaltyImpactFrictionMap.h"
#include "rigidbody2d/PlanarPortal.h"
#include "rigidbody2d/RigidBody2DForce.h"
#include "rigidbody2d/RigidBody2DGeometry.h"
#include "rigidbody2d/RigidBody2DState.h"
#include "rigidbody2d/RigidBody2DStaticPlane.h"

#include "scisim/ConstrainedMaps/ImpactMaps/ImpactMap.h"
#include "scisim/Math/Rational.h"

#include "rigidbody2d/RigidBody2DIntegratorSettings.h"

#include "CameraSettings2D.h"

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

static bool loadCameraSettings(const rapidxml::xml_node<> &node,
                               CameraSettings2D &camera_settings) {
  // If the camera is not specified, we are done
  const rapidxml::xml_node<> *camera_node{node.first_node("camera")};
  if (camera_node == nullptr) {
    return true;
  }

  camera_settings.set = true;

  // Attempt to read the center setting
  {
    const rapidxml::xml_attribute<> *center_attrib{
        camera_node->first_attribute("center")};
    if (center_attrib == nullptr) {
      std::cerr << "Failed to locate center attribute for camera node."
                << std::endl;
      return false;
    }
    if (!StringUtilities::readScalarList(center_attrib->value(), 2, ' ',
                                         camera_settings.center)) {
      std::cerr << "Failed to load center attribute for camera node, must "
                   "provide 2 scalars."
                << std::endl;
      return false;
    }
  }

  // Attempt to read the scale setting
  {
    const rapidxml::xml_attribute<> *scale_attrib{
        camera_node->first_attribute("scale")};
    if (scale_attrib == nullptr) {
      std::cerr << "Failed to locate scale attribute for camera node."
                << std::endl;
      return false;
    }
    if (!StringUtilities::extractFromString(scale_attrib->value(),
                                            camera_settings.scale) ||
        camera_settings.scale <= 0.0) {
      std::cerr << "Failed to load scale attribute for camera node, must "
                   "provide a single positive scalar."
                << std::endl;
      return false;
    }
  }

  // Attempt to read the fps setting
  {
    const rapidxml::xml_attribute<> *fps_attrib{
        camera_node->first_attribute("fps")};
    if (fps_attrib == nullptr) {
      std::cerr << "Failed to locate fps attribute for camera node."
                << std::endl;
      return false;
    }
    int fps;
    if (!StringUtilities::extractFromString(fps_attrib->value(), fps) ||
        fps <= 0) {
      std::cerr << "Failed to parse fps attribute for camera node, must "
                   "provide a non-negative integer."
                << std::endl;
      return false;
    }
    camera_settings.fps = unsigned(fps);
  }

  // Attempt to read the render_at_fps setting
  {
    const rapidxml::xml_attribute<> *render_at_fps_attrib{
        camera_node->first_attribute("render_at_fps")};
    if (render_at_fps_attrib == nullptr) {
      std::cerr << "Failed to locate render_at_fps attribute for camera node."
                << std::endl;
      return false;
    }
    if (!StringUtilities::extractFromString(render_at_fps_attrib->value(),
                                            camera_settings.render_at_fps)) {
      std::cerr << "Failed to parse render_at_fps attribute for camera node, "
                   "must provide a boolean."
                << std::endl;
      return false;
    }
  }

  // Attempt to read the locked setting
  {
    const rapidxml::xml_attribute<> *locked_attrib{
        camera_node->first_attribute("locked")};
    if (locked_attrib == nullptr) {
      std::cerr << "Failed to locate locked attribute for camera node."
                << std::endl;
      return false;
    }
    if (!StringUtilities::extractFromString(locked_attrib->value(),
                                            camera_settings.locked)) {
      std::cerr << "Failed to parse locked attribute for camera node, must "
                   "provide a boolean."
                << std::endl;
      return false;
    }
  }

  return true;
}

static bool loadEndTime(const rapidxml::xml_node<> &node, scalar &end_time) {
  // If the end time is not specified, set it to infinity
  const rapidxml::xml_node<> *end_time_node{node.first_node("end_time")};
  if (end_time_node == nullptr) {
    end_time = SCALAR_INFINITY;
    return true;
  }

  // Attempt to parse the time setting
  const rapidxml::xml_attribute<> *t_attrib{
      end_time_node->first_attribute("t")};
  if (t_attrib == nullptr) {
    std::cerr << "Failed to locate t attribute for end_time node." << std::endl;
    return false;
  }
  if (!StringUtilities::extractFromString(t_attrib->value(), end_time) ||
      end_time <= 0.0) {
    std::cerr << "Failed to parse t attribute for end_time. Value must be a "
                 "positive scalar."
              << std::endl;
    return false;
  }

  return true;
}

static bool loadScriptingSetup(const rapidxml::xml_node<> &node,
                               std::string &scripting_callback) {
  assert(scripting_callback.empty());

  const rapidxml::xml_node<> *scripting_node{node.first_node("scripting")};
  if (!scripting_node) {
    return true;
  }

  const rapidxml::xml_attribute<> *name_node{
      scripting_node->first_attribute("callback")};
  if (name_node) {
    scripting_callback = name_node->value();
  } else {
    return false;
  }

  return true;
}

static bool loadStaticPlanes(const rapidxml::xml_node<> &node,
                             std::vector<RigidBody2DStaticPlane> &planes) {
  for (rapidxml::xml_node<> *nd = node.first_node("static_plane"); nd;
       nd = nd->next_sibling("static_plane")) {
    // Attempt to read the point on the plane
    Vector2s x;
    {
      const rapidxml::xml_attribute<> *x_attrib{nd->first_attribute("x")};
      if (x_attrib == nullptr) {
        std::cerr << "Failed to locate x attribute for static_plane node."
                  << std::endl;
        return false;
      }
      if (!StringUtilities::readScalarList(x_attrib->value(), 2, ' ', x)) {
        std::cerr << "Failed to load x attribute for static_plane node, must "
                     "provide 2 scalars."
                  << std::endl;
        return false;
      }
    }

    // Attempt to read the plane's normal
    Vector2s n;
    {
      const rapidxml::xml_attribute<> *n_attrib{nd->first_attribute("n")};
      if (n_attrib == nullptr) {
        std::cerr << "Failed to locate n attribute for static_plane node."
                  << std::endl;
        return false;
      }
      if (!StringUtilities::readScalarList(n_attrib->value(), 2, ' ', n)) {
        std::cerr << "Failed to load n attribute for static_plane node, must "
                     "provide 2 scalars."
                  << std::endl;
        return false;
      }
      if (n.norm() == 0.0) {
        std::cerr << "Failed to load n attribute for static_plane node, must "
                     "provide a nonzero vector."
                  << std::endl;
        return false;
      }
    }
    planes.emplace_back(x, n.normalized());
  }

  return true;
}

// TODO: minus ones here can underflow
static bool loadPlanarPortals(const rapidxml::xml_node<> &node,
                              std::vector<RigidBody2DStaticPlane> &planes,
                              std::vector<PlanarPortal> &planar_portals) {
  if (!(node.first_node("planar_portal") ||
        node.first_node("lees_edwards_portal"))) {
    return true;
  }

  // If we have a portal we must have at least one plane
  if (planes.size() < 2) {
    std::cerr << "Error, must provide at least two planes before instantiating "
                 "a planar portal."
              << std::endl;
    return false;
  }

  // Pairs of planes to turn into portals
  std::vector<std::pair<unsigned, unsigned>> plane_pairs;
  std::vector<scalar> plane_tangent_velocities;
  std::vector<scalar> plane_bounds;

  // Load planes without kinematic velocities
  for (rapidxml::xml_node<> *nd = node.first_node("planar_portal"); nd;
       nd = nd->next_sibling("planar_portal")) {
    // Read the first plane index
    const rapidxml::xml_attribute<> *const attrib_a{
        nd->first_attribute("planeA")};
    if (!attrib_a) {
      std::cerr << "Failed to locate planeA attribute for planar_portal node."
                << std::endl;
      return false;
    }
    unsigned index_a;
    if (!StringUtilities::extractFromString(attrib_a->value(), index_a)) {
      std::cerr << "Failed to parse planeA attribute for planar_portal node "
                   "with value "
                << attrib_a->value()
                << ". Attribute must be an unsigned integer." << std::endl;
      return false;
    }
    if (index_a >= planes.size()) {
      std::cerr << "Failed to parse planeA attribute for planar_portal node "
                   "with value "
                << attrib_a->value()
                << ". Attribute must be an index of a plane between " << 0
                << " and " << planes.size() - 1 << std::endl;
      return false;
    }

    // Read the second plane index
    const rapidxml::xml_attribute<> *const attrib_b{
        nd->first_attribute("planeB")};
    if (!attrib_b) {
      std::cerr << "Failed to locate planeB attribute for planar_portal node."
                << std::endl;
      return false;
    }
    unsigned index_b;
    if (!StringUtilities::extractFromString(attrib_b->value(), index_b)) {
      std::cerr << "Failed to parse planeB attribute for planar_portal node "
                   "with value "
                << attrib_b->value()
                << ". Attribute must be an unsigned integer." << std::endl;
      return false;
    }
    if (index_b >= planes.size()) {
      std::cerr << "Failed to parse planeB attribute for planar_portal node "
                   "with value "
                << attrib_b->value()
                << ". Attribute must be an index of a plane between " << 0
                << " and " << planes.size() - 1 << std::endl;
      return false;
    }

    // Indices for this portal can't repeat
    if (index_a == index_b) {
      std::cerr << "Failed to parse planeB attribute for planar_portal node "
                   "with value "
                << attrib_b->value()
                << ". Value is a repeat of attribute planeA." << std::endl;
      return false;
    }

    // Indices can not repeat previous portals
    for (std::vector<std::pair<unsigned, unsigned>>::size_type i = 0;
         i < plane_pairs.size(); ++i) {
      if (index_a == plane_pairs[i].first || index_a == plane_pairs[i].second) {
        std::cerr << "Failed to parse planeA attribute for planar_portal node "
                     "with value "
                  << attrib_a->value()
                  << ". Plane index is used by an existing portal."
                  << std::endl;
        return false;
      }
      if (index_b == plane_pairs[i].first || index_b == plane_pairs[i].second) {
        std::cerr << "Failed to parse planeB attribute for planar_portal node "
                     "with value "
                  << attrib_b->value()
                  << ". Plane index is used by an existing portal."
                  << std::endl;
        return false;
      }
    }

    plane_pairs.emplace_back(index_a, index_b);
    plane_tangent_velocities.emplace_back(0.0);
    plane_bounds.emplace_back(0.0);
  }

  // Load planes with kinematic velocities
  for (rapidxml::xml_node<> *nd = node.first_node("lees_edwards_portal"); nd;
       nd = nd->next_sibling("lees_edwards_portal")) {
    // Read the first plane index
    const rapidxml::xml_attribute<> *const attrib_a =
        nd->first_attribute("planeA");
    if (!attrib_a) {
      std::cerr
          << "Failed to locate planeA attribute for lees_edwards_portal node."
          << std::endl;
      return false;
    }
    unsigned index_a;
    if (!StringUtilities::extractFromString(attrib_a->value(), index_a)) {
      std::cerr << "Failed to parse planeA attribute for lees_edwards_portal "
                   "node with value "
                << attrib_a->value()
                << ". Attribute must be an unsigned integer." << std::endl;
      return false;
    }
    if (index_a >= planes.size()) {
      std::cerr << "Failed to parse planeA attribute for lees_edwards_portal "
                   "node with value "
                << attrib_a->value()
                << ". Attribute must be an index of a plane between " << 0
                << " and " << planes.size() - 1 << std::endl;
      return false;
    }

    // Read the second plane index
    const rapidxml::xml_attribute<> *const attrib_b{
        nd->first_attribute("planeB")};
    if (!attrib_b) {
      std::cerr
          << "Failed to locate planeB attribute for lees_edwards_portal node."
          << std::endl;
      return false;
    }
    unsigned index_b;
    if (!StringUtilities::extractFromString(attrib_b->value(), index_b)) {
      std::cerr << "Failed to parse planeB attribute for lees_edwards_portal "
                   "node with value "
                << attrib_b->value()
                << ". Attribute must be an unsigned integer." << std::endl;
      return false;
    }
    if (index_b >= planes.size()) {
      std::cerr << "Failed to parse planeB attribute for lees_edwards_portal "
                   "node with value "
                << attrib_b->value()
                << ". Attribute must be an index of a plane between " << 0
                << " and " << planes.size() - 1 << std::endl;
      return false;
    }

    // Indices for this portal can't repeat
    if (index_a == index_b) {
      std::cerr << "Failed to parse planeB attribute for lees_edwards_portal "
                   "node with value "
                << attrib_b->value()
                << ". Value is a repeat of attribute planeA." << std::endl;
      return false;
    }

    // Indices can not repeat previous portals
    for (std::vector<std::pair<unsigned, unsigned>>::size_type i = 0;
         i < plane_pairs.size(); ++i) {
      if (index_a == plane_pairs[i].first || index_a == plane_pairs[i].second) {
        std::cerr << "Failed to parse planeA attribute for lees_edwards_portal "
                     "node with value "
                  << attrib_a->value()
                  << ". Plane index is used by an existing portal."
                  << std::endl;
        return false;
      }
      if (index_b == plane_pairs[i].first || index_b == plane_pairs[i].second) {
        std::cerr << "Failed to parse planeB attribute for lees_edwards_portal "
                     "node with value "
                  << attrib_b->value()
                  << ". Plane index is used by an existing portal."
                  << std::endl;
        return false;
      }
    }

    plane_pairs.emplace_back(index_a, index_b);

    // Load the velocity of the portal
    scalar v;
    {
      if (nd->first_attribute("v") == nullptr) {
        std::cerr << "Could not locate v attribue for lees_edwards_portal"
                  << std::endl;
        return false;
      }
      const rapidxml::xml_attribute<> &v_attrib{*nd->first_attribute("v")};
      if (!StringUtilities::extractFromString(std::string{v_attrib.value()},
                                              v)) {
        std::cerr << "Could not load v attribue for lees_edwards_portal, value "
                     "must be a scalar"
                  << std::endl;
        return false;
      }
    }

    // Load the bounds on portal a's translation
    scalar bounds;
    {
      if (nd->first_attribute("bounds") == nullptr) {
        std::cerr << "Could not locate bounds attribue for lees_edwards_portal"
                  << std::endl;
        return false;
      }
      const rapidxml::xml_attribute<> &bounds_attrib{
          *nd->first_attribute("bounds")};
      if (!StringUtilities::extractFromString(
              std::string{bounds_attrib.value()}, bounds)) {
        std::cerr << "Could not load bounds attribue for lees_edwards_portal, "
                     "value must be a scalar"
                  << std::endl;
        return false;
      }
      if (bounds < 0) {
        std::cerr << "Failed to load bounds attribute for lees_edwards_portal, "
                     "value must be positive"
                  << std::endl;
        return false;
      }
    }

    plane_tangent_velocities.emplace_back(v);
    plane_bounds.emplace_back(bounds);
  }

  assert(plane_pairs.size() == plane_tangent_velocities.size());
  assert(plane_pairs.size() == plane_bounds.size());
  for (std::vector<std::pair<unsigned, unsigned>>::size_type i = 0;
       i < plane_pairs.size(); ++i) {
    planar_portals.emplace_back(planes[plane_pairs[i].first],
                                planes[plane_pairs[i].second],
                                plane_tangent_velocities[i], plane_bounds[i]);
  }

  // TODO: This could get slow if there are a ton of portals, but probably not
  // too big of a deal for now
  {
    std::vector<unsigned> indices;
    for (std::vector<std::pair<unsigned, unsigned>>::size_type i = 0;
         i < plane_pairs.size(); ++i) {
      indices.emplace_back(plane_pairs[i].first);
      indices.emplace_back(plane_pairs[i].second);
    }
    std::sort(indices.begin(), indices.end());
    for (std::vector<unsigned>::size_type i = indices.size(); i-- > 0;) {
      planes.erase(planes.begin() + indices[i]);
    }
  }

  return true;
}

static bool loadStaticDrums(const rapidxml::xml_node<> &node,
                            std::vector<RigidBody2DStaticDrum> &drums) {
  for (rapidxml::xml_node<> *nd = node.first_node("static_drum"); nd;
       nd = nd->next_sibling("static_drum")) {
    // Attempt to read the drum's center
    Vector2s x;
    {
      const rapidxml::xml_attribute<> *const attrib{nd->first_attribute("x")};
      if (attrib == nullptr) {
        std::cerr << "Failed to locate x attribute for static_drum node."
                  << std::endl;
        return false;
      }
      if (!StringUtilities::readScalarList(attrib->value(), 2, ' ', x)) {
        std::cerr << "Failed to load x attribute for static_drum node, must "
                     "provide 2 scalars."
                  << std::endl;
        return false;
      }
    }

    // Attempt to read the drum's radius
    scalar radius;
    {
      const rapidxml::xml_attribute<> *const attrib{nd->first_attribute("r")};
      if (attrib == nullptr) {
        std::cerr << "Failed to locate r attribute for static_drum node."
                  << std::endl;
        return false;
      }
      if (!StringUtilities::extractFromString(std::string{attrib->value()},
                                              radius) ||
          radius <= 0.0) {
        std::cerr << "Failed to load r attribute for static_drum node, must "
                     "provide a positive scalar."
                  << std::endl;
        return false;
      }
    }

    drums.emplace_back(x, radius);
  }

  return true;
}

static bool
loadGravityForce(const rapidxml::xml_node<> &node,
                 std::vector<std::unique_ptr<RigidBody2DForce>> &forces) {
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
    forces.emplace_back(new NearEarthGravityForce{f});
  }

  return true;
}

static bool
loadGeometry(const rapidxml::xml_node<> &node,
             std::vector<std::unique_ptr<RigidBody2DGeometry>> &geometry) {
  geometry.clear();

  for (rapidxml::xml_node<> *nd = node.first_node("geometry"); nd;
       nd = nd->next_sibling("geometry")) {
    // Load the type of the geometry
    std::string geom_type;
    {
      const rapidxml::xml_attribute<> *const attrib{
          nd->first_attribute("type")};
      if (attrib == nullptr) {
        std::cerr << "Failed to locate type attribute for geometry node."
                  << std::endl;
        return false;
      }
      geom_type = attrib->value();
      // Parse the remaining arguments based on the geometry type
      if (geom_type == "circle") {
        // Read the radius
        const rapidxml::xml_attribute<> *const r_attrib{
            nd->first_attribute("r")};
        if (r_attrib == nullptr) {
          std::cerr << "Failed to locate r attribute for circle geometry node."
                    << std::endl;
          return false;
        }
        scalar r;
        if (!StringUtilities::extractFromString(r_attrib->value(), r) ||
            r <= 0.0) {
          std::cerr << "Failed to read r attribute for circle geometry, must "
                       "provide a positive scalar."
                    << std::endl;
          return false;
        }
        geometry.emplace_back(new CircleGeometry{r});
      } else if (geom_type == "box") {
        // Read the half-widths
        const rapidxml::xml_attribute<> *const r_attrib{
            nd->first_attribute("r")};
        if (r_attrib == nullptr) {
          std::cerr << "Failed to locate r attribute for box geometry node."
                    << std::endl;
          return false;
        }
        Vector2s r;
        if (!StringUtilities::readScalarList(r_attrib->value(), 2, ' ', r) ||
            (r.array() <= 0.0).any()) {
          std::cerr << "Failed to read r attribute for box geometry, must "
                       "provide two positive scalars."
                    << std::endl;
          return false;
        }
        geometry.emplace_back(new BoxGeometry{r});
      } else if (geom_type == "annulus") {
        const rapidxml::xml_attribute<> *const r0_attrib{
            nd->first_attribute("r0")};
        if (r0_attrib == nullptr) {
          std::cerr << "Failed to locate r0 attribute for circle geometry node."
                    << std::endl;
          return false;
        }
        const rapidxml::xml_attribute<> *const r1_attrib{
            nd->first_attribute("r1")};
        if (r1_attrib == nullptr) {
          std::cerr << "Failed to locate r1 attribute for circle geometry node."
                    << std::endl;
          return false;
        }
        scalar r0;
        if (!StringUtilities::extractFromString(r0_attrib->value(), r0) ||
            r0 < 0.0) {
          std::cerr << "Failed to read r0 attribute for circle geometry, must "
                       "provide a positive scalar."
                    << std::endl;
          return false;
        }
        scalar r1;
        if (!StringUtilities::extractFromString(r1_attrib->value(), r1) ||
            r1 <= 0.0) {
          std::cerr << "Failed to read r1 attribute for circle geometry, must "
                       "provide a positive scalar."
                    << std::endl;
          return false;
        }
        if (r0 > r1) {
          std::cerr << "Failed to load annulus geometry, r0 must be less than "
                       "or equal to r1."
                    << std::endl;
          return false;
        }
        geometry.emplace_back(new AnnulusGeometry{r0, r1});
      } else {
        std::cerr << "Invalid type specified for geometry node of " << geom_type
                  << "." << std::endl;
        return false;
      }
    }
  }

  return true;
}

static bool
loadBodies(const rapidxml::xml_node<> &node,
           const std::vector<std::unique_ptr<RigidBody2DGeometry>> &geometry,
           VectorXs &q, VectorXs &v, VectorXs &m, VectorXu &indices,
           std::vector<bool> &fixed) {
  std::vector<Vector2s> xs;
  std::vector<scalar> thetas;
  std::vector<Vector2s> vs;
  std::vector<scalar> omegas;
  std::vector<scalar> densities;
  std::vector<unsigned> geometry_indices;
  assert(fixed.empty());

  for (rapidxml::xml_node<> *nd = node.first_node("rigid_body"); nd;
       nd = nd->next_sibling("rigid_body")) {
    // Load the center of mass' position
    {
      const rapidxml::xml_attribute<> *const x_attrib{nd->first_attribute("x")};
      if (x_attrib == nullptr) {
        std::cerr << "Failed to locate x attribute for rigid_body node."
                  << std::endl;
        return false;
      }
      Vector2s x;
      if (!StringUtilities::readScalarList(x_attrib->value(), 2, ' ', x)) {
        std::cerr << "Failed to load x attribute for rigid_body node, must "
                     "provide two scalars."
                  << std::endl;
        return false;
      }
      xs.emplace_back(x);
    }

    // Load the rotation about the center of mass
    {
      const rapidxml::xml_attribute<> *const theta_attrib{
          nd->first_attribute("theta")};
      if (theta_attrib == nullptr) {
        std::cerr << "Failed to locate theta attribute for rigid_body node."
                  << std::endl;
        return false;
      }
      scalar theta;
      if (!StringUtilities::extractFromString(theta_attrib->value(), theta)) {
        std::cerr << "Failed to load theta attribute for rigid_body node, must "
                     "provide a single scalar."
                  << std::endl;
        return false;
      }
      thetas.emplace_back(theta);
    }

    // Load the center of mass' velocity
    {
      const rapidxml::xml_attribute<> *const v_attrib{nd->first_attribute("v")};
      if (v_attrib == nullptr) {
        std::cerr << "Failed to locate v attribute for rigid_body node."
                  << std::endl;
        return false;
      }
      Vector2s v_body;
      if (!StringUtilities::readScalarList(v_attrib->value(), 2, ' ', v_body)) {
        std::cerr << "Failed to load v attribute for rigid_body node, must "
                     "provide two scalars."
                  << std::endl;
        return false;
      }
      vs.emplace_back(v_body);
    }

    // Load the angular velocity
    {
      const rapidxml::xml_attribute<> *const omega_attrib{
          nd->first_attribute("omega")};
      if (omega_attrib == nullptr) {
        std::cerr << "Failed to locate omega attribute for rigid_body node."
                  << std::endl;
        return false;
      }
      scalar omega;
      if (!StringUtilities::extractFromString(omega_attrib->value(), omega)) {
        std::cerr << "Failed to load omega attribute for rigid_body node, must "
                     "provide a single scalar."
                  << std::endl;
        return false;
      }
      omegas.emplace_back(omega);
    }

    // Load the density
    {
      const rapidxml::xml_attribute<> *const rho_attrib{
          nd->first_attribute("rho")};
      if (rho_attrib == nullptr) {
        std::cerr << "Failed to locate rho attribute for rigid_body node."
                  << std::endl;
        return false;
      }
      scalar rho;
      if (!StringUtilities::extractFromString(rho_attrib->value(), rho) ||
          rho <= 0.0) {
        std::cerr << "Failed to load rho attribute for rigid_body node, must "
                     "provide a single positive scalar."
                  << std::endl;
        return false;
      }
      densities.emplace_back(rho);
    }

    // Load the index of this body's geometry
    {
      const rapidxml::xml_attribute<> *const geo_idx_attrib{
          nd->first_attribute("geo_idx")};
      if (geo_idx_attrib == nullptr) {
        std::cerr << "Failed to locate geo_idx attribute for rigid_body node."
                  << std::endl;
        return false;
      }
      int geometry_index;
      if (!StringUtilities::extractFromString(geo_idx_attrib->value(),
                                              geometry_index) ||
          geometry_index >= int(geometry.size())) {
        std::cerr << "Failed to load geo_idx attribute for rigid_body node, "
                     "must provide an unsigned integer less than the number of "
                     "geometry instances."
                  << std::endl;
        return false;
      }
      if (geometry_index < 0) {
        geometry_index = geometry.size() + geometry_index;
      }
      assert(geometry_index >= 0);
      assert(geometry_index < int(geometry.size()));
      geometry_indices.emplace_back(unsigned(geometry_index));
    }

    // Load the optional fixed attribute
    {
      const rapidxml::xml_attribute<> *const fixed_attrib{
          nd->first_attribute("fixed")};
      if (fixed_attrib == nullptr) {
        fixed.push_back(false);
      } else {
        bool fixed_val;
        if (!StringUtilities::extractFromString(fixed_attrib->value(),
                                                fixed_val)) {
          std::cerr << "Failed to load fixed attribute for rigid_body node, "
                       "fixed must be a boolean."
                    << std::endl;
          return false;
        }
        fixed.push_back(fixed_val);
      }
    }
  }

  const unsigned nbodies{static_cast<unsigned>(xs.size())};

  // Build the generalized configuration
  q.resize(3 * nbodies);
  for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
    q.segment<2>(3 * bdy_idx) = xs[bdy_idx];
    q(3 * bdy_idx + 2) = thetas[bdy_idx];
  }

  // Build the generalized velocity
  v.resize(3 * nbodies);
  for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
    v.segment<2>(3 * bdy_idx) = vs[bdy_idx];
    v(3 * bdy_idx + 2) = omegas[bdy_idx];
  }

  // Copy the geometry indices to the output
  indices = Eigen::Map<VectorXu>(geometry_indices.data(), nbodies);

  // Build the diagonal of the mass matrix
  m.resize(3 * nbodies);
  for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
    scalar mass;
    scalar inertia;
    geometry[indices(bdy_idx)]->computeMassAndInertia(densities[bdy_idx], mass,
                                                      inertia);
    m.segment<2>(3 * bdy_idx).setConstant(mass);
    m(3 * bdy_idx + 2) = inertia;
  }

  return true;
}

static bool loadCollisionDetection(const rapidxml::xml_node<> &node,
                                   scalar &spatial_grid_scale,
                                   bool &reduce_bandwidth) {
  const rapidxml::xml_node<> *const cd_node{
      node.first_node("collision_detection")};
  if (cd_node == nullptr) {
    std::cerr << "Failed to locate collision_detection node." << std::endl;
    return false;
  }

  {
    const rapidxml::xml_attribute<> *const broadphase_attrib{
        cd_node->first_attribute("broadphase")};
    if (broadphase_attrib == nullptr) {
      std::cerr << "Failed to locate broadphase attribute for "
                   "collision_detection node."
                << std::endl;
      return false;
    }
    if (broadphase_attrib->value() != std::string{"spatial_grid"}) {
      std::cerr
          << "Failed to load broadphase attribute for collision_detection "
             "node. Invalid value specified. Valid options are: spatial_grid"
          << std::endl;
      return false;
    }
  }

  {
    const rapidxml::xml_attribute<> *const scale_attrib{
        cd_node->first_attribute("cell_scale")};
    if (scale_attrib == nullptr) {
      std::cerr << "Failed to locate cell_scale attribute for "
                   "collision_detection node."
                << std::endl;
      return false;
    }
    const bool parsed{StringUtilities::extractFromString(scale_attrib->value(),
                                                         spatial_grid_scale)};
    if (!parsed || spatial_grid_scale <= 0.0) {
      std::cerr << "Failed to parse cell_scale attribute for "
                   "collision_detection node, must provide a positive scalar."
                << std::endl;
      return false;
    }
  }

  {
    const rapidxml::xml_attribute<> *bandwidth_attrib{
        cd_node->first_attribute("bandwidth_reduction")};
    if (bandwidth_attrib == nullptr) {
      std::cerr << "Failed to locate bandwidth_reduction attribute for "
                   "collision_detection node."
                << std::endl;
      return false;
    }
    if (bandwidth_attrib->value() == std::string{"none"}) {
      reduce_bandwidth = false;
    } else if (bandwidth_attrib->value() == std::string{"cuthill_mckee"}) {
      reduce_bandwidth = true;
    } else {
      std::cerr << "Failed to load bandwidth_reduction attribute for "
                   "collision_detection node. Invalid value specified. Valid "
                   "options are: none, cuthill_mckee"
                << std::endl;
      return false;
    }
  }

  return true;
}

static bool loadPenaltyImpactFrictionMap(
    const rapidxml::xml_node<> &node, std::string &dt_string,
    Rational<std::intmax_t> &dt, std::unique_ptr<ImpactFrictionMap> &if_map) {
  // Attempt to load the timestep
  Rational<std::intmax_t> new_dt;
  std::string new_dt_string;
  {
    if (!node.first_attribute("dt")) {
      std::cerr
          << "Error, failed to locate dt attribute for penalty_impact_map."
          << std::endl;
      return false;
    }
    if (!RationalTools::extractFromString(node.first_attribute("dt")->value(),
                                          new_dt)) {
      std::cerr << "Error, failed to load dt's value for penalty_impact_map. "
                   "Must provide a number."
                << std::endl;
      return false;
    }
    if (!new_dt.positive()) {
      std::cerr << "Error, failed to load dt's value for penalty_impact_map. "
                   "Must provide a positive number."
                << std::endl;
      return false;
    }
    new_dt_string = node.first_attribute("dt")->value();
  }

  scalar kn;
  {
    if (!node.first_attribute("kn")) {
      std::cerr
          << "Error, failed to locate kn attribute for penalty_impact_map."
          << std::endl;
      return false;
    }
    if (!StringUtilities::extractFromString(node.first_attribute("kn")->value(),
                                            kn)) {
      std::cerr << "Error, failed to load kn's value for penalty_impact_map. "
                   "Must provide a number."
                << std::endl;
      return false;
    }
    if (kn < 0.0) {
      std::cerr << "Error, failed to load kn's value for penalty_impact_map. "
                   "Must provide a positive number."
                << std::endl;
      return false;
    }
  }

  scalar kt;
  {
    if (!node.first_attribute("kt")) {
      std::cerr
          << "Error, failed to locate kt attribute for penalty_impact_map."
          << std::endl;
      return false;
    }
    if (!StringUtilities::extractFromString(node.first_attribute("kt")->value(),
                                            kt)) {
      std::cerr << "Error, failed to load kt's value for penalty_impact_map. "
                   "Must provide a number."
                << std::endl;
      return false;
    }
    if (kt < 0.0) {
      std::cerr << "Error, failed to load kt's value for penalty_impact_map. "
                   "Must provide a positive number."
                << std::endl;
      return false;
    }
  }

  scalar gamman;
  {
    if (!node.first_attribute("gamman")) {
      std::cerr
          << "Error, failed to locate gamman attribute for penalty_impact_map."
          << std::endl;
      return false;
    }
    if (!StringUtilities::extractFromString(
            node.first_attribute("gamman")->value(), gamman)) {
      std::cerr << "Error, failed to load gamman's value for "
                   "penalty_impact_map. Must provide a number."
                << std::endl;
      return false;
    }
    if (gamman < 0.0) {
      std::cerr << "Error, failed to load gamman's value for "
                   "penalty_impact_map. Must provide a positive number."
                << std::endl;
      return false;
    }
  }

  scalar gammat;
  {
    if (!node.first_attribute("gammat")) {
      std::cerr
          << "Error, failed to locate gammat attribute for penalty_impact_map."
          << std::endl;
      return false;
    }
    if (!StringUtilities::extractFromString(
            node.first_attribute("gammat")->value(), gammat)) {
      std::cerr << "Error, failed to load gammat's value for "
                   "penalty_impact_map. Must provide a number."
                << std::endl;
      return false;
    }
    if (gammat < 0.0) {
      std::cerr << "Error, failed to load gammat's value for "
                   "penalty_impact_map. Must provide a positive number."
                << std::endl;
      return false;
    }
  }

  scalar mu;
  {
    if (!node.first_attribute("mu")) {
      std::cerr
          << "Error, failed to locate mu attribute for penalty_impact_map."
          << std::endl;
      return false;
    }
    if (!StringUtilities::extractFromString(node.first_attribute("mu")->value(),
                                            mu)) {
      std::cerr << "Error, failed to load mu's value for penalty_impact_map. "
                   "Must provide a number."
                << std::endl;
      return false;
    }
    if (mu < 0.0) {
      std::cerr << "Error, failed to load mu's value for penalty_impact_map. "
                   "Must provide a positive number."
                << std::endl;
      return false;
    }
  }

  if_map.reset(new PenaltyImpactFrictionMap(kn, kt, gamman, gammat, mu));
  std::swap(new_dt, dt);
  std::swap(new_dt_string, dt_string);

  return true;
}

static bool
basicParse(const rapidxml::xml_document<> &doc, const std::string &file_name,
           std::string &scripting_callback, RigidBody2DState &sim_state,
           std::unique_ptr<UnconstrainedMap> &unconstrained_map,
           std::string &dt_string, Rational<std::intmax_t> &dt,
           scalar &end_time, std::unique_ptr<ImpactOperator> &impact_operator,
           std::unique_ptr<ImpactMap> &impact_map, scalar &CoR,
           std::unique_ptr<FrictionSolver> &friction_solver, scalar &mu,
           std::unique_ptr<ImpactFrictionMap> &if_map,
           scalar &spatial_grid_scale, bool &reduce_bandwidth,
           CameraSettings2D &camera_settings) {
  // Attempt to locate the root node
  if (doc.first_node("rigidbody2d_scene") == nullptr) {
    std::cerr
        << "Failed to locate root node rigidbody2d_scene in xml scene file: "
        << file_name << std::endl;
    return false;
  }
  const rapidxml::xml_node<> &root_node{*doc.first_node("rigidbody2d_scene")};

  // Attempt to load an optional scripting callback
  if (!loadScriptingSetup(root_node, scripting_callback)) {
    return false;
  }

  // Attempt to load the optional end time, if present
  if (!loadEndTime(root_node, end_time)) {
    return false;
  }

  // Attempt to load the optional camera settings, if present
  if (!loadCameraSettings(root_node, camera_settings)) {
    return false;
  }

  // Attempt to load an impact operator
  impact_operator.reset(nullptr);
  impact_map.reset(nullptr);
  CoR = SCALAR_NAN;

  friction_solver.reset(nullptr);
  mu = SCALAR_NAN;
  if_map.reset(nullptr);

  if (root_node.first_node("penalty_impact_map") != nullptr) {
    if (impact_operator != nullptr) {
      std::cerr
          << "Error loading penalty_impact_map, impact_operator specified: "
          << impact_operator->name() << std::endl;
      return false;
    }
    if (impact_map != nullptr) {
      std::cerr << "Error loading penalty_impact_map, impact_map specified."
                << std::endl;
      return false;
    }
    if (friction_solver != nullptr) {
      std::cerr
          << "Error loading penalty_impact_map, friction_solver specified: "
          << friction_solver->name() << std::endl;
      return false;
    }
    if (if_map != nullptr) {
      std::cerr << "Error loading penalty_impact_map, if_map of type "
                << if_map->name() << " already specified" << std::endl;
      return false;
    }
    if (unconstrained_map != nullptr) {
      std::cerr
          << "Error loading penalty_impact_map, unconstrained_map of type "
          << unconstrained_map->name() << " already specified" << std::endl;
      return false;
    }
    assert(impact_operator == nullptr);
    assert(impact_map == nullptr);
    assert(friction_solver == nullptr);
    assert(if_map == nullptr);
    assert(unconstrained_map == nullptr);

    if (!loadPenaltyImpactFrictionMap(
            *root_node.first_node("penalty_impact_map"), dt_string, dt,
            if_map)) {
      std::cerr << "Failed to load penalty_impact_map in xml scene file: "
                << file_name << std::endl;
      return false;
    }

    // unconstrained_map.reset( new NullUnconstrainedMap );
    // friction_solver.reset( new NullFrictionSolver );

    assert(impact_operator == nullptr);
    assert(impact_map == nullptr);
    assert(friction_solver == nullptr);
    assert(if_map != nullptr);
    assert(unconstrained_map == nullptr);
  }

  if (!loadCollisionDetection(root_node, spatial_grid_scale,
                              reduce_bandwidth)) {
    std::cerr << "Failed to load collision_detection node in xml scene file: "
              << file_name << std::endl;
    return false;
  }

  // Load forces
  std::vector<std::unique_ptr<RigidBody2DForce>> forces;
  // Attempt to load a gravity force
  if (!loadGravityForce(root_node, forces)) {
    return false;
  }

  // Attempt to load any user-provided static drums
  std::vector<RigidBody2DStaticPlane> planes;
  if (!loadStaticPlanes(root_node, planes)) {
    return false;
  }

  // Attempt to load planar portals
  std::vector<PlanarPortal> planar_portals;
  if (!loadPlanarPortals(root_node, planes, planar_portals)) {
    return false;
  }

  std::vector<RigidBody2DStaticDrum> drums;
  if (!loadStaticDrums(root_node, drums)) {
    return false;
  }

  // Load geometry to attatch to bodies
  std::vector<std::unique_ptr<RigidBody2DGeometry>> geometry;
  if (!loadGeometry(root_node, geometry)) {
    return false;
  }

  // TODO: Load state directly into RigidBody2DState
  // Load the state of the bodies
  VectorXs q;
  VectorXs v;
  VectorXs m;
  VectorXu indices;
  std::vector<bool> fixed;
  if (!loadBodies(root_node, geometry, q, v, m, indices, fixed)) {
    return false;
  }

  sim_state = RigidBody2DState{
      q, v, m, fixed, indices, geometry, forces, planes, planar_portals, drums};

  return true;
}

bool RigidBody2DSceneParser::parseXMLSceneFile(
    const std::string &file_name, std::string &scripting_callback,
    RigidBody2DState &sim_state,
    std::unique_ptr<UnconstrainedMap> &unconstrained_map,
    std::string &dt_string, Rational<std::intmax_t> &dt, scalar &end_time,
    std::unique_ptr<ImpactOperator> &impact_operator,
    std::unique_ptr<ImpactMap> &impact_map, scalar &CoR,
    std::unique_ptr<FrictionSolver> &friction_solver, scalar &mu,
    std::unique_ptr<ImpactFrictionMap> &if_map, scalar &spatial_grid_scale,
    bool &reduce_bandwidth, CameraSettings2D &camera_settings) {
  // Attempt to load the xml document
  std::vector<char> xmlchars;
  rapidxml::xml_document<> doc;
  if (!loadXMLFile(file_name, xmlchars, doc)) {
    return false;
  }
  return basicParse(doc, file_name, scripting_callback, sim_state,
                    unconstrained_map, dt_string, dt, end_time, impact_operator,
                    impact_map, CoR, friction_solver, mu, if_map,
                    spatial_grid_scale, reduce_bandwidth, camera_settings);
}

bool RigidBody2DSceneParser::parseXMLSceneFile(
    const std::string &file_name, std::string &scripting_callback,
    RigidBody2DState &sim_state,
    RigidBody2DIntegratorSettings &integrator_settings, std::string &dt_string,
    CameraSettings2D &rendering_state) {
  return parseXMLSceneFile(
      file_name, scripting_callback, sim_state,
      integrator_settings.unconstrained_map, dt_string, integrator_settings.dt,
      integrator_settings.end_time, integrator_settings.impact_operator,
      integrator_settings.impact_map, integrator_settings.CoR,
      integrator_settings.friction_solver, integrator_settings.mu,
      integrator_settings.if_map, integrator_settings.spatial_grid_scale,
      integrator_settings.reduce_bandwidth, rendering_state);
}

bool RigidBody2DSceneParser::parseXMLSceneFile(
    const std::string &file_name, std::string &scripting_callback,
    RigidBody2DState &sim_state,
    RigidBody2DIntegratorSettings &integrator_settings,
    std::string &dt_string) {
  CameraSettings2D rendering_state;
  return parseXMLSceneFile(file_name, scripting_callback, sim_state,
                           integrator_settings, dt_string, rendering_state);
}

#ifdef USE_CAIRO
static bool parseCairoSettings(const rapidxml::xml_document<> &doc,
                               CairoRenderSettings &cairo_settings) {
  // Attempt to locate the root node
  if (doc.first_node("rigidbody2d_scene") == nullptr) {
    std::cerr
        << "Failed to locate root node rigidbody2d_scene in xml scene file."
        << std::endl;
    return false;
  }
  const rapidxml::xml_node<> &root_node{*doc.first_node("rigidbody2d_scene")};

  const rapidxml::xml_node<> *const nd = root_node.first_node("cairo_settings");
  if (nd == nullptr) {
    return true;
  }

  // Attempt to read the image dimensions
  Vector2i dims;
  {
    const rapidxml::xml_attribute<> *const attrib{
        nd->first_attribute("image_dimensions")};
    if (attrib == nullptr) {
      std::cerr << "Failed to locate image_dimensions attribute for "
                   "cairo_settings node."
                << std::endl;
      return false;
    }
    if (!StringUtilities::readScalarList(attrib->value(), 2, ' ', dims)) {
      std::cerr << "Failed to load image_dimensions attribute for "
                   "cairo_settings node, must provide 2 positive integers."
                << std::endl;
      return false;
    }
    if ((dims.array() <= 0).any()) {
      std::cerr << "Failed to load image_dimensions attribute for "
                   "cairo_settings node, must provide 2 positive integers."
                << std::endl;
      return false;
    }
  }

  // Attempt to read the camera scale
  double camera_scale;
  {
    const rapidxml::xml_attribute<> *const attrib{
        nd->first_attribute("camera_scale")};
    if (attrib == nullptr) {
      std::cerr
          << "Failed to locate camera_scale attribute for cairo_settings node."
          << std::endl;
      return false;
    }
    if (!StringUtilities::extractFromString(attrib->value(), camera_scale)) {
      std::cerr << "Failed to load camera_scale attribute for cairo_settings "
                   "node, must provide a scalar."
                << std::endl;
      return false;
    }
  }

  // Attempt to read the camera center
  Eigen::Vector2d camera_center;
  {
    const rapidxml::xml_attribute<> *const attrib{
        nd->first_attribute("camera_center")};
    if (attrib == nullptr) {
      std::cerr
          << "Failed to locate camera_center attribute for cairo_settings node."
          << std::endl;
      return false;
    }
    if (!StringUtilities::readScalarList(attrib->value(), 2, ' ',
                                         camera_center)) {
      std::cerr << "Failed to load camera_center attribute for cairo_settings "
                   "node, must provide 2 scalars."
                << std::endl;
      return false;
    }
  }

  // Attempt to read the background color
  Eigen::Vector4d bg_color;
  {
    const rapidxml::xml_attribute<> *const attrib{
        nd->first_attribute("background_color")};
    if (attrib == nullptr) {
      std::cerr << "Failed to locate background_color attribute for "
                   "cairo_settings node."
                << std::endl;
      return false;
    }
    if (!StringUtilities::readScalarList(attrib->value(), 4, ' ', bg_color)) {
      std::cerr << "Failed to load background_color attribute for "
                   "cairo_settings node, must provide 4 scalars."
                << std::endl;
      return false;
    }
    if ((bg_color.array() < 0.0).any() || (bg_color.array() > 1.0).any()) {
      std::cerr << "Failed to load background_color attribute for "
                   "cairo_settings node, values must be in [0, 1]."
                << std::endl;
      return false;
    }
  }

  // Attempt to read the circle color
  Eigen::Vector4d circle_color;
  {
    const rapidxml::xml_attribute<> *const attrib{
        nd->first_attribute("circle_color")};
    if (attrib == nullptr) {
      std::cerr
          << "Failed to locate circle_color attribute for cairo_settings node."
          << std::endl;
      return false;
    }
    if (!StringUtilities::readScalarList(attrib->value(), 4, ' ',
                                         circle_color)) {
      std::cerr << "Failed to load circle_color attribute for cairo_settings "
                   "node, must provide 4 scalars."
                << std::endl;
      return false;
    }
    if ((circle_color.array() < 0.0).any() ||
        (circle_color.array() > 1.0).any()) {
      std::cerr << "Failed to load circle_color attribute for cairo_settings "
                   "node, values must be in [0, 1]."
                << std::endl;
      return false;
    }
  }

  // Attempt to read the fixed circle color
  Eigen::Vector4d fixed_circle_color;
  {
    const rapidxml::xml_attribute<> *const attrib{
        nd->first_attribute("fixed_circle_color")};
    if (attrib == nullptr) {
      std::cerr << "Failed to locate fixed_circle_color attribute for "
                   "cairo_settings node."
                << std::endl;
      return false;
    }
    if (!StringUtilities::readScalarList(attrib->value(), 4, ' ',
                                         fixed_circle_color)) {
      std::cerr << "Failed to load fixed_circle_color attribute for "
                   "cairo_settings node, must provide 4 scalars."
                << std::endl;
      return false;
    }
    if ((fixed_circle_color.array() < 0.0).any() ||
        (fixed_circle_color.array() > 1.0).any()) {
      std::cerr << "Failed to load fixed_circle_color attribute for "
                   "cairo_settings node, values must be in [0, 1]."
                << std::endl;
      return false;
    }
  }

  cairo_settings.setImageDimensions(dims);
  cairo_settings.setCameraScale(camera_scale);
  cairo_settings.setCameraCenter(camera_center);
  cairo_settings.setBackgroundColor(bg_color);
  cairo_settings.setCircleColor(circle_color);
  cairo_settings.setFixedCircleColor(fixed_circle_color);

  return true;
}

bool RigidBody2DSceneParser::parseXMLSceneFile(
    const std::string &file_name, std::string &scripting_callback,
    RigidBody2DState &sim_state,
    std::unique_ptr<UnconstrainedMap> &unconstrained_map,
    std::string &dt_string, Rational<std::intmax_t> &dt, scalar &end_time,
    std::unique_ptr<ImpactOperator> &impact_operator,
    std::unique_ptr<ImpactMap> &impact_map, scalar &CoR,
    std::unique_ptr<FrictionSolver> &friction_solver, scalar &mu,
    std::unique_ptr<ImpactFrictionMap> &if_map, scalar &spatial_grid_scale,
    bool &reduce_bandwidth, CairoRenderSettings &cairo_settings) {
  // Attempt to load the xml document
  std::vector<char> xmlchars;
  rapidxml::xml_document<> doc;
  if (!loadXMLFile(file_name, xmlchars, doc)) {
    return false;
  }
  if (!parseCairoSettings(doc, cairo_settings)) {
    return false;
  }
  CameraSettings2D rendering_state;
  return basicParse(doc, file_name, scripting_callback, sim_state,
                    unconstrained_map, dt_string, dt, end_time, impact_operator,
                    impact_map, CoR, friction_solver, mu, if_map,
                    spatial_grid_scale, reduce_bandwidth, rendering_state);
}
#endif
