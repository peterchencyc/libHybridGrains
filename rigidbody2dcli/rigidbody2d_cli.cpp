// rigidbody2d_cli.cpp
//
// Breannan Smith
// Last updated: 10/01/2015

#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>

#include "scisim/CompileDefinitions.h"
#include "scisim/ConstrainedMaps/FrictionSolver.h"
#include "scisim/ConstrainedMaps/ImpactFrictionMap.h"
#include "scisim/ConstrainedMaps/ImpactMaps/ImpactMap.h"
#include "scisim/ConstrainedMaps/ImpactMaps/ImpactOperator.h"
#include "scisim/ConstrainedMaps/ImpactMaps/ImpactSolution.h"
#include "scisim/ConstrainedMaps/NullFrictionSolver.h"
#include "scisim/HDF5File.h"
#include "scisim/Math/MathDefines.h"
#include "scisim/Math/MathUtilities.h"
#include "scisim/Math/Rational.h"
#include "scisim/PythonTools.h"
#include "scisim/StringUtilities.h"
#include "scisim/Timer/TimeUtils.h"
#include "scisim/UnconstrainedMaps/NullUnconstrainedMap.h"
#include "scisim/UnconstrainedMaps/UnconstrainedMap.h"
#include "scisim/Utilities.h"

#include "rigidbody2d/CircleGeometry.h"
#include "rigidbody2d/PythonScripting.h"
#include "rigidbody2d/RigidBody2DSim.h"
#include "rigidbody2d/RigidBody2DUtilities.h"

#include "rigidbody2dutils/CameraSettings2D.h"
#include "rigidbody2dutils/RigidBody2DSceneParser.h"

#ifdef USE_CAIRO
#include "rigidbody2dutils/CairoImage.h"
#include "rigidbody2dutils/CairoRenderSettings.h"
#endif

static RigidBody2DSim g_sim;
static unsigned g_iteration = 0;
static std::unique_ptr<UnconstrainedMap> g_unconstrained_map{nullptr};
static Rational<std::intmax_t> g_dt;
static scalar g_end_time = SCALAR_NAN;
static std::unique_ptr<ImpactOperator> g_impact_operator{nullptr};
static scalar g_CoR = SCALAR_NAN;
static std::unique_ptr<FrictionSolver> g_friction_solver{nullptr};
static scalar g_mu = SCALAR_NAN;
static std::unique_ptr<ImpactMap> g_impact_map{nullptr};
static std::unique_ptr<ImpactFrictionMap> g_impact_friction_map{nullptr};
static PythonScripting g_scripting;
static bool g_reduce_bandwidth;

static std::string g_output_dir_name;
static bool g_output_forces{false};
// Number of timesteps between saves
static unsigned g_steps_per_save{0};
// Number of saves that been conducted so far
static unsigned g_output_frame{0};
static unsigned g_dt_string_precision{0};
static unsigned g_save_number_width{0};

static bool g_overwrite_snapshots{true};

#ifdef USE_CAIRO
static bool g_render_images = false;
static CairoRenderSettings g_render_settings;
#endif

// Magic number to print in front of binary output to aid in debugging
static const unsigned MAGIC_BINARY_NUMBER{8675309};

static std::string
generateOutputConfigurationDataFileName(const std::string &prefix,
                                        const std::string &extension) {
  std::stringstream ss;
  if (!g_output_dir_name.empty()) {
    ss << g_output_dir_name << "/";
  }
  ss << prefix << "_" << std::setfill('0') << std::setw(g_save_number_width)
     << g_output_frame << "." << extension;
  return ss.str();
}

static void printCompileInfo(std::ostream &output_stream) {
  output_stream << "Git Revision:     " << CompileDefinitions::GitSHA1
                << std::endl;
  output_stream << "Build Mode:       " << CompileDefinitions::BuildMode
                << std::endl;
  output_stream << "C Compiler:       " << CompileDefinitions::CCompiler
                << std::endl;
  output_stream << "C++ Compiler:     " << CompileDefinitions::CXXCompiler
                << std::endl;
#ifdef FORTRAN_FOUND
  output_stream << "Fortran Compiler: " << CompileDefinitions::FortranCompiler
                << std::endl;
#endif
}

static unsigned
computeTimestepDisplayPrecision(const Rational<std::intmax_t> &dt,
                                const std::string &dt_string) {
  if (dt_string.find('.') != std::string::npos) {
    return unsigned(
        StringUtilities::computeNumCharactersToRight(dt_string, '.'));
  } else {
    std::string converted_dt_string;
    std::stringstream ss;
    ss << std::fixed << scalar(dt);
    ss >> converted_dt_string;
    return unsigned(
        StringUtilities::computeNumCharactersToRight(converted_dt_string, '.'));
  }
}

static std::string xmlFilePath(const std::string &xml_file_name) {
  std::string path;
  std::string file_name;
  StringUtilities::splitAtLastCharacterOccurence(xml_file_name, path, file_name,
                                                 '/');
  if (file_name.empty()) {
    using std::swap;
    swap(path, file_name);
  }
  return path;
}

static bool loadXMLScene(const std::string &xml_file_name) {
  std::string scripting_callback_name;
  RigidBody2DState state;
  std::string dt_string;
  scalar spatial_grid_cell_scale;
  RigidBody2DState new_state;

#ifndef USE_CAIRO
  CameraSettings2D unused_camera_settings;
  const bool loaded_successfully{RigidBody2DSceneParser::parseXMLSceneFile(
      xml_file_name, scripting_callback_name, new_state, g_unconstrained_map,
      dt_string, g_dt, g_end_time, g_impact_operator, g_impact_map, g_CoR,
      g_friction_solver, g_mu, g_impact_friction_map, spatial_grid_cell_scale,
      g_reduce_bandwidth, unused_camera_settings)};
#else
  const bool loaded_successfully{RigidBody2DSceneParser::parseXMLSceneFile(
      xml_file_name, scripting_callback_name, new_state, g_unconstrained_map,
      dt_string, g_dt, g_end_time, g_impact_operator, g_impact_map, g_CoR,
      g_friction_solver, g_mu, g_impact_friction_map, spatial_grid_cell_scale,
      g_reduce_bandwidth, g_render_settings)};
#endif

  if (!loaded_successfully) {
    return false;
  }

  g_sim.state() = std::move(new_state);
  g_sim.grid().setCellScale(spatial_grid_cell_scale);

  g_dt_string_precision = computeTimestepDisplayPrecision(g_dt, dt_string);

  // Configure the scripting
  PythonScripting new_scripting{xmlFilePath(xml_file_name),
                                scripting_callback_name};
  swap(g_scripting, new_scripting);

  // User-provided start of simulation python callback
  g_scripting.setState(g_sim.state());
  g_scripting.startOfSimCallback();
  g_scripting.forgetState();

  return true;
}

static std::string generateSimulationTimeString() {
  std::stringstream time_stream;
  time_stream << std::fixed << std::setprecision(g_dt_string_precision)
              << g_iteration * scalar(g_dt);
  return time_stream.str();
}

static int saveState() {
#ifdef USE_HDF5
  // Generate a base filename
  const std::string output_file_name{
      generateOutputConfigurationDataFileName("config", "h5")};

  // Print a status message with the simulation time and output number
  std::cout << "Saving state at time " << generateSimulationTimeString()
            << " to " << output_file_name;
  std::cout << "        " << TimeUtils::currentTime() << std::endl;

  // Save the simulation state
  try {
    HDF5File output_file{output_file_name, HDF5AccessType::READ_WRITE};
    // Save the iteration and time step and time
    output_file.writeScalar("", "timestep", scalar(g_dt));
    output_file.writeScalar("", "iteration", g_iteration);
    output_file.writeScalar("", "time", scalar(g_dt) * g_iteration);
    // Save out the git hash
    output_file.writeString("", "git_hash", CompileDefinitions::GitSHA1);
    // Tag this as 2d rigid body simulation
    output_file.writeString("", "sim_type", "rigidbody2d");
    // Save the real time
    // output_file.createGroup( "/run_stats" );
    // output_file.writeString( "/run_stats", "real_time",
    // TimeUtils::currentTime() );
    // Write out the simulation data
    g_sim.writeBinaryState(output_file);
  } catch (const std::string &error) {
    std::cerr << error << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
#else
  std::cerr << "Error, state output requires HDF5 support." << std::endl;
  std::exit(EXIT_FAILURE);
#endif
}

#ifdef USE_CAIRO
static int renderCairoImage() {
  // Generate a base filename
  const std::string output_file_name{
      generateOutputConfigurationDataFileName("image", "png")};
  std::cout << "Saving image to: " << output_file_name << std::endl;

  CairoImage cairo_image(
      g_render_settings.imageDimensions(), g_render_settings.backgroundColor(),
      g_render_settings.cameraScale(), g_render_settings.cameraCenter());

  for (const RigidBody2DStaticDrum &drum : g_sim.getState().drums()) {
    cairo_image.drawDrum(drum.x(), drum.theta(), drum.r(),
                         Eigen::Vector4d(0.0, 0.0, 0.0, 1.0));
  }

  const unsigned nbodies = g_sim.getState().nbodies();
  for (unsigned bdy_idx = 0; bdy_idx < nbodies; bdy_idx++) {
    const std::unique_ptr<RigidBody2DGeometry> &geo =
        g_sim.getState().bodyGeometry(bdy_idx);
    switch (geo->type()) {
    case RigidBody2DGeometryType::CIRCLE: {
      const CircleGeometry &circle = static_cast<const CircleGeometry &>(*geo);
      const Vector2s center = g_sim.getState().q().segment<2>(3 * bdy_idx);
      // TODO: Handle fixed bodies here
      if (!g_sim.getState().fixed(bdy_idx)) {
        cairo_image.drawFilledCircle(center, circle.r(),
                                     g_render_settings.circleColor());
      } else {
        cairo_image.drawFilledCircle(center, circle.r(),
                                     g_render_settings.fixedCircleColor());
      }
      break;
    }
    case RigidBody2DGeometryType::BOX: {
      std::cerr << "Error, box rendering not coded up yet!" << std::endl;
      std::exit(EXIT_FAILURE);
      // break;
    }
    }
  }

  cairo_image.saveImage(output_file_name);

  return EXIT_SUCCESS;
}
#endif

static int exportConfigurationData() {
  assert(g_steps_per_save != 0);
  if (g_iteration % g_steps_per_save == 0) {
    if (!g_output_dir_name.empty()) {
      if (saveState() == EXIT_FAILURE) {
        return EXIT_FAILURE;
      }
#ifdef USE_CAIRO
      if (g_render_images) {
        if (renderCairoImage() == EXIT_FAILURE) {
          return EXIT_FAILURE;
        }
      }
#endif
    }

    ++g_output_frame;
  }
  return EXIT_SUCCESS;
}

#ifdef USE_HDF5
static std::string generateOutputConstraintForceDataFileName() {
  std::stringstream ss;
  assert(g_output_frame > 0);
  ss << g_output_dir_name << "/forces_" << std::setfill('0')
     << std::setw(g_save_number_width) << g_output_frame - 1 << ".h5";
  return ss.str();
}
#endif

static int stepSystem() {
  const unsigned next_iter{g_iteration + 1};

  HDF5File force_file;
  assert(g_steps_per_save != 0);
  if (g_output_forces && g_iteration % g_steps_per_save == 0) {
#ifdef USE_HDF5
    assert(!g_output_dir_name.empty());
    const std::string constraint_force_file_name{
        generateOutputConstraintForceDataFileName()};
    std::cout << "Saving forces at time " << generateSimulationTimeString()
              << " to " << constraint_force_file_name << std::endl;
    try {
      force_file.open(constraint_force_file_name, HDF5AccessType::READ_WRITE);
      // Save the iteration and time step and time
      force_file.writeScalar("", "timestep", scalar(g_dt));
      force_file.writeScalar("", "iteration", g_iteration);
      force_file.writeScalar("", "time", scalar(g_dt) * g_iteration);
      // Save out the git hash
      force_file.writeString("", "git_hash", CompileDefinitions::GitSHA1);
      // Save the real time
      // force_file.createGroup( "/run_stats" );
      // force_file.writeString( "/run_stats", "real_time",
      // TimeUtils::currentTime() );
    } catch (const std::string &error) {
      std::cerr << error << std::endl;
      return EXIT_FAILURE;
    }
#else
    std::cerr << "Error, force output requires HDF5 support." << std::endl;
    std::exit(EXIT_FAILURE);
#endif
  }

  if (g_impact_friction_map != nullptr &&
      g_impact_friction_map->name() == "penalty_impact_friction_map") {
    static NullUnconstrainedMap null_map;
    static NullFrictionSolver null_friction_solver;
    g_sim.flow(g_scripting, next_iter, g_dt, null_map, g_CoR, g_mu,
               null_friction_solver, *g_impact_friction_map,
               g_reduce_bandwidth);
  } else if (g_unconstrained_map == nullptr && g_impact_operator == nullptr &&
             g_impact_map == nullptr && g_friction_solver == nullptr &&
             g_impact_friction_map == nullptr) {
    // Nothing to do
  } else if (g_unconstrained_map != nullptr && g_impact_operator == nullptr &&
             g_impact_map == nullptr && g_friction_solver == nullptr &&
             g_impact_friction_map == nullptr) {
    g_sim.flow(g_scripting, next_iter, g_dt, *g_unconstrained_map);
  } else if (g_unconstrained_map != nullptr && g_impact_operator != nullptr &&
             g_impact_map != nullptr && g_friction_solver == nullptr &&
             g_impact_friction_map == nullptr) {
    assert(g_impact_map != nullptr);
    ImpactSolution impact_solution;
    if (force_file.is_open()) {
      g_impact_map->exportForcesNextStep(impact_solution);
    }
    g_sim.flow(g_scripting, next_iter, g_dt, *g_unconstrained_map,
               *g_impact_operator, g_CoR, *g_impact_map, g_reduce_bandwidth);
    if (force_file.is_open()) {
      try {
        impact_solution.writeSolution(force_file);
      } catch (const std::string &error) {
        std::cerr << error << std::endl;
        return EXIT_FAILURE;
      }
    }
  } else if (g_unconstrained_map != nullptr && g_impact_operator == nullptr &&
             g_impact_map == nullptr && g_friction_solver != nullptr &&
             g_impact_friction_map != nullptr) {
    if (force_file.is_open()) {
      g_impact_friction_map->exportForcesNextStep(force_file);
    }
    g_sim.flow(g_scripting, next_iter, g_dt, *g_unconstrained_map, g_CoR, g_mu,
               *g_friction_solver, *g_impact_friction_map, g_reduce_bandwidth);
  } else {
    std::cerr
        << "Impossible code path hit in stepSystem. This is a bug. Exiting."
        << std::endl;
    return EXIT_FAILURE;
  }

  ++g_iteration;

  return exportConfigurationData();
}

static int executeSimLoop() {
  if (exportConfigurationData() == EXIT_FAILURE) {
    return EXIT_FAILURE;
  }

  while (true) {
    // N.B. this will ocassionaly not trigger at the *exact* equal time due to
    // floating point errors
    if (g_iteration * scalar(g_dt) >= g_end_time) {
      // Take one final step to ensure we have force data for end time
      if (g_output_forces) {
        if (stepSystem() == EXIT_FAILURE) {
          return EXIT_FAILURE;
        }
      }
      // User-provided end of simulation python callback
      g_scripting.setState(g_sim.state());
      g_scripting.endOfSimCallback();
      g_scripting.forgetState();
      std::cout << "Simulation complete at time " << g_iteration * scalar(g_dt)
                << ". Exiting." << std::endl;
      return EXIT_SUCCESS;
    }
    TimingTools timing_tools;
    timing_tools.start();
    if (stepSystem() == EXIT_FAILURE) {
      return EXIT_FAILURE;
    }
    timing_tools.stop("");
  }
}

static void printUsage(const std::string &executable_name) {
  std::cout << "Usage: " << executable_name << " xml_scene_file_name [options]"
            << std::endl;
  std::cout << "Options are:" << std::endl;
  std::cout
      << "   -h/--help                : prints this help message and exits"
      << std::endl;
  std::cout << "   -i/--impulses            : saves impulses in addition to "
               "configuration if an output directory is set"
            << std::endl;
  std::cout << "   -e/--end scalar          : overrides the end time specified "
               "in the scene file"
            << std::endl;
  std::cout << "   -o/--output_dir dir      : saves simulation state to the "
               "given directory"
            << std::endl;
  std::cout << "   -f/--frequency integer   : rate at which to save simulation "
               "data, in Hz; ignored if no output directory specified"
            << std::endl;
#ifdef USE_CAIRO
  std::cout << "   -c/--cairo               : enable image rendering with Cairo"
            << std::endl;
#endif
}

static bool parseCommandLineOptions(int *argc, char ***argv,
                                    bool &help_mode_enabled,
                                    scalar &end_time_override,
                                    unsigned &output_frequency) {
  const struct option long_options[] = {
      {"help", no_argument, nullptr, 'h'},
      {"impulses", no_argument, nullptr, 'i'},
      {"end", required_argument, nullptr, 'e'},
      {"output_dir", required_argument, nullptr, 'o'},
      {"frequency", required_argument, nullptr, 'f'},
#ifdef USE_CAIRO
      {"cairo", no_argument, nullptr, 'c'},
#endif
      {nullptr, 0, nullptr, 0}};

  while (true) {
    int option_index = 0;
#ifdef USE_CAIRO
    const int c{getopt_long(*argc, *argv, "hips:r:e:o:f:c", long_options,
                            &option_index)};
#else
    const int c{getopt_long(*argc, *argv, "hips:r:e:o:f:", long_options,
                            &option_index)};
#endif
    if (c == -1) {
      break;
    }
    switch (c) {
    case 'h': {
      help_mode_enabled = true;
      break;
    }
    case 'i': {
      g_output_forces = true;
      break;
    }
    case 'e': {
      if (!StringUtilities::extractFromString(optarg, end_time_override)) {
        std::cerr << "Failed to read value for argument for -e/--end. Value "
                     "must be a positive scalar."
                  << std::endl;
        return false;
      }
      if (end_time_override <= 0) {
        std::cerr << "Failed to read value for argument for -e/--end. Value "
                     "must be a positive scalar."
                  << std::endl;
        return false;
      }
      break;
    }
    case 'o': {
      g_output_dir_name = optarg;
      break;
    }
    case 'f': {
      if (!StringUtilities::extractFromString(optarg, output_frequency)) {
        std::cerr << "Failed to read value for argument for -f/--frequency. "
                     "Value must be an unsigned integer."
                  << std::endl;
        return false;
      }
      break;
    }
#ifdef USE_CAIRO
    case 'c': {
      g_render_images = true;
      break;
    }
#endif
    case '?': {
      return false;
    }
    default: {
      std::cerr
          << "This is a bug in the command line parser. Please file a report."
          << std::endl;
      return false;
    }
    }
  }

  return true;
}

#ifdef USE_PYTHON
static void exitCleanup() { Py_Finalize(); }
#endif

int main(int argc, char **argv) {
  // Command line options
  bool help_mode_enabled{false};
  scalar end_time_override{-1};
  unsigned output_frequency{0};

  // Attempt to load command line options
  if (!parseCommandLineOptions(&argc, &argv, help_mode_enabled,
                               end_time_override, output_frequency)) {
    return EXIT_FAILURE;
  }

  // If the user requested help, print help and exit
  if (help_mode_enabled) {
    printUsage(argv[0]);
    return EXIT_SUCCESS;
  }

  // Check for impossible combinations of options
  if (g_output_forces && g_output_dir_name.empty()) {
    std::cerr << "Impulse output requires an output directory." << std::endl;
    return EXIT_FAILURE;
  }

#ifdef USE_PYTHON
  // Initialize the Python interpreter
  Py_SetProgramName(argv[0]);
  Py_Initialize();

  // Initialize a callback that will close down the interpreter
  atexit(exitCleanup);

  // Allow subsequent Python commands to use the sys module
  PythonTools::pythonCommand("import sys");

  // Prevent Python from intercepting the interrupt signal
  PythonTools::pythonCommand("import signal");
  PythonTools::pythonCommand("signal.signal( signal.SIGINT, signal.SIG_DFL )");

  // Initialize the callbacks
  PythonScripting::initializeCallbacks();
#endif

  // The user must provide the path to an xml scene file
  if (argc != optind + 1) {
    std::cerr << "Invalid arguments. Must provide a single xml scene file name."
              << std::endl;
    return EXIT_FAILURE;
  }

  // Attempt to load the user-provided scene
  if (!loadXMLScene(std::string{argv[optind]})) {
    return EXIT_FAILURE;
  }

  // Override the default end time with the requested one, if provided
  if (end_time_override > 0.0) {
    g_end_time = end_time_override;
  }

  // Compute the data output rate
  assert(g_dt.positive());
  // If the user provided an output frequency
  if (output_frequency != 0) {
    const Rational<std::intmax_t> potential_steps_per_frame{
        std::intmax_t(1) / (g_dt * std::intmax_t(output_frequency))};
    if (!potential_steps_per_frame.isInteger()) {
      std::cerr << "Timestep and output frequency do not yield an integer "
                   "number of timesteps for data output. Exiting."
                << std::endl;
      return EXIT_FAILURE;
    }
    g_steps_per_save = unsigned(potential_steps_per_frame.numerator());
  }
  // Otherwise default to dumping every frame
  else {
    g_steps_per_save = 1;
  }
  assert(g_end_time > 0.0);
  g_save_number_width = MathUtilities::computeNumDigits(
      1 + unsigned(ceil(g_end_time / scalar(g_dt))) / g_steps_per_save);

  printCompileInfo(std::cout);
  assert(g_sim.state().q().size() % 3 == 0);
  std::cout << "Body count: " << g_sim.state().q().size() / 3 << std::endl;

  // If there are any intitial collisions, warn the user
  //{
  //  std::map<std::string,unsigned> collision_counts;
  //  std::map<std::string,scalar> collision_depths;
  //  g_sim.computeNumberOfCollisions( collision_counts, collision_depths );
  //  assert( collision_counts.size() == collision_depths.size() );
  //  if( !collision_counts.empty() ) { std::cout << "Warning, initial
  //  collisions detected (name : count : total_depth):" << std::endl; } for(
  //  const auto& count_pair : collision_counts )
  //  {
  //    const std::string& constraint_name = count_pair.first;
  //    const unsigned& constraint_count = count_pair.second;
  //    assert( collision_depths.find( constraint_name ) !=
  //    collision_depths.end() ); const scalar& constraint_depth =
  //    collision_depths[constraint_name]; std::string depth_string; if(
  //    !std::isnan( constraint_depth ) )
  //    {
  //      depth_string = StringUtilities::convertToString( constraint_depth );
  //    }
  //    else
  //    {
  //      depth_string = "depth_computation_not_supported";
  //    }
  //    std::cout << "   " << constraint_name << " : " << constraint_count << "
  //    : " << depth_string << std::endl;
  //  }
  //}

  if (g_end_time == SCALAR_INFINITY) {
    std::cout << "No end time specified. Simulation will run indefinitely."
              << std::endl;
  }

  return executeSimLoop();
}
