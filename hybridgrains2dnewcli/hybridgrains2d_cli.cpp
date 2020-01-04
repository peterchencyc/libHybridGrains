#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>

#include "scisim/CompileDefinitions.h"
#include "scisim/HDF5File.h"
#include "scisim/Math/Rational.h"
#include "scisim/PythonTools.h"
#include "scisim/StringUtilities.h"
#include "scisim/Timer/TimeUtils.h"

#include "hybridgrains2dnewutils/HybridGrains2DSceneParser.h"

#include "hybridgrains2dnew/DiscreteIntegrator.h"
#include "hybridgrains2dnew/HybridDefinitions.h"
#include "hybridgrains2dnew/HybridGrains2DSim.h"
#include "hybridgrains2dnew/MPMIntegrator.h"

#include "rigidbody2d/RigidBody2DIntegratorSettings.h"
#include "rigidbody2d/RigidBody2DState.h"

#include "rigidbody2dutils/CameraSettings2D.h"
#include "rigidbody2dutils/RigidBody2DSceneParser.h"

#include "mpmgrains2d/InitialSimulationState.h"

#include "mpmgrains2dutils/MPMGrains2DParser.h"

static HybridGrains2DSim g_sim;

static std::string g_output_dir_name;
static bool g_output_forces{false};
// Number of timesteps between saves
static unsigned g_steps_per_save{0};
// Number of saves that been conducted so far
static unsigned g_output_frame{0};
static unsigned g_dt_string_precision{0};
static unsigned g_save_number_width{0};

static bool g_serialize_snapshots{false};
static bool g_overwrite_snapshots{true};

// Magic number to print in front of binary output to aid in debugging
static constexpr unsigned MAGIC_BINARY_NUMBER{8675309};

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
  output_stream << "CXX Compiler:     " << CompileDefinitions::CXXCompiler
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

static std::string generateSimulationTimeString() {
  std::stringstream time_stream;
  time_stream << std::fixed << std::setprecision(g_dt_string_precision)
              << scalar(
                     std::intmax_t(g_sim.integratorState().overallIteration()) *
                     g_sim.integratorState().overallTimestep());
  return time_stream.str();
}

#ifdef USE_HDF5
static int saveState() {
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
    output_file.writeScalar("", "global_timestep",
                            scalar(g_sim.integratorState().overallTimestep()));
    output_file.writeScalar("", "global_iteration",
                            g_sim.integratorState().overallIteration());
    output_file.writeScalar(
        "", "global_time",
        scalar(std::intmax_t(g_sim.integratorState().overallIteration()) *
               g_sim.integratorState().overallTimestep()));
    // Save out the git hash
    output_file.writeString("", "git_hash", CompileDefinitions::GitSHA1);
    // Tag this as 2d hybrid simulation
    output_file.writeString("", "sim_type", "hybrid2d");
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
}
#endif

static bool loadXMLScene(const std::string &xml_file_name) {
  // Load the names of the discrete and continuum configuration files
  HybridIntegratorSettings integrator_settings;
  if (!HybridGrains2DSceneParser::parseXMLSceneFile(xml_file_name,
                                                    integrator_settings)) {
    return false;
  }

  // Load the discrete state
  std::cout << "Loading discrete config file: "
            << integrator_settings.discrete_file_name << std::endl;
  std::string dt_string;
  std::string scripting_callback_name;
  RigidBody2DState new_state;
  RigidBody2DIntegratorSettings new_integrator_settings;
  const bool discrete_loaded_successfully{
      RigidBody2DSceneParser::parseXMLSceneFile(
          integrator_settings.discrete_file_name, scripting_callback_name,
          new_state, new_integrator_settings, dt_string)};
  if (!discrete_loaded_successfully) {
    std::cerr << "Failed to load discrete configuration file: "
              << integrator_settings.discrete_file_name << std::endl;
    return false;
  }
  std::cout
      << "!!!!!!!!!! DEBUG INFORMATION: Initial total simulated discrete mass: "
      << new_state.totalSimulatedMass() << std::endl;

  // Remove masked-out discrete bodies
  new_state.removeBodiesIntersectingBoxes(
      integrator_settings.discrete_body_masks, nullptr);

  // Load the continuum state
  std::cout << "Loading continuum config file: "
            << integrator_settings.continuum_file_name << std::endl;
  InitialSimulationState initial_continuum_state;
  const bool continuum_loaded_successfully{MPMGrains2DParser::readXMLFile(
      integrator_settings.continuum_file_name, initial_continuum_state)};
  if (!continuum_loaded_successfully) {
    std::cerr << "Failed to load continuum configuration file: "
              << integrator_settings.continuum_file_name << std::endl;
    return false;
  }

  // Set the simulation state
  {
    const DiscreteIntegrator discrete_integrator(
        new_integrator_settings.dt, new_integrator_settings.unconstrained_map,
        new_integrator_settings.impact_operator,
        new_integrator_settings.friction_solver, new_integrator_settings.if_map,
        new_integrator_settings.impact_map, new_integrator_settings.CoR,
        new_integrator_settings.mu, new_integrator_settings.reduce_bandwidth);
    const MPMIntegrator mpm_integrator(initial_continuum_state.dt);
    const SimulationState test_state{
        initial_continuum_state.generateSimulationState()};
    g_sim = HybridGrains2DSim(
        HybridIntegratorState(
            integrator_settings.overall_dt, integrator_settings.end_time,
            discrete_integrator, mpm_integrator,
            integrator_settings.integrator_style,
            integrator_settings.kinematically_scripted_hybrid_fronts,
            integrator_settings.poorman_settings),
        new_state, initial_continuum_state.generateSimulationState());

#ifdef USE_HDF5
    saveState();
#endif

    g_sim.initializeAvoidAVoid();

    // Initialize the discrete scripting callback
    {
      std::string path;
      std::string file_name;
      StringUtilities::splitAtLastCharacterOccurence(xml_file_name, path,
                                                     file_name, '/');
      if (file_name.empty()) {
        using std::swap;
        swap(path, file_name);
      }
      g_sim.integratorState().discreteIntegrator().setPythonCallback(
          path, scripting_callback_name);
      g_sim.integratorState().discreteIntegrator().pythonStartOfSim(
          g_sim.discreteSim());
    }
  }

  // Compute the number of characters after the decimal point in the timestep
  // string
  g_dt_string_precision = computeTimestepDisplayPrecision(
      integrator_settings.overall_dt, integrator_settings.time_step_string);

  // checkForInitialCollisions( g_sim.discreteSim() );

  std::cout << "!!!!!!!!!! DEBUG INFORMATION: Initialized discrete mass: "
            << g_sim.discreteSim().getState().totalSimulatedMass() << std::endl;
  std::cout << "!!!!!!!!!! DEBUG INFORMATION: Initialized continuum mass: "
            << g_sim.continuumState().material_points.totalMass() << std::endl;

  return true;
}

static int serializeSystem() {
  // Generate a base filename
  const std::string serialized_file_name{
      g_overwrite_snapshots
          ? "serial.bin"
          : generateOutputConfigurationDataFileName("serial", "bin")};

  // Print a message to the user that the state is being written
  std::cout << "Serializing: " << generateSimulationTimeString() << " to "
            << serialized_file_name;
  std::cout << "        " << TimeUtils::currentTime() << std::endl;

  // Attempt to open the output file
  std::ofstream serial_stream{serialized_file_name, std::ios::binary};
  if (!serial_stream.is_open()) {
    std::cerr << "Failed to open serialization file: " << serialized_file_name
              << std::endl;
    std::cerr << "Exiting." << std::endl;
    return EXIT_FAILURE;
  }

  // Write the magic number
  Utilities::serializeBuiltInType(MAGIC_BINARY_NUMBER, serial_stream);

  // Write the git revision
  {
    const std::string git_revision{CompileDefinitions::GitSHA1};
    StringUtilities::serializeString(git_revision, serial_stream);
  }

  // Write the actual state
  g_sim.serialize(serial_stream);
  StringUtilities::serializeString(g_output_dir_name, serial_stream);
  Utilities::serializeBuiltInType(g_output_forces, serial_stream);
  Utilities::serializeBuiltInType(g_steps_per_save, serial_stream);
  Utilities::serializeBuiltInType(g_output_frame, serial_stream);
  Utilities::serializeBuiltInType(g_dt_string_precision, serial_stream);
  Utilities::serializeBuiltInType(g_save_number_width, serial_stream);
  Utilities::serializeBuiltInType(g_serialize_snapshots, serial_stream);
  Utilities::serializeBuiltInType(g_overwrite_snapshots, serial_stream);

  return EXIT_SUCCESS;
}

static int deserializeSystem(const std::string &file_name) {
  std::cout << "Loading serialized simulation state file: " << file_name
            << std::endl;

  // Attempt to open the input file
  std::ifstream serial_stream{file_name, std::ios::binary};
  if (!serial_stream.is_open()) {
    std::cerr << "Failed to open serialized state in file: " << file_name
              << std::endl;
    std::cerr << "Exiting." << std::endl;
    return EXIT_FAILURE;
  }

  // Verify the magic number
  if (Utilities::deserialize<unsigned>(serial_stream) != MAGIC_BINARY_NUMBER) {
    std::cerr << "File " << file_name
              << " does not appear to be a serialized 2D SCISim rigid body "
                 "simulation. Exiting."
              << std::endl;
    return EXIT_FAILURE;
  }

  // Read the git revision
  {
    const std::string git_revision{
        StringUtilities::deserializeString(serial_stream)};
    if (CompileDefinitions::GitSHA1 != git_revision) {
      std::cerr
          << "Warning, resuming from data file for a different git revision."
          << std::endl;
      std::cerr << "   Serialized Git Revision: " << git_revision << std::endl;
      std::cerr << "      Current Git Revision: " << CompileDefinitions::GitSHA1
                << std::endl;
    }
    std::cout << "Git Revision: " << git_revision << std::endl;
  }

  g_sim.deserialize(serial_stream);
  g_output_dir_name = StringUtilities::deserializeString(serial_stream);
  g_output_forces = Utilities::deserialize<bool>(serial_stream);
  g_steps_per_save = Utilities::deserialize<unsigned>(serial_stream);
  g_output_frame = Utilities::deserialize<unsigned>(serial_stream);
  g_dt_string_precision = Utilities::deserialize<unsigned>(serial_stream);
  g_save_number_width = Utilities::deserialize<unsigned>(serial_stream);
  g_serialize_snapshots = Utilities::deserialize<bool>(serial_stream);
  g_overwrite_snapshots = Utilities::deserialize<bool>(serial_stream);

  return EXIT_SUCCESS;
}

static int exportConfigurationData() {
  assert(g_steps_per_save != 0);
  if (g_sim.integratorState().overallIteration() % g_steps_per_save == 0) {
#ifdef USE_HDF5
    if (!g_output_dir_name.empty()) {
      if (saveState() == EXIT_FAILURE) {
        return EXIT_FAILURE;
      }
    }
#endif
    if (g_serialize_snapshots) {
      if (serializeSystem() == EXIT_FAILURE) {
        return EXIT_FAILURE;
      }
    }
    g_output_frame++;
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
  HDF5File force_file;
  assert(g_steps_per_save != 0);
  if (g_output_forces &&
      g_sim.integratorState().overallIteration() % g_steps_per_save == 0) {
#ifdef USE_HDF5
    assert(!g_output_dir_name.empty());
    const std::string constraint_force_file_name{
        generateOutputConstraintForceDataFileName()};
    std::cout << "Saving forces at time " << generateSimulationTimeString()
              << " to " << constraint_force_file_name << std::endl;
    try {
      force_file.open(constraint_force_file_name, HDF5AccessType::READ_WRITE);
      // Save the iteration and time step and time
      force_file.writeScalar("", "discrete_timestep",
                             scalar(g_sim.integratorState().overallTimestep()));
      force_file.writeScalar("", "global_iteration",
                             g_sim.integratorState().overallIteration());
      force_file.writeScalar(
          "", "global_time",
          scalar(g_sim.integratorState().overallTimestep() *
                 std::intmax_t(g_sim.integratorState().overallIteration())));
      // Save out the git hash
      force_file.writeString("", "git_hash", CompileDefinitions::GitSHA1);
    } catch (const std::string &error) {
      std::cerr << error << std::endl;
      return EXIT_FAILURE;
    }
#else
    std::cerr << "Error, force output requires HDF5 support." << std::endl;
    std::exit(EXIT_FAILURE);
#endif
  }

  if (force_file.is_open()) {
    if (g_sim.integratorState().discreteIntegrator().frictionIsEnabled()) {
      g_sim.integratorState()
          .discreteIntegrator()
          .impactFrictionMap()
          ->exportForcesNextStep(force_file);
    } else {
      std::cerr << "Error, must have an impact friction map for force output."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  g_sim.stepSystem();

  return exportConfigurationData();
}

static int executeSimLoop() {
  if (exportConfigurationData() == EXIT_FAILURE) {
    return EXIT_FAILURE;
  }

  while (true) {
    if (std::intmax_t(g_sim.integratorState().overallIteration()) *
            g_sim.integratorState().overallTimestep() >=
        g_sim.integratorState().endTime()) {
      // Take one final step to ensure we have force data for end time
      if (g_output_forces) {

        if (stepSystem() == EXIT_FAILURE) {
          return EXIT_FAILURE;
        }
      }
      // User-provided end of simulation python callback
      g_sim.integratorState().discreteIntegrator().pythonEndOfSim(
          g_sim.discreteSim());
      std::cout << "Simulation complete at time "
                << scalar(std::intmax_t(
                              g_sim.integratorState().overallIteration()) *
                          g_sim.integratorState().overallTimestep())
                << ". Exiting." << std::endl;
      return EXIT_SUCCESS;
    }

    TimingTools timing_tools;
    timing_tools.start();
    if (stepSystem() == EXIT_FAILURE) {
      return EXIT_FAILURE;
    }
    std::cout << "overall" << std::endl;
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
  std::cout << "   -r/--resume file         : resumes the simulation from a "
               "serialized file"
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
  std::cout << "   -s/--serialize_snapshots bool : save a bit identical, "
               "resumable snapshot; if 0 overwrites the snapshot each "
               "timestep, if 1 saves a new snapshot for each timestep"
            << std::endl;
}

static bool parseCommandLineOptions(int *argc, char ***argv,
                                    bool &help_mode_enabled,
                                    Rational<std::intmax_t> &end_time_override,
                                    unsigned &output_frequency,
                                    std::string &serialized_file_name) {
  const struct option long_options[] = {
      {"help", no_argument, 0, 'h'},
      {"impulses", no_argument, 0, 'i'},
      {"serialize_snapshots", required_argument, 0, 's'},
      {"resume", required_argument, 0, 'r'},
      {"end", required_argument, 0, 'e'},
      {"output_dir", required_argument, 0, 'o'},
      {"frequency", required_argument, 0, 'f'},
      {0, 0, 0, 0}};

  while (true) {
    int option_index = 0;
    const int c{
        getopt_long(*argc, *argv, "his:r:e:o:f:", long_options, &option_index)};
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
    case 's': {
      g_serialize_snapshots = true;
      if (!StringUtilities::extractFromString(optarg, g_overwrite_snapshots)) {
        std::cerr << "Failed to read value for argument for "
                     "-s/--serialize_snapshots. Value must be a boolean."
                  << std::endl;
        return false;
      }
      g_overwrite_snapshots = !g_overwrite_snapshots;
      break;
    }
    case 'r': {
      serialized_file_name = optarg;
      break;
    }
    case 'e': {
      if (!RationalTools::extractFromString(optarg, end_time_override)) {
        std::cerr << "Failed to read value for argument for -e/--end. Value "
                     "must be a positive scalar."
                  << std::endl;
        return false;
      }
      if (!end_time_override.positive()) {
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
  Rational<std::intmax_t> end_time_override{-1, 1};
  unsigned output_frequency{0};
  std::string serialized_file_name;

  // Attempt to load command line options
  if (!parseCommandLineOptions(&argc, &argv, help_mode_enabled,
                               end_time_override, output_frequency,
                               serialized_file_name)) {
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

  if (!serialized_file_name.empty()) {
    if (deserializeSystem(serialized_file_name) == EXIT_FAILURE) {
      return EXIT_FAILURE;
    }
    return executeSimLoop();
  }

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
  if (end_time_override.positive()) {
    g_sim.integratorState().setEndTime(end_time_override);
  }

  // Compute the data output rate
  assert(g_sim.integratorState().overallTimestep().positive());
  // If the user provided an output frequency
  if (output_frequency != 0) {
    const Rational<std::intmax_t> potential_steps_per_frame{
        std::intmax_t(1) / (g_sim.integratorState().overallTimestep() *
                            std::intmax_t(output_frequency))};
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
  assert(g_sim.integratorState().endTime().positive());
  g_save_number_width = MathUtilities::computeNumDigits(
      1 + unsigned(ceil(scalar(g_sim.integratorState().endTime() /
                               g_sim.integratorState().overallTimestep()))) /
              g_steps_per_save);

  std::cout << "Discrete time step: "
            << g_sim.integratorState().discreteIntegrator().timestep()
            << std::endl;
  std::cout << "Continuum time step: "
            << g_sim.integratorState().continuumIntegrator().timestep()
            << std::endl;
  printCompileInfo(std::cout);
  assert(g_sim.discreteState().q().size() % 3 == 0);
  std::cout << "Body count: " << g_sim.discreteState().q().size() / 3
            << std::endl;

  return executeSimLoop();
}
