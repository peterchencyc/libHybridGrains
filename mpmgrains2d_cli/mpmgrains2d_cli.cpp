#include <iostream>
#include <iomanip>
#include <fstream>
#include <getopt.h>

#include "mpmgrains2dutils/MPMGrains2DParser.h"

#include "scisim/Math/MathDefines.h"
#include "scisim/Math/MathUtilities.h"
#include "scisim/Math/Rational.h"
#include "scisim/StringUtilities.h"
#include "mpmgrains2d/InitialSimulationState.h"
#include "mpmgrains2d/SimulationState.h"
#include "mpmgrains2d/ExplicitIntegrator.h"
#include "scisim/Timer/TimeUtils.h"
#include "scisim/Utilities.h"
#include "mpmgrains2d/CompileDefinitions.h"

#ifdef USE_HDF5
#include "scisim/HDF5File.h"
#endif

static SimulationState g_state;

static unsigned g_iteration;
static Rational<std::intmax_t> g_dt;
static scalar g_end_time;

#ifdef USE_HDF5
static std::string g_output_dir_name;
#endif

// Number of timesteps between saves
static unsigned g_steps_per_save;
// Number of saves that been conducted so far
static unsigned g_output_frame;
static unsigned g_dt_string_precision;
static unsigned g_save_number_width;

static bool g_serialize_snapshots;
static bool g_overwrite_snapshots;

// Magic number to print in front of binary output to aid in debugging
static constexpr unsigned MAGIC_BINARY_NUMBER{ 9035768 };

static std::string generateOutputFileName( const std::string& prefix, const std::string& extension )
{
  std::stringstream ss;
  #ifdef USE_HDF5
  if( !g_output_dir_name.empty() )
  {
    ss << g_output_dir_name << "/";
  }
  #endif
  ss << prefix << "_" << std::setfill('0') << std::setw( int(g_save_number_width) ) << g_output_frame << "." << extension;
  return ss.str();
}

//#ifdef USE_HDF5
//static std::string generateOutputFileNameWithoutDir( const std::string& prefix, const std::string& extension )
//{
//  std::stringstream ss;
//  ss << prefix << "_" << std::setfill('0') << std::setw( int(g_save_number_width) ) << g_output_frame << "." << extension;
//  return ss.str();
//}
//#endif

static void printCompileInfo( std::ostream& output_stream )
{
  output_stream << "Git Revision:     " << CompileDefinitions::GitSHA1 << std::endl;
  output_stream << "Build Mode:       " << CompileDefinitions::BuildMode << std::endl;
  output_stream << "C Compiler:       " << CompileDefinitions::CCompiler << std::endl;
  output_stream << "C++ Compiler:     " << CompileDefinitions::CXXCompiler << std::endl;
  output_stream << "Build Date:       " << CompileDefinitions::BuildDateTime << std::endl;
}

static unsigned computeTimestepDisplayPrecision( const Rational<std::intmax_t>& dt, const std::string& dt_string )
{
  if( dt_string.find( '.' ) != std::string::npos )
  {
    return unsigned( StringUtilities::computeNumCharactersToRight( dt_string, '.' ) );
  }
  else
  {
    std::string converted_dt_string;
    std::stringstream ss;
    ss << std::fixed << scalar( dt );
    ss >> converted_dt_string;
    return unsigned( StringUtilities::computeNumCharactersToRight( converted_dt_string, '.' ) );
  }
}

static bool loadXMLScene( const std::string& xml_file_name )
{
  InitialSimulationState initial_state;

  if( !MPMGrains2DParser::readXMLFile( xml_file_name, initial_state ) )
  {
    return false;
  }

  // Current timestep
  g_dt = initial_state.dt;
  // End time of the simulation
  g_end_time = initial_state.end_time;

  g_state = initial_state.generateSimulationState();

  // Compute the number of characters after the decimal point in the timestep string
  g_dt_string_precision = computeTimestepDisplayPrecision( g_dt, initial_state.dt_string );

  return true;
}

static std::string generateSimulationTimeString()
{
  std::stringstream time_stream;
  time_stream << std::fixed << std::setprecision( int(g_dt_string_precision) ) << g_iteration * scalar( g_dt );
  return time_stream.str();
}

#ifdef USE_HDF5
static int saveState()
{
  // Generate a base filename
  const std::string output_file_name{ generateOutputFileName( "config", "h5" ) };

  // Print a status message with the simulation time and output number
  std::cout << "Saving state at time " << generateSimulationTimeString() << " to " << output_file_name;
  std::cout << "        " << TimeUtils::currentTime() << std::endl;

  // Save the simulation state
  try
  {
    HDF5File output_file{ output_file_name, HDF5AccessType::READ_WRITE };
    // Save the iteration and time step and time
    output_file.writeScalar( "", "timestep", scalar( g_dt ) );
    output_file.writeScalar( "", "iteration", g_iteration );
    output_file.writeScalar( "", "time", scalar( g_dt ) * g_iteration );
    // Save out the git hash
    output_file.writeString( "", "git_hash", CompileDefinitions::GitSHA1 );
    // Write out the simulation data
    g_state.writeBinaryState( "", output_file );
  }
  catch( const std::string& error )
  {
    std::cerr << error << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
#endif

static int serializeSystem()
{
  // Generate a base filename
  const std::string serialized_file_name{ g_overwrite_snapshots ? "serial.bin" : generateOutputFileName( "serial", "bin" ) };

  // Print a message to the user that the state is being written
  std::cout << "Serializing: " << generateSimulationTimeString() << " to " << serialized_file_name;
  std::cout << "        " << TimeUtils::currentTime() << std::endl;

  // Attempt to open the output file
  std::ofstream serial_stream( serialized_file_name, std::ios::binary );
  if( !serial_stream.is_open() )
  {
    std::cerr << "Failed to open serialization file: " << serialized_file_name << std::endl;
    std::cerr << "Exiting." << std::endl;
    return EXIT_FAILURE;
  }

  // Write the magic number
  Utilities::serializeBuiltInType( MAGIC_BINARY_NUMBER, serial_stream );

  // Write the git revision
  {
    StringUtilities::serializeString( CompileDefinitions::GitSHA1, serial_stream );
  }

  // Write the actual state
  g_state.serialize( serial_stream );
  Utilities::serializeBuiltInType( g_iteration, serial_stream );
  RationalTools::serialize( g_dt, serial_stream );
  Utilities::serializeBuiltInType( g_end_time, serial_stream );
  #ifdef USE_HDF5
  StringUtilities::serializeString( g_output_dir_name, serial_stream );
  #endif
  Utilities::serializeBuiltInType( g_steps_per_save, serial_stream );
  Utilities::serializeBuiltInType( g_output_frame, serial_stream );
  Utilities::serializeBuiltInType( g_dt_string_precision, serial_stream );
  Utilities::serializeBuiltInType( g_save_number_width, serial_stream );
  Utilities::serializeBuiltInType( g_serialize_snapshots, serial_stream );
  Utilities::serializeBuiltInType( g_overwrite_snapshots, serial_stream );

  return EXIT_SUCCESS;
}

static int deserializeSystem( const std::string& file_name )
{
  std::cout << "Loading serialized simulation state file: " << file_name << std::endl;

  // Attempt to open the input file
  std::ifstream serial_stream{ file_name, std::ios::binary };
  if( !serial_stream.is_open() )
  {
    std::cerr << "Failed to open serialized state in file: " << file_name << std::endl;
    std::cerr << "Exiting." << std::endl;
    return EXIT_FAILURE;
  }

  // Verify the magic number
  if( Utilities::deserialize<unsigned>( serial_stream ) != MAGIC_BINARY_NUMBER )
  {
    std::cerr << "File " << file_name << " does not appear to be a serialized 2D MPM simulation. Exiting." << std::endl;
    return EXIT_FAILURE;
  }

  // Read the git revision
  {
    const std::string git_revision{ StringUtilities::deserializeString( serial_stream ) };
    if( CompileDefinitions::GitSHA1 != git_revision )
    {
      std::cerr << "Warning, resuming from data file for a different git revision." << std::endl;
      std::cerr << "   Serialized Git Revision: " << git_revision << std::endl;
      std::cerr << "      Current Git Revision: " << CompileDefinitions::GitSHA1 << std::endl;
    }
    std::cout << "Git Revision: " << git_revision << std::endl;
  }

  g_state.deserialize( serial_stream );
  g_iteration = Utilities::deserialize<unsigned>( serial_stream );
  RationalTools::deserialize( g_dt, serial_stream );
  assert( g_dt.positive() );
  g_end_time = Utilities::deserialize<scalar>( serial_stream );
  assert( g_end_time > 0.0 );
  #ifdef USE_HDF5
  g_output_dir_name = StringUtilities::deserializeString( serial_stream );
  #endif
  g_steps_per_save = Utilities::deserialize<unsigned>( serial_stream );
  g_output_frame = Utilities::deserialize<unsigned>( serial_stream );
  g_dt_string_precision = Utilities::deserialize<unsigned>( serial_stream );
  g_save_number_width = Utilities::deserialize<unsigned>( serial_stream );
  g_serialize_snapshots = Utilities::deserialize<bool>( serial_stream );
  g_overwrite_snapshots = Utilities::deserialize<bool>( serial_stream );

  return EXIT_SUCCESS;
}

static int exportConfigurationData()
{
   assert( g_steps_per_save != 0 );
   if( g_iteration % g_steps_per_save == 0 )
   {
     #ifdef USE_HDF5
     if( !g_output_dir_name.empty() )
     {
       if( saveState() == EXIT_FAILURE )
       {
         return EXIT_FAILURE;
       }
     }
     #endif
     if( g_serialize_snapshots )
     {
       if( serializeSystem() == EXIT_FAILURE )
       {
         return EXIT_FAILURE;
       }
     }
     g_output_frame++;
   }
  return EXIT_SUCCESS;
}

static int stepSystem()
{
  const unsigned next_iter{ g_iteration + 1 };

  std::cout << "Stepping to: " << next_iter * scalar(g_dt) << " of " << g_end_time << std::endl;

  ExplicitIntegrator::flow( scalar(g_dt), g_state );

  g_iteration++;

  return exportConfigurationData();
}

static int executeSimLoop()
{
  if( exportConfigurationData() == EXIT_FAILURE )
  {
    return EXIT_FAILURE;
  }

  while( true )
  {
    // N.B. this will ocassionaly not trigger at the *exact* equal time due to floating point errors
    if( g_iteration * scalar( g_dt ) >= g_end_time )
    {
      std::cout << "Simulation complete at time " << g_iteration * scalar( g_dt ) << ". Exiting." << std::endl;
      return EXIT_SUCCESS;
    }

    if( stepSystem() == EXIT_FAILURE )
    {
      return EXIT_FAILURE;
    }
  }
}

static void printUsage( const std::string& executable_name )
{
  std::cout << "Usage: " << executable_name << " xml_scene_file_name [options]" << std::endl;
  std::cout << "Options are:" << std::endl;
  std::cout << "   -h/--help                : prints this help message and exits" << std::endl;
  std::cout << "   -r/--resume file         : resumes the simulation from a serialized file" << std::endl;
  std::cout << "   -e/--end scalar          : overrides the end time specified in the scene file" << std::endl;
  #ifdef USE_HDF5
  std::cout << "   -o/--output_dir dir      : saves simulation state to the given directory" << std::endl;
  #endif
  std::cout << "   -f/--frequency integer   : rate at which to save simulation data, in Hz; ignored if no output directory specified" << std::endl;
  std::cout << "   -s/--serialize_snapshots bool : save a bit identical, resumable snapshot; if 0 overwrites the snapshot each timestep, if 1 saves a new snapshot for each timestep" << std::endl;
}

static bool parseCommandLineOptions( int* argc, char*** argv, bool& help_mode_enabled, scalar& end_time_override, unsigned& output_frequency, std::string& serialized_file_name )
{
  const struct option long_options[] =
  {
    { "help", no_argument, nullptr, 'h' },
    { "serialize_snapshots", required_argument, nullptr, 's' },
    { "resume", required_argument, nullptr, 'r' },
    { "end", required_argument, nullptr, 'e' },
    #ifdef USE_HDF5
    { "output_dir", required_argument, nullptr, 'o' },
    #endif
    { "frequency", required_argument, nullptr, 'f' },
    { nullptr, 0, nullptr, 0 }
  };

  while( true )
  {
    int option_index = 0;
    #ifndef USE_HDF5
    const int c{ getopt_long( *argc, *argv, "hs:r:e:f:", long_options, &option_index ) };
    #else
    const int c{ getopt_long( *argc, *argv, "hs:r:e:o:f:", long_options, &option_index ) };
    #endif
    if( c == -1 )
    {
      break;
    }
    switch( c )
    {
      case 'h':
      {
        help_mode_enabled = true;
        break;
      }
      case 's':
      {
        g_serialize_snapshots = true;
        if( !StringUtilities::extractFromString( optarg, g_overwrite_snapshots ) )
        {
          std::cerr << "Failed to read value for argument for -s/--serialize_snapshots. Value must be a boolean." << std::endl;
          return false;
        }
        g_overwrite_snapshots = !g_overwrite_snapshots;
        break;
      }
      case 'r':
      {
        serialized_file_name = optarg;
        break;
      }
      case 'e':
      {
        if( !StringUtilities::extractFromString( optarg, end_time_override ) )
        {
          std::cerr << "Failed to read value for argument for -e/--end. Value must be a positive scalar." << std::endl;
          return false;
        }
        if( end_time_override <= 0 )
        {
          std::cerr << "Failed to read value for argument for -e/--end. Value must be a positive scalar." << std::endl;
          return false;
        }
        break;
      }
      #ifdef USE_HDF5
      case 'o':
      {
        g_output_dir_name = optarg;
        break;
      }
      #endif
      case 'f':
      {
        if( !StringUtilities::extractFromString( optarg, output_frequency ) )
        {
          std::cerr << "Failed to read value for argument for -f/--frequency. Value must be an unsigned integer." << std::endl;
          return false;
        }
        break;
      }
      case '?':
      {
        return false;
      }
      default:
      {
        std::cerr << "This is a bug in the command line parser. Please file a report." << std::endl;
        return false;
      }
    }
  }

  return true;
}

static std::string basisFunctionTypeToString( const BasisFunctionType& bft )
{
  switch( bft )
  {
    case BasisFunctionType::Linear:
      return "Standard 1st order";
    case BasisFunctionType::ThirdOrder:
      return "Stardard 3rd order";
    case BasisFunctionType::uGIMPLinear:
      return "uGIMP 1st order";
    #ifndef CMAKE_DETECTED_CLANG_COMPILER
    default:
    {
      std::cerr << "Invalid basis function order in basisFunctionTypeToString, this is a bug." << std::endl;
      std::exit( EXIT_FAILURE );
    }
    #endif
  }
}

int main( int argc, char** argv )
{
  // Command line options
  bool help_mode_enabled{ false };
  scalar end_time_override{ -1.0 };
  unsigned output_frequency{ 0 };
  std::string serialized_file_name;

  g_serialize_snapshots = false;
  g_overwrite_snapshots = true;

  // Attempt to load command line options
  if( !parseCommandLineOptions( &argc, &argv, help_mode_enabled, end_time_override, output_frequency, serialized_file_name ) )
  {
    return EXIT_FAILURE;
  }

  // If the user requested help, print help and exit
  if( help_mode_enabled )
  {
    printUsage( argv[0] );
    return EXIT_SUCCESS;
  }

  if( !serialized_file_name.empty() )
  {
    if( deserializeSystem( serialized_file_name ) == EXIT_FAILURE )
    {
      return EXIT_FAILURE;
    }
    return executeSimLoop();
  }

  // The user must provide the path to an xml scene file
  if( argc != optind + 1 )
  {
    printUsage( argv[0] );
    return EXIT_FAILURE;
  }

  // Attempt to load the user-provided scene
  if( !loadXMLScene( std::string{ argv[optind] } ) )
  {
    return EXIT_FAILURE;
  }

  // Override the default end time with the requested one, if provided
  if( end_time_override > 0.0 )
  {
    g_end_time = end_time_override;
  }

  // Compute the data output rate
  assert( g_dt.positive() );
  // If the user provided an output frequency
  if( output_frequency != 0 )
  {
    const Rational<std::intmax_t> potential_steps_per_frame{ std::intmax_t( 1 ) / ( g_dt * std::intmax_t( output_frequency ) ) };
    if( !potential_steps_per_frame.isInteger() )
    {
      std::cerr << "Timestep and output frequency do not yield an integer number of timesteps for data output. Exiting." << std::endl;
      return EXIT_FAILURE;
    }
    g_steps_per_save = unsigned( potential_steps_per_frame.numerator() );
  }
  // Otherwise default to dumping every frame
  else
  {
    g_steps_per_save = 1;
  }
  g_iteration = 0;
  g_output_frame = 0;
  assert( g_end_time > 0.0 );
  g_save_number_width = MathUtilities::computeNumDigits( 1 + unsigned( ceil( g_end_time / scalar( g_dt ) ) ) / g_steps_per_save );

  printCompileInfo( std::cout );

  std::cout << "Simulation stats:" << std::endl;
  std::cout << "  Material Point Count: " << g_state.material_points.npoints << std::endl;
  std::cout << "  Basis function order: " << basisFunctionTypeToString( g_state.basis_functions->type() ) << std::endl;
  std::cout << "  Shear Modulus: " << g_state.shear_modulus << std::endl;
  std::cout << "  Bulk Modulus: " << g_state.bulk_modulus << std::endl;
  std::cout << "  Drucker-Prager Alpha: " << g_state.Drucker_Prager_alpha << std::endl;
  std::cout << "  Gravity: " << g_state.near_earth_gravity.transpose() << std::endl;
  std::cout << "  PIC/FLIP Alpha: " << g_state.alpha << std::endl;
  std::cout << "  Time-Step: " << g_dt << std::endl;
  std::cout << "  End Time: " << g_end_time << std::endl;

  return executeSimLoop();
}
