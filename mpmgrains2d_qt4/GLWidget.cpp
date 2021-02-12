#include "GLWidget.h"

#include <QtGui>
#include <QtOpenGL>

#include <iostream>
#include <limits>

#include "mpmgrains2dutils/MPMGrains2DParser.h"
#include "mpmgrains2d/InitialSimulationState.h"
#include "mpmgrains2d/ExplicitIntegrator.h"
#include "mpmgrains2d/StaticPlane.h"

#include "scisim/HDF5File.h"
#include "scisim/CompileDefinitions.h"

#ifndef NDEBUG
static std::string glErrorToString( const GLenum error_code )
{
  switch( error_code )
  {
    case GL_NO_ERROR:
      return "GL_NO_ERROR";
    case GL_INVALID_ENUM:
      return "GL_INVALID_ENUM";
    case GL_INVALID_VALUE:
      return "GL_INVALID_VALUE";
    case GL_INVALID_OPERATION:
      return "GL_INVALID_OPERATION";
    case GL_INVALID_FRAMEBUFFER_OPERATION:
      return "GL_INVALID_FRAMEBUFFER_OPERATION";
    case GL_OUT_OF_MEMORY:
      return "GL_OUT_OF_MEMORY";
    case GL_STACK_UNDERFLOW:
      return "GL_STACK_UNDERFLOW";
    case GL_STACK_OVERFLOW:
      return "GL_STACK_OVERFLOW";
    default:
      return "Unknown error. Please contact the maintainer of this code.";
  }
}

static bool checkGLErrors()
{
  const GLenum error_code = glGetError();
  if( error_code != GL_NO_ERROR )
  {
    std::cerr << "OpenGL error: " << glErrorToString( error_code ) << std::endl;
    return false;
  }
  return true;
}
#endif

static void getViewportDimensions( GLint& width, GLint& height )
{
  GLint viewport[4];
  glGetIntegerv( GL_VIEWPORT, viewport );
  width = viewport[2];
  height = viewport[3];
}

GLWidget::GLWidget( QWidget* parent )
: QGLWidget( QGLFormat( QGL::SampleBuffers ), parent )
, m_camera_controller()
, m_render_at_fps( false )
, m_lock_camera( false )
, m_last_pos()
, m_left_mouse_button_pressed( false )
, m_right_mouse_button_pressed( false )
, m_circle_renderer( 32 )
, m_display_precision( 9 )
, m_display_HUD( true )
, m_movie_dir_name()
, m_movie_dir()
, m_output_frame( 0 )
, m_output_fps()
, m_steps_per_frame()
, m_iteration()
, m_dt()
, m_end_time()
, m_state()
{}

GLWidget::~GLWidget()
{}

QSize GLWidget::minimumSizeHint() const
{
  return QSize( 50, 50 );
}

QSize GLWidget::sizeHint() const
{
  return QSize( int(m_camera_controller.width()), int(m_camera_controller.height()) );
}

bool GLWidget::openScene( const QString& xml_scene_file_name, unsigned& fps, bool& render_at_fps, bool& lock_camera )
{
  InitialSimulationState initial_state;
  const bool settings_loaded = MPMGrains2DParser::readXMLFile( xml_scene_file_name.toStdString(), initial_state );
  if( !settings_loaded )
  {
    std::cerr << "Failed to load file: " << xml_scene_file_name.toStdString() << std::endl;
    return false;
  }
  
  // Current timestep
  m_dt = initial_state.dt;
  // End time of the simulation
  m_end_time = initial_state.end_time;
    
  m_state = initial_state.generateSimulationState();
  std::cout << "Total MPM mass: " << m_state.material_points.totalMass() << std::endl;

  m_iteration = 0;

  /*
  if( input_settings.m_camera_settings.isSet() )
  {
    m_camera_controller.setScaleFactor( input_settings.m_camera_settings.scale() );
    m_camera_controller.setCenter( input_settings.m_camera_settings.center().x(), input_settings.m_camera_settings.center().y() );
    GLint width;
    GLint height;
    getViewportDimensions( width, height );
    m_camera_controller.reshape( unsigned(width), unsigned(height) );
  }
  else
  {
   //*/
    centerCamera( false );
  //}
  
  // Set the camera
  {
    const bool lock_backup{ m_lock_camera };
    m_lock_camera = false;
    
    /*
    // Set the camera
    if( camera_settings.set )
    {
      m_camera_controller.setCenter( camera_settings.center.x(), camera_settings.center.y() );
      m_camera_controller.setScaleFactor( camera_settings.scale );
      if( render_on_load )
      {
        updateGL();
      }
    }
    else
    {
     //*/
      centerCamera( false );
      // updateGL is embedded in centerCamera
    //}
    
    m_lock_camera = lock_backup;
  }


  assert( m_output_fps > 0 );
  setMovieFPS( m_output_fps );

  updateGL();

  // For the parent to update the UI
  fps = m_output_fps;
  render_at_fps = m_render_at_fps;
  lock_camera = m_lock_camera;

  return true;
}

void GLWidget::stepSystem()
{
  if( m_iteration * scalar( m_dt ) >= m_end_time )
  {
    std::cout << "Simulation complete at time " << m_iteration * scalar( m_dt );
    std::cout << " of " << m_end_time << ". Exiting." << std::endl;
    std::exit( EXIT_SUCCESS );
  }

  //const unsigned next_iter{ m_iteration + 1 };

  ExplicitIntegrator::flow( scalar(m_dt), m_state );

  ++m_iteration;

  if( !m_render_at_fps || m_iteration % m_steps_per_frame == 0 )
  {
    updateGL();
  }

  if( m_movie_dir_name.size() != 0 )
  {
    assert( m_steps_per_frame > 0 );
    if( m_iteration % m_steps_per_frame == 0 )
    {
      // Save a screenshot of the current state
      QString output_image_name = QString( tr("frame%1.png") ).arg( m_output_frame, 10, 10, QLatin1Char('0') );
      saveScreenshot( m_movie_dir.filePath( output_image_name ) );
      #ifdef USE_HDF5
      QString output_hdf5_name{ QString{ tr( "frame%1.h5" ) }.arg( m_output_frame, 10, 10, QLatin1Char('0') ) };
      saveHDF5File( m_movie_dir.filePath( output_hdf5_name ) );
      #endif
      ++m_output_frame;
    }
  }
}

void GLWidget::resetSystem()
{
  std::cerr << "GLWidget::resetSystem not implemented" << std::endl;
}

void GLWidget::initializeGL()
{
  qglClearColor( QColor( 255, 255, 255, 255 ) );
  assert( checkGLErrors() );
}

void GLWidget::resizeGL( int width, int height )
{
  assert( width >= 0 ); assert( height >= 0 );

  m_camera_controller.reshape( unsigned(width), unsigned(height) );

  assert( checkGLErrors() );
}

void GLWidget::paintGL()
{
  glMatrixMode( GL_MODELVIEW );

  glClear( GL_COLOR_BUFFER_BIT );

  if( axesDrawingIsEnabled() ) { paintAxes(); }

  paintSystem();

  if( m_display_HUD ) { paintHUD(); }

  assert( autoBufferSwap() );

  assert( checkGLErrors() );
}

bool GLWidget::axesDrawingIsEnabled() const
{
  return m_left_mouse_button_pressed;
}

void GLWidget::paintAxes() const
{
  // Draw the positive x axis
  qglColor( QColor( 255, 0, 0 ) );
  glLineWidth( 2.0 );
  glBegin( GL_LINES );
  glVertex4f( 0.0, 0.0, 0.0, 1.0 );
  glVertex4f( 1.0, 0.0, 0.0, 0.0 );
  glEnd();

  // Draw the negative x axis
  qglColor( QColor( 255, 0, 0 ) );
  glLineWidth( 2.0 );
  glLineStipple( 8, 0xAAAA );
  glEnable( GL_LINE_STIPPLE );
  glBegin( GL_LINES );
  glVertex4f( 0.0, 0.0, 0.0, 1.0 );
  glVertex4f( -1.0, 0.0, 0.0, 0.0 );
  glEnd();
  glDisable( GL_LINE_STIPPLE );

  // Draw the positive y axis
  qglColor( QColor( 0, 255, 0 ) );
  glLineWidth( 2.0 );
  glBegin( GL_LINES );
  glVertex4f( 0.0, 0.0, 0.0, 1.0 );
  glVertex4f( 0.0, 1.0, 0.0, 0.0 );
  glEnd();

  // Draw the negative y axis
  qglColor( QColor( 0, 255, 0 ) );
  glLineWidth( 2.0 );
  glLineStipple( 8, 0xAAAA );
  glEnable( GL_LINE_STIPPLE );
  glBegin( GL_LINES );
  glVertex4f( 0.0, 0.0, 0.0, 1.0 );
  glVertex4f( 0.0, -1.0, 0.0, 0.0 );
  glEnd();
  glDisable( GL_LINE_STIPPLE );
}

void GLWidget::renderAtFPS( const bool render_at_fps )
{
  m_render_at_fps = render_at_fps;
}

void GLWidget::lockCamera( const bool lock_camera )
{
  m_lock_camera = lock_camera;
}

void GLWidget::toggleHUD()
{
  m_display_HUD = !m_display_HUD;

  updateGL();
}

void GLWidget::centerCamera( const bool update_gl )
{
  if( m_lock_camera ) { return; }

  Eigen::Vector4d bbox;
  m_state.computeBoundingBox( bbox );
  const double minx = bbox( 0 );
  const double maxx = bbox( 1 );
  const double miny = bbox( 2 );
  const double maxy = bbox( 3 );

  const double cx = minx + 0.5 * ( maxx - minx );
  const double rx = maxx - cx;
  const double cy = miny + 0.5 * ( maxy - miny );
  const double ry = maxy - cy;

  const scalar ratio{ scalar( m_camera_controller.height() ) / scalar( m_camera_controller.width() ) };
  const scalar size{ 1.2 * std::max( ratio * rx, ry ) };
  
  m_camera_controller.setCenter( cx, cy );
  m_camera_controller.setScaleFactor( size );

  if( update_gl )
  {
    GLint width;
    GLint height;
    getViewportDimensions( width, height );
    m_camera_controller.reshape( width, height );
    updateGL();
  }
}

void GLWidget::saveScreenshot( const QString& file_name )
{
  std::cout << "Saving screenshot of time " << m_iteration * scalar( m_dt ) << " to " << file_name.toStdString() << std::endl;
  const QImage frame_buffer = grabFrameBuffer();
  frame_buffer.save( file_name );
}

#ifdef USE_HDF5
void GLWidget::saveHDF5File( const QString& file_name )
{
  std::cout << "Saving H5 state of time " << m_iteration * scalar( m_dt ) << " to " << file_name.toStdString() << std::endl;

  // Save the simulation state
  try
  {
    HDF5File output_file{ file_name.toStdString(), HDF5AccessType::READ_WRITE };
    // Save the iteration and time step and time
    output_file.writeScalar( "", "timestep", scalar( m_dt ) );
    output_file.writeScalar( "", "iteration", m_iteration );
    output_file.writeScalar( "", "time", scalar( m_dt ) * m_iteration );
    // Save out the git hash
    output_file.writeString( "", "git_hash", CompileDefinitions::GitSHA1 );
    // Write out the simulation data
    m_state.writeBinaryState( "", output_file );
  }
  catch( const std::string& error )
  {
    std::cerr << error << std::endl;
    std::exit( EXIT_FAILURE );
  }
}
#endif

void GLWidget::setMovieDir( const QString& dir_name )
{
  m_movie_dir_name = dir_name;
  m_output_frame = 0;

  // Save a screenshot of the current state
  if( m_movie_dir_name.size() != 0 )
  {
    m_movie_dir.setPath( m_movie_dir_name );
    assert( m_movie_dir.exists() );

    QString output_image_name = QString( tr("frame%1.png") ).arg( m_output_frame, 10, 10, QLatin1Char('0') );
    saveScreenshot( m_movie_dir.filePath( output_image_name ) );
    #ifdef USE_HDF5
    QString output_hdf5_name{ QString{ tr( "frame%1.h5" ) }.arg( m_output_frame, 10, 10, QLatin1Char('0') ) };
    saveHDF5File( m_movie_dir.filePath( output_hdf5_name ) );
    #endif
    ++m_output_frame;
  }
}

void GLWidget::setMovieFPS( const unsigned fps )
{
  assert( fps > 0 );
  m_output_fps = fps;
  m_output_frame = 0;
  if( 1.0 < scalar( m_dt * std::intmax_t( m_output_fps ) ) )
  {
    std::cerr << "Warning, requested movie frame rate faster than timestep. Dumping at timestep rate." << std::endl;
    m_steps_per_frame = 1;
  }
  else
  {
    const Rational<std::intmax_t> potential_steps_per_frame = std::intmax_t( 1 ) / ( m_dt * std::intmax_t( m_output_fps ) );
    if( !potential_steps_per_frame.isInteger() )
    {
      if( m_dt != Rational<std::intmax_t>( 0 ) )
      {
        std::cerr << "Warning, timestep and output frequency do not yield an integer number of timesteps for data output. Dumping at timestep rate." << std::endl;
      }
      m_steps_per_frame = 1;
    }
    else
    {
      m_steps_per_frame = unsigned( potential_steps_per_frame.numerator() );
    }
  }
}

void GLWidget::exportCameraSettings()
{
  std::cout << "<camera center=\"" << m_camera_controller.centerX() << " " << m_camera_controller.centerY() << "\" scale=\"" << m_camera_controller.scaleFactor() << "\" fps=\"" << m_output_fps << "\" render_at_fps=\"" << m_render_at_fps << "\" locked=\"" << m_lock_camera << "\"/>" << std::endl;
}

//static void paintInfiniteLine( const Vector2s& x, const Vector2s& n )
//{
//  const scalar theta = - scalar( 180.0 ) * atan2( n.x(), n.y() ) / MathUtilities::PI<scalar>();
//
//  glPushMatrix();
//
//  glTranslated( GLdouble( x.x() ), GLdouble( x.y() ), GLdouble( 0.0 ) );
//  glRotated( GLdouble( theta ), GLdouble( 0.0 ), GLdouble( 0.0 ), GLdouble( 1.0 ) );
//
//  glBegin( GL_LINES );
//  glVertex4d(  0.0,  0.0, 0.0, 1.0 );
//  glVertex4d(  1.0,  0.0, 0.0, 0.0 );
//  glVertex4d(  0.0,  0.0, 0.0, 1.0 );
//  glVertex4d( -1.0,  0.0, 0.0, 0.0 );
//  glEnd();
//
//  glPopMatrix();
//}

//static void paintSolidHalfPlane( const Vector2s& x, const Vector2s& n )
//{
//  paintInfiniteLine( x, n );
//
//  //const scalar theta = -180.0 * atan2( n.x(), n.y() ) / 3.14159265359;
//  //
//  //glPushMatrix();
//  //
//  //glTranslated( (GLdouble) x.x(), (GLdouble) x.y(), (GLdouble) 0.0 );
//  //glRotated( (GLdouble) theta, 0.0, (GLdouble) 0.0, (GLdouble) 1.0 );
//  //
//  //glBegin( GL_TRIANGLES );
//  //glVertex4d(  0.0,  0.0, 0.0, 1.0 );
//  //glVertex4d(  1.0,  0.0, 0.0, 0.0 );
//  //glVertex4d(  0.0, -1.0, 0.0, 0.0 );
//  //glVertex4d(  0.0,  0.0, 0.0, 1.0 );
//  //glVertex4d( -1.0,  0.0, 0.0, 0.0 );
//  //glVertex4d(  0.0, -1.0, 0.0, 0.0 );
//  //glEnd();
//  //
//  //glPopMatrix();
//}

// TODO: Abstract the shared code in here into its own function
//static void paintPlanarPortal( const PlanarPortal& planar_portal )
//{
//  // Draw the first plane of the portal
//  {
//    const scalar theta = -180.0 * atan2( planar_portal.planeA().n().x(), planar_portal.planeA().n().y() ) / MathDefines::PI<scalar>();
//
//    glPushMatrix();
//    glTranslated( GLdouble( planar_portal.planeA().x().x() ), GLdouble( planar_portal.planeA().x().y() ), 0.0 );
//    glRotated( GLdouble( theta ), 0.0, 0.0, 1.0 );
//
//    glLineStipple( 8, 0xAAAA );
//    glEnable( GL_LINE_STIPPLE );
//    glBegin( GL_LINES );
//    glVertex4d( 0.0, 0.0, 0.0, 1.0 );
//    glVertex4d( -1.0, 0.0, 0.0, 0.0 );
//    glEnd();
//    glDisable( GL_LINE_STIPPLE );
//
//    glLineStipple( 8, 0x5555 );
//    glEnable( GL_LINE_STIPPLE );
//    glBegin( GL_LINES );
//    glVertex4d( 0.0, 0.0, 0.0, 1.0 );
//    glVertex4d( 1.0, 0.0, 0.0, 0.0 );
//    glEnd();
//    glDisable( GL_LINE_STIPPLE );
//
//    glPopMatrix();
//
//    assert( planar_portal.boundsA()(0) < planar_portal.boundsA()(1) );
//    assert( ( planar_portal.boundsA()(0) != -SCALAR_INFINITY && planar_portal.boundsA()(1) != SCALAR_INFINITY ) || ( planar_portal.boundsA()(0) == -SCALAR_INFINITY && planar_portal.boundsA()(1) == SCALAR_INFINITY ) );
//    if( planar_portal.boundsA()(0) != -SCALAR_INFINITY )
//    {
//      // Draw the lower bound
//      glPushMatrix();
//      glTranslated( GLdouble( planar_portal.planeA().x().x() ), GLdouble( planar_portal.planeA().x().y() ), 0.0 );
//      glRotated( GLdouble( theta ), 0.0, 0.0, 1.0 );
//      glTranslated( planar_portal.boundsA()(0), 0.0, 0.0 );
//      glLineStipple( 8, 0x5555 );
//      glEnable( GL_LINE_STIPPLE );
//      glBegin( GL_LINES );
//      glVertex4d( 0.0, 0.0, 0.0, 1.0 );
//      glVertex4d( 0.0, -1.0, 0.0, 0.0 );
//      glEnd();
//      glDisable( GL_LINE_STIPPLE );
//      glPopMatrix();
//
//      // Draw the upper bound
//      glPushMatrix();
//      glTranslated( GLdouble( planar_portal.planeA().x().x() ), GLdouble( planar_portal.planeA().x().y() ), 0.0 );
//      glRotated( GLdouble( theta ), 0.0, 0.0, 1.0 );
//      glTranslated( planar_portal.boundsA()(1), 0.0, 0.0 );
//      glLineStipple( 8, 0x5555 );
//      glEnable( GL_LINE_STIPPLE );
//      glBegin( GL_LINES );
//      glVertex4d( 0.0, 0.0, 0.0, 1.0 );
//      glVertex4d( 0.0, -1.0, 0.0, 0.0 );
//      glEnd();
//      glDisable( GL_LINE_STIPPLE );
//      glPopMatrix();
//
//      // Draw a short line to indicate the current position of the center of the portal
//      glPushMatrix();
//      glTranslated( GLdouble( planar_portal.planeAx() ), GLdouble( planar_portal.planeAy() ), 0.0 );
//      glRotated( GLdouble( theta ), 0.0, 0.0, 1.0 );
//
//      // Draw an infinite line to show what half of portal is free
//      glLineStipple( 8, 0x5555 );
//      glEnable( GL_LINE_STIPPLE );
//      glBegin( GL_LINES );
//      glVertex2d( 0.0, 0.0 );
//      glVertex2d( 0.0, - 0.1 * ( planar_portal.boundsA()(1) - planar_portal.boundsA()(0) ) );
//      glEnd();
//      glDisable( GL_LINE_STIPPLE );
//
//      glPopMatrix();
//    }
//    else
//    {
//      glPushMatrix();
//      glTranslated( GLdouble( planar_portal.planeA().x().x() ), GLdouble( planar_portal.planeA().x().y() ), 0.0 );
//      glRotated( GLdouble( theta ), 0.0, 0.0, 1.0 );
//
//      // Draw an infinite line to show what half of portal is free
//      glLineStipple( 8, 0x5555 );
//      glEnable( GL_LINE_STIPPLE );
//      glBegin( GL_LINES );
//      glVertex4d( 0.0, 0.0, 0.0, 1.0 );
//      glVertex4d( 0.0, -1.0, 0.0, 0.0 );
//      glEnd();
//      glDisable( GL_LINE_STIPPLE );
//
//      glPopMatrix();
//    }
//  }
//
//  // Draw the second plane of the portal
//  {
//    const scalar theta = -180.0 * atan2( planar_portal.planeB().n().x(), planar_portal.planeB().n().y() ) / MathDefines::PI<scalar>();
//
//    glPushMatrix();
//    glTranslated( GLdouble( planar_portal.planeB().x().x() ), GLdouble( planar_portal.planeB().x().y() ), 0.0 );
//    glRotated( GLdouble( theta ), 0.0, 0.0, 1.0 );
//
//    glLineStipple( 8, 0xAAAA );
//    glEnable( GL_LINE_STIPPLE );
//    glBegin( GL_LINES );
//    glVertex4d( 0.0, 0.0, 0.0, 1.0 );
//    glVertex4d( -1.0, 0.0, 0.0, 0.0 );
//    glEnd();
//    glDisable( GL_LINE_STIPPLE );
//
//    glLineStipple( 8, 0x5555 );
//    glEnable( GL_LINE_STIPPLE );
//    glBegin( GL_LINES );
//    glVertex4d( 0.0, 0.0, 0.0, 1.0 );
//    glVertex4d( 1.0, 0.0, 0.0, 0.0 );
//    glEnd();
//    glDisable( GL_LINE_STIPPLE );
//
//    glPopMatrix();
//
//    assert( planar_portal.boundsB()(0) < planar_portal.boundsB()(1) );
//    assert( ( planar_portal.boundsB()(0) != -SCALAR_INFINITY && planar_portal.boundsB()(1) != SCALAR_INFINITY ) || ( planar_portal.boundsB()(0) == -SCALAR_INFINITY && planar_portal.boundsB()(1) == SCALAR_INFINITY ) );
//    if( planar_portal.boundsA()(0) != -SCALAR_INFINITY )
//    {
//      // Draw the lower bound
//      glPushMatrix();
//      glTranslated( GLdouble( planar_portal.planeB().x().x() ), GLdouble( planar_portal.planeB().x().y() ), 0.0 );
//      glRotated( GLdouble( theta ), 0.0, 0.0, 1.0 );
//      glTranslated( planar_portal.boundsB()(0), 0.0, 0.0 );
//      glLineStipple( 8, 0x5555 );
//      glEnable( GL_LINE_STIPPLE );
//      glBegin( GL_LINES );
//      glVertex4d( 0.0, 0.0, 0.0, 1.0 );
//      glVertex4d( 0.0, -1.0, 0.0, 0.0 );
//      glEnd();
//      glDisable( GL_LINE_STIPPLE );
//      glPopMatrix();
//
//      // Draw the upper bound
//      glPushMatrix();
//      glTranslated( GLdouble( planar_portal.planeB().x().x() ), GLdouble( planar_portal.planeB().x().y() ), 0.0 );
//      glRotated( GLdouble( theta ), 0.0, 0.0, 1.0 );
//      glTranslated( planar_portal.boundsB()(1), 0.0, 0.0 );
//      glLineStipple( 8, 0x5555 );
//      glEnable( GL_LINE_STIPPLE );
//      glBegin( GL_LINES );
//      glVertex4d( 0.0, 0.0, 0.0, 1.0 );
//      glVertex4d( 0.0, -1.0, 0.0, 0.0 );
//      glEnd();
//      glDisable( GL_LINE_STIPPLE );
//      glPopMatrix();
//
//      // Draw a short line to indicate the current position of the center of the portal
//      glPushMatrix();
//      glTranslated( GLdouble( planar_portal.planeBx() ), GLdouble( planar_portal.planeBy() ), 0.0 );
//      glRotated( GLdouble( theta ), 0.0, 0.0, 1.0 );
//
//      // Draw an infinite line to show what half of portal is free
//      glLineStipple( 8, 0x5555 );
//      glEnable( GL_LINE_STIPPLE );
//      glBegin( GL_LINES );
//      glVertex2d( 0.0, 0.0 );
//      glVertex2d( 0.0, - 0.1 * ( planar_portal.boundsB()(1) - planar_portal.boundsB()(0) ) );
//      glEnd();
//      glDisable( GL_LINE_STIPPLE );
//
//      glPopMatrix();
//    }
//    else
//    {
//      glPushMatrix();
//      glTranslated( GLdouble( planar_portal.planeB().x().x() ), GLdouble( planar_portal.planeB().x().y() ), 0.0 );
//      glRotated( GLdouble( theta ), 0.0, 0.0, 1.0 );
//
//      // Draw an infinite line to show what half of portal is free
//      glLineStipple( 8, 0x5555 );
//      glEnable( GL_LINE_STIPPLE );
//      glBegin( GL_LINES );
//      glVertex4d( 0.0, 0.0, 0.0, 1.0 );
//      glVertex4d( 0.0, -1.0, 0.0, 0.0 );
//      glEnd();
//      glDisable( GL_LINE_STIPPLE );
//
//      glPopMatrix();
//    }
//  }
//}

void GLWidget::drawPoints() const
{
  glPushAttrib( GL_COLOR );
  glColor3d( 0.0, 0.8, 0.0 );
  for( int i = 0; i < int(m_state.material_points.npoints); ++i )
  {
    glPushMatrix();
    glTranslated( m_state.material_points.q( 0, i ), m_state.material_points.q( 1, i ), 0.0 );
    glScaled( 0.05, 0.05, 1.0 );
    m_circle_renderer.renderCircle();
    glPopMatrix();
  }
  glPopAttrib();
}

static void paintInfiniteLine( const Vector2s& x, const Vector2s& n )
{
  const scalar theta{ - scalar( 180.0 ) * atan2( n.x(), n.y() ) / MathDefines::PI<scalar>() };

  glPushMatrix();

  glTranslated( GLdouble( x.x() ), GLdouble( x.y() ), GLdouble( 0.0 ) );
  glRotated( GLdouble( theta ), GLdouble( 0.0 ), GLdouble( 0.0 ), GLdouble( 1.0 ) );

  glBegin( GL_LINES );
  glVertex4d(  0.0,  0.0, 0.0, 1.0 );
  glVertex4d(  1.0,  0.0, 0.0, 0.0 );
  glVertex4d(  0.0,  0.0, 0.0, 1.0 );
  glVertex4d( -1.0,  0.0, 0.0, 0.0 );
  glEnd();

  glPopMatrix();
}

static void drawPlanes( const std::vector<MPMStaticPlane>& static_planes )
{
  glPushAttrib( GL_COLOR );
  glPushAttrib( GL_LINE_WIDTH );
  glColor3d( 0.0, 0.0, 0.0 );
  glLineWidth( 2.0 );
  for( const MPMStaticPlane& plane : static_planes )
  {
    paintInfiniteLine( plane.x, plane.n );
  }
  glPopAttrib();
  glPopAttrib();
}

static void drawGrid( const PhysicsGrid& grid )
{
  if( ( grid.cell_count.array() == 0 ).all() )
  {
    return;
  }

  glPushAttrib( GL_LINE_WIDTH );
  glPushAttrib( GL_COLOR );

  glLineWidth( 0.5 );
  glColor3d( 0, 0, 0 );

  glBegin( GL_LINES );
  for( int vertical_line_num = 0; vertical_line_num < int(grid.cell_count.x()) + 1; ++vertical_line_num )
  {
    const Vector2s x0 = Vector2s(grid.min) + vertical_line_num * grid.cell_width * Vector2s::UnitX();
    const Vector2s x1 = Vector2s(grid.min(0), grid.max(1)) + vertical_line_num * grid.cell_width * Vector2s::UnitX();
    glVertex2d( x0.x(), x0.y() );
    glVertex2d( x1.x(), x1.y() );
  }
  for( int horizontal_line_num = 0; horizontal_line_num < int(grid.cell_count.y()) + 1; ++horizontal_line_num )
  {
    const Vector2s x0 = Vector2s(grid.min) + horizontal_line_num * grid.cell_width * Vector2s::UnitY();
    const Vector2s x1 = Vector2s(grid.max(0), grid.min(1)) + horizontal_line_num * grid.cell_width * Vector2s::UnitY();
    glVertex2d( x0.x(), x0.y() );
    glVertex2d( x1.x(), x1.y() );
  }
  glEnd();

  glPopAttrib();
  glPopAttrib();
}

//static void paintLine( const Vector2s& p0, const Vector2s& p1 )
//{
//  glBegin( GL_LINES );
//  glVertex2d( p0.x(), p0.y() );
//  glVertex2d( p1.x(), p1.y() );
//  glEnd();
//}

void GLWidget::paintSystem() const
{
  drawPoints();
  drawGrid( m_state.physics_grid );
  drawPlanes( m_state.static_planes );

  /*
  if( m_simulation.periodicBoundaryEnabled() )
  {
    const double xmin{ m_simulation.grid().simulationRegion()(0) };
    const double xmax{ m_simulation.grid().simulationRegion()(2) };
    const double ymin{ m_simulation.grid().simulationRegion()(1) };
    const double ymax{ m_simulation.grid().simulationRegion()(3) };
    glPushAttrib( GL_LINE_WIDTH );
    glPushAttrib( GL_COLOR );
    glLineWidth( 2.0 );
    glColor3d( 0.0, 0.0, 0.0 );
    glEnable( GL_LINE_STIPPLE );
    glLineStipple( 8, 0xAAAA );
    // Paint the left border
    paintLine( { xmin, ymin }, { xmin, ymax } );
    // Paint the right border
    paintLine( { xmax, ymin }, { xmax, ymax } );
    // Paint the top border
    paintLine( { xmin, ymax }, { xmax, ymax } );
    // Paint the bottom border
    paintLine( { xmin, ymin }, { xmax, ymin } );
    glDisable( GL_LINE_STIPPLE );
    // Paint an indicators to show the vertical border's location
    paintLine( { m_simulation.leesEdwardsPos(), ymax }, { m_simulation.leesEdwardsPos(), ymax + .03 * ( ymax - ymin ) } );
    paintLine( { -m_simulation.leesEdwardsPos(), ymin }, { -m_simulation.leesEdwardsPos(), ymin - .03 * ( ymax - ymin ) } );
    glPopAttrib();
    glPopAttrib();
  }
  //*/

//  m_renderer.draw( m_simulation, m_camera_controller.width(), m_camera_controller.height() );
  //for( std::vector<PlaneObstacle>::size_type i = 0; i < m_simulation.m_plane_obstacles.size(); ++i ) paintInfiniteLine( m_simulation.m_plane_obstacles[i].x(), m_simulation.m_plane_obstacles[i].n() );

  //const SimulationStateBalls2D& state = m_sim.state();
  //
  //// Draw each static drum
  //glPushAttrib( GL_COLOR );
  //qglColor( QColor( 0, 0, 0 ) );
  //{
  //  const std::vector<StaticDrum>& drums = state.staticDrums();
  //  for( std::vector<StaticDrum>::size_type i = 0; i < drums.size(); ++i )
  //  {
  //    m_circle_renderer.renderSolidCircle( drums[i].x(), drums[i].r() );
  //  }
  //}
  //glPopAttrib();

  //// Draw each planar portal
  //glPushAttrib( GL_COLOR );
  //glPushAttrib( GL_LINE_WIDTH );
  //glLineWidth( 2.0 );
  //{
  //  // TODO: Create a set number of nice looking colors for the portal ahead of time instead of regenerating them
  //  std::mt19937_64 mt( 123456 );
  //  std::uniform_int_distribution<int> color_gen( 0, 255 );
  //  const std::vector<PlanarPortal>& planar_portals = m_sim.state().planarPortals();
  //  for( const PlanarPortal& planar_portal : planar_portals )
  //  {
  //    const int r = color_gen( mt );
  //    const int g = color_gen( mt );
  //    const int b = color_gen( mt );
  //    qglColor( QColor( r, g, b ) );
  //    paintPlanarPortal( planar_portal );
  //  }
  //}
  //glPopAttrib();
  //glPopAttrib();

  //// Draw each static plane
  //glPushAttrib( GL_COLOR );
  //glPushAttrib( GL_LINE_WIDTH );
  //glLineWidth( 2.0 );
  //qglColor( QColor( 0, 0, 0 ) );
  //{
  //  const std::vector<RigidBody2DStaticPlane>& planes = m_sim.state().planes();
  //  for( const RigidBody2DStaticPlane& plane : planes )
  //  {
  //    paintInfiniteLine( plane.x(), plane.n() );
  //  }
  //}
  //glPopAttrib();
  //glPopAttrib();

  //// Draw each body
  //{
  //  const VectorXs& q = m_sim.state().q();
  //  assert( q.size() % 3 == 0 );
  //  const unsigned nbodies = q.size() / 3;
  //  for( unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx )
  //  {
  //    glPushMatrix();
  //    glTranslated( GLdouble( q( 3 * bdy_idx + 0 ) ), GLdouble( q( 3 * bdy_idx + 1 ) ), GLdouble( 0.0 ) );
  //    const scalar theta_degrees = 180.0 * q( 3 * bdy_idx + 2 ) / MathDefines::PI<scalar>();
  //    glRotated( theta_degrees, 0.0, 0.0, 1.0 );
  //    assert( bdy_idx < m_body_colors.size() / 3 );
  //    assert( bdy_idx < m_sim.state().geometryIndices().size() );
  //    assert( m_sim.state().geometryIndices()(bdy_idx) < m_body_renderers.size() );
  //    m_body_renderers[ m_sim.state().geometryIndices()(bdy_idx) ]->render( m_body_colors.segment<3>( 3 * bdy_idx ) );
  //    glPopMatrix();
  //  }
  //}
}

static QString generateTimeString( const unsigned iteration, const Rational<std::intmax_t>& dt, const int display_precision, const scalar& end_time )
{
  QString time_string( QObject::tr("t: ") );
  time_string += QString::number( iteration * scalar( dt ), 'f', display_precision );
  if( std::isinf( end_time ) )
  {
    time_string += QString( QObject::tr(" / ") );
    time_string += QString::number( end_time );
  }
  return time_string;
}

void GLWidget::paintHUD()
{
  static int text_width = 0;

  glPushAttrib( GL_MATRIX_MODE );
  glMatrixMode( GL_PROJECTION );
  glPushMatrix();
  glMatrixMode( GL_MODELVIEW );
  glPushMatrix();

  // Set an orthographic projection with height and width equal to window height and width
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  GLint width;
  GLint height;
  getViewportDimensions( width, height );
  glOrtho( 0, width, 0, height, -1, 1 );
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();

  // Enable blending for transparent HUD elements
  glPushAttrib( GL_BLEND );
  glEnable( GL_BLEND );
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

  // Draw a semi-transparent overlay so text is visible regardless of background color
  const Eigen::Matrix<GLdouble,2,1> overlay_start( 0, height - 1 * 12 - 2 );
  const Eigen::Matrix<GLdouble,2,1> overlay_extnt( text_width + 2 + 2, height );
  glColor4d( 0.0, 0.0, 0.0, 0.5 );
  glBegin( GL_QUADS );
  glVertex2d( GLdouble( overlay_start.x() ), GLdouble( overlay_start.y() ) );
  glVertex2d( GLdouble( overlay_start.x() + overlay_extnt.x() ), GLdouble( overlay_start.y() ) );
  glVertex2d( GLdouble( overlay_start.x() + overlay_extnt.x() ), GLdouble( overlay_start.y() + overlay_extnt.y() ) );
  glVertex2d( GLdouble( overlay_start.x() ), GLdouble( overlay_start.y() + overlay_extnt.y() ) );
  glEnd();

  glDisable( GL_BLEND );
  glPopAttrib();

	glMatrixMode( GL_MODELVIEW );
  glPopMatrix();
  glMatrixMode( GL_PROJECTION );
  glPopMatrix();
  glPopAttrib();

  // String to display in upper left corner
  const QString time_string{ generateTimeString( m_iteration, m_dt, m_display_precision, m_end_time ) };
  //const QString delta_H = generateNumericString( " dH: ", m_delta_H0 );
  //const QString delta_px = generateNumericString( "dpx: ", m_delta_p0.x() );
  //const QString delta_py = generateNumericString( "dpy: ", m_delta_p0.y() );
  //const QString delta_L = generateNumericString( " dL: ", m_delta_L0 );
  {
    const QFontMetrics font_metrics( QFont( "Courier", 12 ) );
    text_width = std::max( text_width, font_metrics.boundingRect( time_string ).width() );
    //text_width = std::max( text_width, font_metrics.boundingRect( delta_H ).width() );
    //text_width = std::max( text_width, font_metrics.boundingRect( delta_px ).width() );
    //text_width = std::max( text_width, font_metrics.boundingRect( delta_py ).width() );
    //text_width = std::max( text_width, font_metrics.boundingRect( delta_L ).width() );
  }

  qglColor( QColor( 255, 255, 255 ) );
  const QFont font( "Courier", 12 );
  renderText( 2, font.pointSize(), time_string, font );
  //renderText( 2, 2 * font.pointSize(), delta_H, font );
  //renderText( 2, 3 * font.pointSize(), delta_px, font );
  //renderText( 2, 4 * font.pointSize(), delta_py, font );
  //renderText( 2, 5 * font.pointSize(), delta_L, font );

  assert( checkGLErrors() );
}

void GLWidget::mousePressEvent( QMouseEvent* event )
{
  if( m_lock_camera ) { return; }
  
  bool repaint_needed = false;

  if( event->buttons() & Qt::LeftButton )
  {
    m_left_mouse_button_pressed = true;
    repaint_needed = true;
  }
  if( event->buttons() & Qt::RightButton )
  {
    m_right_mouse_button_pressed = true;
  }

  if( repaint_needed ) updateGL();

  m_last_pos = event->pos();
}

void GLWidget::mouseReleaseEvent( QMouseEvent* event )
{
  if( m_lock_camera ) { return; }

  bool repaint_needed = false;

  if( !( event->buttons() & Qt::LeftButton ) && m_left_mouse_button_pressed )
  {
    m_left_mouse_button_pressed = false;
    repaint_needed = true;
  }
  if( !( event->buttons() & Qt::RightButton ) && m_right_mouse_button_pressed )
  {
    m_right_mouse_button_pressed = false;
  }

  if( repaint_needed ) { updateGL(); }
}

void GLWidget::mouseMoveEvent( QMouseEvent* event )
{
  if( m_lock_camera ) { return; }

  const int dx = event->x() - m_last_pos.x();
  const int dy = event->y() - m_last_pos.y();
  m_last_pos = event->pos();

  bool repaint_needed = false;

  if( event->buttons() & Qt::LeftButton )
  {
    assert( m_left_mouse_button_pressed );
    m_camera_controller.translateView( dx, dy );
    repaint_needed = true;
  }

  if( event->buttons() & Qt::RightButton )
  {
    assert( m_right_mouse_button_pressed );
    m_camera_controller.zoomView( 0.02 * m_camera_controller.scaleFactor() * dy );
    repaint_needed = true;
  }

  if( repaint_needed ) { updateGL(); }

  assert( checkGLErrors() );
}

void GLWidget::wheelEvent( QWheelEvent* event )
{
  if( m_lock_camera ) { return; }

  m_camera_controller.zoomView( -0.002 * m_camera_controller.scaleFactor() * event->delta() );
  updateGL();
  assert( checkGLErrors() );
}
