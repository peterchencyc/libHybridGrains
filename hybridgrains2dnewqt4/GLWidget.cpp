#include "GLWidget.h"

#include <QtGui>
#include <iostream>

#include "rigidbody2d/BoxGeometry.h"
#include "rigidbody2d/CircleGeometry.h"
#include "rigidbody2d/RigidBody2DIntegratorSettings.h"

#include "rigidbody2dutils/CameraSettings2D.h"
#include "rigidbody2dutils/RigidBody2DSceneParser.h"

#include "hybridgrains2dnew/DiscreteIntegrator.h"
#include "hybridgrains2dnew/HybridDefinitions.h"

#include "hybridgrains2dnewutils/HybridGrains2DSceneParser.h"

#include "mpmgrains2dutils/MPMGrains2DParser.h"

#include "mpmgrains2d/InitialSimulationState.h"

#include "scisim/CompileDefinitions.h"
#include "scisim/HDF5File.h"

#ifndef NDEBUG
static std::string glErrorToString(const GLenum error_code) {
  switch (error_code) {
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

static bool checkGLErrors() {
  const GLenum error_code{glGetError()};
  if (error_code != GL_NO_ERROR) {
    std::cerr << "OpenGL error: " << glErrorToString(error_code) << std::endl;
    return false;
  }
  return true;
}
#endif

static void getViewportDimensions(GLint &width, GLint &height) {
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  width = viewport[2];
  height = viewport[3];
}

static std::vector<Vector3s> generateTemplateColors(const int ncolors) {
  std::vector<Vector3s> colors(ncolors);

  std::mt19937_64 color_gen_mt;
  color_gen_mt.seed(1337);
  {
    std::uniform_real_distribution<scalar> color_gen{0.0, 1.0};
    for (int color_num = 0; color_num < ncolors; color_num++) {
      scalar r = 1.0;
      scalar g = 1.0;
      scalar b = 1.0;
      while ((0.2126 * r + 0.7152 * g + 0.0722 * b) > 0.9 ||
             (0.2126 * r + 0.7152 * g + 0.0722 * b) < 0.4) {
        r = color_gen(color_gen_mt);
        g = color_gen(color_gen_mt);
        b = color_gen(color_gen_mt);
      }
      colors[color_num] << r, g, b;
    }
  }

  return colors;
}

GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(QGLFormat{QGL::SampleBuffers}, parent), m_camera_controller(),
      m_render_at_fps(false), m_lock_camera(false), m_last_pos(),
      m_left_mouse_button_pressed(false), m_right_mouse_button_pressed(false),
      m_circle_renderer(32), m_template_colors(generateTemplateColors(100)),
      m_display_precision(0), m_display_HUD(true), m_display_grid(true),
      m_display_mpm(true), m_display_discrete(true),
      m_display_zone_indicators(true), m_movie_dir_name(), m_movie_dir(),
      m_output_frame(0), m_output_fps(), m_steps_per_frame(), m_sim0(),
      m_sim() {}

GLWidget::~GLWidget() {}

QSize GLWidget::minimumSizeHint() const { return QSize{50, 50}; }

QSize GLWidget::sizeHint() const {
  return QSize{static_cast<int>(m_camera_controller.width()),
               static_cast<int>(m_camera_controller.height())};
}

static int computeTimestepDisplayPrecision(const Rational<std::intmax_t> &dt,
                                           const std::string &dt_string) {
  // Specified as a float
  if (dt_string.find('.') != std::string::npos) {
    return int(StringUtilities::computeNumCharactersToRight(dt_string, '.'));
  }
  // Specified as a rational
  else if (dt_string.find('/') != std::string::npos) {
    std::string converted_dt_string;
    std::stringstream ss;
    ss << std::fixed << scalar(dt);
    ss >> converted_dt_string;
    return int(
        StringUtilities::computeNumCharactersToRight(converted_dt_string, '.'));
  }
  // Specified as an int
  else {
    return int(dt_string.length());
  }
}

bool GLWidget::openScene(const QString &xml_scene_file_name,
                         const bool &render_on_load, unsigned &fps,
                         bool &render_at_fps, bool &lock_camera) {
  // Load the names of the discrete and continuum configuration files
  HybridIntegratorSettings integrator_settings;
  if (!HybridGrains2DSceneParser::parseXMLSceneFile(
          xml_scene_file_name.toStdString(), integrator_settings)) {
    return false;
  }

  // Compute the number of characters after the decimal point in the timestep
  // string
  m_display_precision = computeTimestepDisplayPrecision(
      integrator_settings.overall_dt, integrator_settings.time_step_string);

  // Load the discrete state
  std::cout << "Loading discrete config file: "
            << integrator_settings.discrete_file_name << std::endl;
  std::string dt_string;
  std::string scripting_callback_name;
  RigidBody2DState new_state;
  RigidBody2DIntegratorSettings new_integrator_settings;
  CameraSettings2D camera_settings;
  const bool discrete_loaded_successfully{
      RigidBody2DSceneParser::parseXMLSceneFile(
          integrator_settings.discrete_file_name, scripting_callback_name,
          new_state, new_integrator_settings, dt_string, camera_settings)};
  if (!discrete_loaded_successfully) {
    std::cerr << "Failed to load discrete configuration file: "
              << integrator_settings.discrete_file_name << std::endl;
    return false;
  }

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

  std::cout << "Discrete timestep: " << new_integrator_settings.dt << std::endl;
  std::cout << "Continuum timestep: " << initial_continuum_state.dt
            << std::endl;
  std::cout << "Overall timestep: " << integrator_settings.overall_dt
            << std::endl;

  if (new_integrator_settings.dt != integrator_settings.overall_dt) {
    std::cerr << "Error, discrete and overall timestep disagree." << std::endl;
    return false;
  }

  if (initial_continuum_state.dt != integrator_settings.overall_dt) {
    std::cerr << "Error, continuum and overall timestep disagree." << std::endl;
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
    m_sim = HybridGrains2DSim(
        HybridIntegratorState(
            integrator_settings.overall_dt, integrator_settings.end_time,
            discrete_integrator, mpm_integrator,
            integrator_settings.integrator_style,
            integrator_settings.kinematically_scripted_hybrid_fronts,
            integrator_settings.poorman_settings),
        new_state, initial_continuum_state.generateSimulationState());

    m_sim.initializeAvoidAVoid();

    // Initialize the discrete scripting callback
    {
      std::cout << "Initializing python callback..." << std::flush;
      std::string path;
      std::string file_name;
      StringUtilities::splitAtLastCharacterOccurence(
          xml_scene_file_name.toStdString(), path, file_name, '/');
      if (file_name.empty()) {
        using std::swap;
        swap(path, file_name);
      }
      m_sim.integratorState().discreteIntegrator().setPythonCallback(
          path, scripting_callback_name);
      m_sim.integratorState().discreteIntegrator().pythonStartOfSim(
          m_sim.discreteSim());
      std::cout << " done." << std::endl;
    }
  }

  m_sim0 = m_sim;

  // Set the camera
  {
    const bool lock_backup{m_lock_camera};
    m_lock_camera = false;

    // Set the camera
    if (camera_settings.set) {
      m_camera_controller.setCenter(camera_settings.center.x(),
                                    camera_settings.center.y());
      m_camera_controller.setScaleFactor(camera_settings.scale);
      if (render_on_load) {
        updateGL();
      }
    } else {
      centerCamera(render_on_load);
      // updateGL is embedded in centerCamera
    }

    m_lock_camera = lock_backup;
  }

  m_output_fps = camera_settings.fps;
  m_render_at_fps = camera_settings.render_at_fps;
  m_lock_camera = camera_settings.locked;
  setMovieFPS(m_output_fps);

  // For the parent to update the UI
  fps = m_output_fps;
  render_at_fps = m_render_at_fps;
  lock_camera = m_lock_camera;

  return true;
}

void GLWidget::stepSystem() {
  if (std::intmax_t(m_sim.integratorState().overallIteration()) *
          m_sim.integratorState().overallTimestep() >=
      m_sim.integratorState().endTime()) {
    // User-provided end of simulation python callback
    m_sim.integratorState().discreteIntegrator().pythonEndOfSim(
        m_sim.discreteSim());
    std::cout << "Simulation complete. Exiting." << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  m_sim.stepSystem();

  if (!m_render_at_fps ||
      m_sim.integratorState().overallIteration() % m_steps_per_frame == 0) {
    updateGL();
  }

  if (m_movie_dir_name.size() != 0) {
    assert(m_steps_per_frame > 0);
    if (m_sim.integratorState().overallIteration() % m_steps_per_frame == 0) {
      // Save a screenshot of the current state
      QString output_image_name{QString{tr("frame%1.png")}.arg(
          m_output_frame, 10, 10, QLatin1Char{'0'})};
      saveScreenshot(m_movie_dir.filePath(output_image_name));
#ifdef USE_HDF5
      QString output_hdf5_name{QString{tr("frame%1.h5")}.arg(
          m_output_frame, 10, 10, QLatin1Char('0'))};
      saveHDF5File(m_movie_dir.filePath(output_hdf5_name));
#endif
      m_output_frame++;
    }
  }
}

void GLWidget::resetSystem() {
  std::cout << "Resetting the simulation." << std::endl;

  m_sim = m_sim0;

  // m_H0 = m_sim.computeTotalEnergy();
  // m_p0 = m_sim.computeTotalMomentum();
  // m_L0 = m_sim.computeTotalAngularMomentum();
  // m_delta_H0 = 0.0;
  // m_delta_p0.setZero();
  // m_delta_L0 = 0.0;

  // Reset the output movie option
  m_movie_dir_name = QString{};
  m_movie_dir = QDir{};
  m_output_frame = 0;

  updateGL();
}

void GLWidget::initializeGL() {
  qglClearColor(QColor{255, 255, 255, 255});
  assert(checkGLErrors());
}

void GLWidget::resizeGL(int width, int height) {
  assert(width >= 0);
  assert(height >= 0);

  m_camera_controller.reshape(width, height);

  assert(checkGLErrors());
}

void GLWidget::paintGL() {
  glMatrixMode(GL_MODELVIEW);

  glClear(GL_COLOR_BUFFER_BIT);

  if (axesDrawingIsEnabled()) {
    paintAxes();
  }

  paintSystem();

  if (m_display_HUD) {
    paintHUD();
  }

  assert(autoBufferSwap());
  assert(checkGLErrors());
}

bool GLWidget::axesDrawingIsEnabled() const {
  return m_left_mouse_button_pressed;
}

void GLWidget::paintAxes() const {
  // Draw the positive x axis
  qglColor(QColor{255, 0, 0});
  glLineWidth(2.0);
  glBegin(GL_LINES);
  glVertex4f(0.0, 0.0, 0.0, 1.0);
  glVertex4f(1.0, 0.0, 0.0, 0.0);
  glEnd();

  // Draw the negative x axis
  qglColor(QColor{255, 0, 0});
  glLineWidth(2.0);
  glLineStipple(8, 0xAAAA);
  glEnable(GL_LINE_STIPPLE);
  glBegin(GL_LINES);
  glVertex4f(0.0, 0.0, 0.0, 1.0);
  glVertex4f(-1.0, 0.0, 0.0, 0.0);
  glEnd();
  glDisable(GL_LINE_STIPPLE);

  // Draw the positive y axis
  qglColor(QColor{0, 255, 0});
  glLineWidth(2.0);
  glBegin(GL_LINES);
  glVertex4f(0.0, 0.0, 0.0, 1.0);
  glVertex4f(0.0, 1.0, 0.0, 0.0);
  glEnd();

  // Draw the negative y axis
  qglColor(QColor{0, 255, 0});
  glLineWidth(2.0);
  glLineStipple(8, 0xAAAA);
  glEnable(GL_LINE_STIPPLE);
  glBegin(GL_LINES);
  glVertex4f(0.0, 0.0, 0.0, 1.0);
  glVertex4f(0.0, -1.0, 0.0, 0.0);
  glEnd();
  glDisable(GL_LINE_STIPPLE);
}

void GLWidget::renderAtFPS(const bool render_at_fps) {
  m_render_at_fps = render_at_fps;
}

void GLWidget::lockCamera(const bool lock_camera) {
  m_lock_camera = lock_camera;
}

void GLWidget::toggleHUD() {
  m_display_HUD = !m_display_HUD;

  updateGL();
}

void GLWidget::toggleGridDisplay() {
  m_display_grid = !m_display_grid;

  updateGL();
}

void GLWidget::toggleMPMDisplay() {
  m_display_mpm = !m_display_mpm;

  updateGL();
}

void GLWidget::toggleDiscreteDisplay() {
  m_display_discrete = !m_display_discrete;

  updateGL();
}

void GLWidget::toggleZoneIndicatorsDisplay() {
  m_display_zone_indicators = !m_display_zone_indicators;

  updateGL();
}

void GLWidget::centerCamera(const bool update_gl) {
  if (m_lock_camera) {
    return;
  }

  if (m_sim.discreteState().q().size() == 0) {
    m_camera_controller.reset();
  } else {
    const Array4s bbox{m_sim.discreteState().computeBoundingBox()};
    const scalar &minx{bbox(0)};
    const scalar &maxx{bbox(2)};
    assert(minx < maxx);
    const scalar &miny{bbox(1)};
    const scalar &maxy{bbox(3)};
    assert(miny < maxy);

    const scalar cx{minx + 0.5 * (maxx - minx)};
    const scalar rx{maxx - cx};
    const scalar cy{miny + 0.5 * (maxy - miny)};
    const scalar ry{maxy - cy};

    const scalar ratio{scalar(m_camera_controller.height()) /
                       scalar(m_camera_controller.width())};
    const scalar size{1.2 * std::max(ratio * rx, ry)};

    m_camera_controller.setCenter(cx, cy);
    m_camera_controller.setScaleFactor(size);
  }

  if (update_gl) {
    GLint width;
    GLint height;
    getViewportDimensions(width, height);
    m_camera_controller.reshape(width, height);
    updateGL();
  }
}

void GLWidget::saveScreenshot(const QString &file_name) {
  std::cout << "Saving screenshot of time "
            << scalar(
                   std::intmax_t(m_sim.integratorState().overallIteration()) *
                   m_sim.integratorState().overallTimestep())
            << " to " << file_name.toStdString() << std::endl;
  const QImage frame_buffer{grabFrameBuffer()};
  frame_buffer.save(file_name);
}

#ifdef USE_HDF5
void GLWidget::saveHDF5File(const QString &file_name) {
  std::cout << "Saving H5 state of time "
            << scalar(
                   std::intmax_t(m_sim.integratorState().overallIteration()) *
                   m_sim.integratorState().overallTimestep())
            << " to " << file_name.toStdString() << std::endl;

  // Save the simulation state
  try {
    HDF5File output_file{file_name.toStdString(), HDF5AccessType::READ_WRITE};
    // Save the iteration and time step and time
    output_file.writeScalar("", "global_timestep",
                            scalar(m_sim.integratorState().overallTimestep()));
    output_file.writeScalar("", "global_iteration",
                            m_sim.integratorState().overallIteration());
    output_file.writeScalar(
        "", "global_time",
        scalar(std::intmax_t(m_sim.integratorState().overallIteration()) *
               m_sim.integratorState().overallTimestep()));
    // Save out the git hash
    output_file.writeString("", "git_hash", CompileDefinitions::GitSHA1);
    // Tag this as 2d hybrid simulation
    output_file.writeString("", "sim_type", "hybrid2d");
    m_sim.writeBinaryState(output_file);
  } catch (const std::string &error) {
    std::cerr << error << std::endl;
    std::exit(EXIT_FAILURE);
  }
}
#endif

void GLWidget::setMovieDir(const QString &dir_name) {
  m_movie_dir_name = dir_name;
  m_output_frame = 0;

  // Save a screenshot of the current state
  if (m_movie_dir_name.size() != 0) {
    m_movie_dir.setPath(m_movie_dir_name);
    assert(m_movie_dir.exists());

    QString output_image_name{QString{tr("frame%1.png")}.arg(
        m_output_frame, 10, 10, QLatin1Char('0'))};
    saveScreenshot(m_movie_dir.filePath(output_image_name));
#ifdef USE_HDF5
    QString output_hdf5_name{QString{tr("frame%1.h5")}.arg(
        m_output_frame, 10, 10, QLatin1Char('0'))};
    saveHDF5File(m_movie_dir.filePath(output_hdf5_name));
#endif
    ++m_output_frame;
  }
}

void GLWidget::setMovieFPS(const unsigned fps) {
  assert(fps > 0);
  m_output_fps = fps;
  m_output_frame = 0;

  if (Rational<std::intmax_t>{1, std::intmax_t(m_output_fps)} <
      m_sim.integratorState().overallTimestep()) {
    std::cerr << "Warning, requested movie frame rate faster than timestep. "
                 "Dumping at timestep rate."
              << std::endl;
    m_steps_per_frame = 1;
  } else {
    const Rational<std::intmax_t> potential_steps_per_frame{
        std::intmax_t(1) / (m_sim.integratorState().overallTimestep() *
                            std::intmax_t(m_output_fps))};
    if (!potential_steps_per_frame.isInteger()) {
      if (m_sim.integratorState().overallTimestep() !=
          Rational<std::intmax_t>{0}) {
        std::cerr
            << "Warning, timestep and output frequency do not yield an integer "
               "number of timesteps for data output. Dumping at timestep rate."
            << std::endl;
      }
      m_steps_per_frame = 1;
    } else {
      m_steps_per_frame = unsigned(potential_steps_per_frame.numerator());
    }
  }
}

void GLWidget::exportCameraSettings() {
  std::cout << "<camera center=\"" << m_camera_controller.centerX() << " "
            << m_camera_controller.centerY() << "\" scale=\""
            << m_camera_controller.scaleFactor() << "\" fps=\"" << m_output_fps
            << "\" render_at_fps=\"" << m_render_at_fps << "\" locked=\""
            << m_lock_camera << "\"/>" << std::endl;
}

static void paintInfiniteLine(const Vector2s &x, const Vector2s &n) {
  const scalar theta{-scalar(180.0) * atan2(n.x(), n.y()) /
                     MathDefines::PI<scalar>()};

  glPushMatrix();

  glTranslated(GLdouble(x.x()), GLdouble(x.y()), GLdouble(0.0));
  glRotated(GLdouble(theta), GLdouble(0.0), GLdouble(0.0), GLdouble(1.0));

  glBegin(GL_LINES);
  glVertex4d(0.0, 0.0, 0.0, 1.0);
  glVertex4d(1.0, 0.0, 0.0, 0.0);
  glVertex4d(0.0, 0.0, 0.0, 1.0);
  glVertex4d(-1.0, 0.0, 0.0, 0.0);
  glEnd();

  glPopMatrix();
}

// TODO: Abstract the shared code in here into its own function
static void paintPlanarPortal(const PlanarPortal &planar_portal) {
  // Draw the first plane of the portal
  {
    const scalar theta{
        -180.0 *
        atan2(planar_portal.planeA().n().x(), planar_portal.planeA().n().y()) /
        MathDefines::PI<scalar>()};

    glPushMatrix();
    glTranslated(GLdouble(planar_portal.planeA().x().x()),
                 GLdouble(planar_portal.planeA().x().y()), 0.0);
    glRotated(GLdouble(theta), 0.0, 0.0, 1.0);

    glLineStipple(8, 0xAAAA);
    glEnable(GL_LINE_STIPPLE);
    glBegin(GL_LINES);
    glVertex4d(0.0, 0.0, 0.0, 1.0);
    glVertex4d(-1.0, 0.0, 0.0, 0.0);
    glEnd();
    glDisable(GL_LINE_STIPPLE);

    glLineStipple(8, 0x5555);
    glEnable(GL_LINE_STIPPLE);
    glBegin(GL_LINES);
    glVertex4d(0.0, 0.0, 0.0, 1.0);
    glVertex4d(1.0, 0.0, 0.0, 0.0);
    glEnd();
    glDisable(GL_LINE_STIPPLE);

    glPopMatrix();

    assert(planar_portal.bounds() >= 0.0);
    if (planar_portal.isLeesEdwards()) {
      // Draw the lower bound
      glPushMatrix();
      glTranslated(GLdouble(planar_portal.planeA().x().x()),
                   GLdouble(planar_portal.planeA().x().y()), 0.0);
      glRotated(GLdouble(theta), 0.0, 0.0, 1.0);
      glTranslated(-planar_portal.bounds(), 0.0, 0.0);
      glLineStipple(8, 0x5555);
      glEnable(GL_LINE_STIPPLE);
      glBegin(GL_LINES);
      glVertex4d(0.0, 0.0, 0.0, 1.0);
      glVertex4d(0.0, -1.0, 0.0, 0.0);
      glEnd();
      glDisable(GL_LINE_STIPPLE);
      glPopMatrix();

      // Draw the upper bound
      glPushMatrix();
      glTranslated(GLdouble(planar_portal.planeA().x().x()),
                   GLdouble(planar_portal.planeA().x().y()), 0.0);
      glRotated(GLdouble(theta), 0.0, 0.0, 1.0);
      glTranslated(planar_portal.bounds(), 0.0, 0.0);
      glLineStipple(8, 0x5555);
      glEnable(GL_LINE_STIPPLE);
      glBegin(GL_LINES);
      glVertex4d(0.0, 0.0, 0.0, 1.0);
      glVertex4d(0.0, -1.0, 0.0, 0.0);
      glEnd();
      glDisable(GL_LINE_STIPPLE);
      glPopMatrix();

      // Draw a short line to indicate the rest position of the center of the
      // portal
      glPushMatrix();
      glTranslated(planar_portal.planeA().x().x(),
                   GLdouble(planar_portal.planeA().x().y()), 0.0);
      glRotated(GLdouble(theta), 0.0, 0.0, 1.0);

      glLineStipple(8, 0x5555);
      glEnable(GL_LINE_STIPPLE);
      glBegin(GL_LINES);
      glVertex2d(0.0, 0.0);
      glVertex2d(0.0, -0.2 * planar_portal.bounds());
      glEnd();
      glDisable(GL_LINE_STIPPLE);

      glPopMatrix();

      // Draw a short line to indicate the current position of the center of the
      // portal
      glPushMatrix();
      glTranslated(GLdouble(planar_portal.transformedAx().x()),
                   GLdouble(planar_portal.transformedAx().y()), 0.0);
      glRotated(GLdouble(theta), 0.0, 0.0, 1.0);

      glLineStipple(8, 0x5555);
      glEnable(GL_LINE_STIPPLE);
      glBegin(GL_LINES);
      glVertex2d(0.0, 0.0);
      glVertex2d(0.0, -0.2 * planar_portal.bounds());
      glEnd();
      glDisable(GL_LINE_STIPPLE);

      glPopMatrix();
    } else {
      glPushMatrix();
      glTranslated(GLdouble(planar_portal.planeA().x().x()),
                   GLdouble(planar_portal.planeA().x().y()), 0.0);
      glRotated(GLdouble(theta), 0.0, 0.0, 1.0);

      // Draw an infinite line to show what half of portal is free
      glLineStipple(8, 0x5555);
      glEnable(GL_LINE_STIPPLE);
      glBegin(GL_LINES);
      glVertex4d(0.0, 0.0, 0.0, 1.0);
      glVertex4d(0.0, -1.0, 0.0, 0.0);
      glEnd();
      glDisable(GL_LINE_STIPPLE);

      glPopMatrix();
    }
  }

  // Draw the second plane of the portal
  {
    const scalar theta{
        -180.0 *
        atan2(planar_portal.planeB().n().x(), planar_portal.planeB().n().y()) /
        MathDefines::PI<scalar>()};

    glPushMatrix();
    glTranslated(GLdouble(planar_portal.planeB().x().x()),
                 GLdouble(planar_portal.planeB().x().y()), 0.0);
    glRotated(GLdouble(theta), 0.0, 0.0, 1.0);

    glLineStipple(8, 0xAAAA);
    glEnable(GL_LINE_STIPPLE);
    glBegin(GL_LINES);
    glVertex4d(0.0, 0.0, 0.0, 1.0);
    glVertex4d(-1.0, 0.0, 0.0, 0.0);
    glEnd();
    glDisable(GL_LINE_STIPPLE);

    glLineStipple(8, 0x5555);
    glEnable(GL_LINE_STIPPLE);
    glBegin(GL_LINES);
    glVertex4d(0.0, 0.0, 0.0, 1.0);
    glVertex4d(1.0, 0.0, 0.0, 0.0);
    glEnd();
    glDisable(GL_LINE_STIPPLE);

    glPopMatrix();

    assert(planar_portal.bounds() >= 0.0);
    if (planar_portal.isLeesEdwards()) {
      // Draw the lower bound
      glPushMatrix();
      glTranslated(GLdouble(planar_portal.planeB().x().x()),
                   GLdouble(planar_portal.planeB().x().y()), 0.0);
      glRotated(GLdouble(theta), 0.0, 0.0, 1.0);
      glTranslated(-planar_portal.bounds(), 0.0, 0.0);
      glLineStipple(8, 0x5555);
      glEnable(GL_LINE_STIPPLE);
      glBegin(GL_LINES);
      glVertex4d(0.0, 0.0, 0.0, 1.0);
      glVertex4d(0.0, -1.0, 0.0, 0.0);
      glEnd();
      glDisable(GL_LINE_STIPPLE);
      glPopMatrix();

      // Draw the upper bound
      glPushMatrix();
      glTranslated(GLdouble(planar_portal.planeB().x().x()),
                   GLdouble(planar_portal.planeB().x().y()), 0.0);
      glRotated(GLdouble(theta), 0.0, 0.0, 1.0);
      glTranslated(planar_portal.bounds(), 0.0, 0.0);
      glLineStipple(8, 0x5555);
      glEnable(GL_LINE_STIPPLE);
      glBegin(GL_LINES);
      glVertex4d(0.0, 0.0, 0.0, 1.0);
      glVertex4d(0.0, -1.0, 0.0, 0.0);
      glEnd();
      glDisable(GL_LINE_STIPPLE);
      glPopMatrix();

      // Draw a short line to indicate the rest position of the center of the
      // portal
      glPushMatrix();
      glTranslated(planar_portal.planeB().x().x(),
                   GLdouble(planar_portal.planeB().x().y()), 0.0);
      glRotated(GLdouble(theta), 0.0, 0.0, 1.0);

      glLineStipple(8, 0x5555);
      glEnable(GL_LINE_STIPPLE);
      glBegin(GL_LINES);
      glVertex2d(0.0, 0.0);
      glVertex2d(0.0, -0.2 * planar_portal.bounds());
      glEnd();
      glDisable(GL_LINE_STIPPLE);

      glPopMatrix();

      // Draw a short line to indicate the current position of the center of the
      // portal
      glPushMatrix();
      glTranslated(GLdouble(planar_portal.transformedBx().x()),
                   GLdouble(planar_portal.transformedBx().y()), 0.0);
      glRotated(GLdouble(theta), 0.0, 0.0, 1.0);

      glLineStipple(8, 0x5555);
      glEnable(GL_LINE_STIPPLE);
      glBegin(GL_LINES);
      glVertex2d(0.0, 0.0);
      glVertex2d(0.0, -0.2 * planar_portal.bounds());
      glEnd();
      glDisable(GL_LINE_STIPPLE);

      glPopMatrix();
    } else {
      glPushMatrix();
      glTranslated(GLdouble(planar_portal.planeB().x().x()),
                   GLdouble(planar_portal.planeB().x().y()), 0.0);
      glRotated(GLdouble(theta), 0.0, 0.0, 1.0);

      // Draw an infinite line to show what half of portal is free
      glLineStipple(8, 0x5555);
      glEnable(GL_LINE_STIPPLE);
      glBegin(GL_LINES);
      glVertex4d(0.0, 0.0, 0.0, 1.0);
      glVertex4d(0.0, -1.0, 0.0, 0.0);
      glEnd();
      glDisable(GL_LINE_STIPPLE);

      glPopMatrix();
    }
  }
}

static void drawMaterialPoints(const GLCircleRenderer2D &circle_renderer,
                               const MaterialPoints &points) {
  glPushAttrib(GL_COLOR);
  for (unsigned i = 0; i < points.npoints; ++i) {
    glPushMatrix();
    glTranslated(points.q(0, i), points.q(1, i), 0.0);
    glScaled(0.05, 0.05, 1.0);
    if (points.kinematically_scripted[i]) {
      circle_renderer.renderSolidCircle(
          Eigen::Matrix<GLdouble, 3, 1>{0.5, 0.5, 0.5});
    } else {
      circle_renderer.renderSolidCircle(
          Eigen::Matrix<GLdouble, 3, 1>{1.0, 0.54, 0.9});
    }
    glPopMatrix();
  }
  glPopAttrib();
}

static void drawGrid(const PhysicsGrid &grid) {
  if ((grid.cell_count.array() == 0).all()) {
    return;
  }

  const Array2s upper_left =
      grid.min + Array2s{0.0, grid.cell_width * grid.cell_count.y()};
  const Array2s lower_right =
      grid.min + Array2s{grid.cell_width * grid.cell_count.x(), 0.0};

  glPushAttrib(GL_LINE_WIDTH);
  glPushAttrib(GL_COLOR);

  glLineWidth(0.5);
  glColor3d(0, 0, 0);

  glBegin(GL_LINES);
  for (unsigned vertical_line_num = 0;
       vertical_line_num < grid.cell_count.x() + 1; vertical_line_num++) {
    const Array2s x0 =
        grid.min + Array2s{vertical_line_num * grid.cell_width, 0.0};
    const Array2s x1 =
        upper_left + Array2s{vertical_line_num * grid.cell_width, 0.0};
    glVertex2d(x0.x(), x0.y());
    glVertex2d(x1.x(), x1.y());
  }
  for (unsigned horizontal_line_num = 0;
       horizontal_line_num < grid.cell_count.y() + 1; horizontal_line_num++) {
    const Array2s x0 =
        grid.min + Array2s{0.0, horizontal_line_num * grid.cell_width};
    const Vector2s x1 =
        lower_right + Array2s{0.0, horizontal_line_num * grid.cell_width};
    glVertex2d(x0.x(), x0.y());
    glVertex2d(x1.x(), x1.y());
  }
  glEnd();

  glPopAttrib();
  glPopAttrib();
}

static void renderBox(const Eigen::Matrix<GLdouble, 3, 1> &color,
                      const Vector2s &r) {
  glPushMatrix();
  glScaled(GLdouble(r.x()), GLdouble(r.y()), GLdouble(1.0));
  glPushAttrib(GL_COLOR);
  glColor3d(color.x(), color.y(), color.z());
  glBegin(GL_TRIANGLE_STRIP);
  glVertex2d(-1.0, 1.0);
  glVertex2d(-1.0, -1.0);
  glVertex2d(1.0, -1.0);
  glVertex2d(1.0, 1.0);
  glVertex2d(-1.0, 1.0);
  glEnd();
  glPopAttrib();
  glPopMatrix();
}

static void renderTeleportedBox(const Eigen::Matrix<GLdouble, 3, 1> &color,
                                const Vector2s &r) {
  glPushMatrix();
  glScaled(GLdouble(r.x()), GLdouble(r.y()), GLdouble(1.0));
  glPushAttrib(GL_COLOR);
  glColor3d(color.x(), color.y(), color.z());
  glPushAttrib(GL_LINE_WIDTH);
  glLineWidth(3.0);
  glBegin(GL_LINE_LOOP);
  glVertex2d(1.0, 1.0);
  glVertex2d(-1.0, 1.0);
  glVertex2d(-1.0, -1.0);
  glVertex2d(1.0, -1.0);
  glEnd();
  glPopAttrib();
  glPopAttrib();
  glPopMatrix();
}

static uint32_t hash(uint32_t a) {
  a = (a + 0x7ed55d16) + (a << 12);
  a = (a ^ 0xc761c23c) ^ (a >> 19);
  a = (a + 0x165667b1) + (a << 5);
  a = (a + 0xd3a2646c) ^ (a << 9);
  a = (a + 0xfd7046c5) + (a << 3);
  a = (a ^ 0xb55a4f09) ^ (a >> 16);
  return a;
}

static int hashInt(const int a) { return int(hash(uint32_t(a))); }

void GLWidget::paintSystem() const {
  if (m_display_discrete) {

    // Draw each static drum
    glPushAttrib(GL_COLOR);
    glPushAttrib(GL_LINE_WIDTH);
    glLineWidth(2.0);
    qglColor(QColor{0, 0, 0});
    {
      const std::vector<RigidBody2DStaticDrum> &drums =
          m_sim.discreteState().drums();
      for (const RigidBody2DStaticDrum &drum : drums) {
        glPushMatrix();
        glTranslated(GLdouble(drum.x().x()), GLdouble(drum.x().y()),
                     GLdouble(0.0));
        glScaled(GLdouble(drum.r()), GLdouble(drum.r()), GLdouble(1.0));
        const scalar theta_degrees{180.0 * drum.theta() /
                                   MathDefines::PI<scalar>()};
        glRotated(theta_degrees, 0.0, 0.0, 1.0);
        m_circle_renderer.renderCutCircle(
            Eigen::Matrix<GLdouble, 3, 1>(0.3, 0.3, 0.3));
        glPopMatrix();
      }
    }
    glPopAttrib();
    glPopAttrib();

    // Draw each body
    {
      const VectorXs &q{m_sim.discreteState().q()};
      assert(q.size() % 3 == 0);
      const unsigned nbodies{static_cast<unsigned>(q.size() / 3)};

      for (unsigned bdy_idx = 0; bdy_idx < nbodies; bdy_idx++)
      // for( const unsigned bdy_idx : m_sim.bodiesToStabilize() )
      {
        glPushMatrix();
        glTranslated(GLdouble(q(3 * bdy_idx)), GLdouble(q(3 * bdy_idx + 1)),
                     GLdouble(0.0));
        const scalar theta_degrees{180.0 * q(3 * bdy_idx + 2) /
                                   MathDefines::PI<scalar>()};
        glRotated(theta_degrees, 0.0, 0.0, 1.0);
        assert(bdy_idx < m_sim.discreteState().geometryIndices().size());

        const Vector3s color =
            m_sim.discreteSim().isKinematicallyScripted(bdy_idx)
                ? Vector3s{0.5, 0.5, 0.5}
                : m_template_colors[hashInt(m_sim.discreteState()
                                                .getUniqueBodyIndex(bdy_idx)) %
                                    m_template_colors.size()];

        const std::unique_ptr<RigidBody2DGeometry> &geo =
            m_sim.discreteState().bodyGeometry(bdy_idx);

        switch (geo->type()) {
        case RigidBody2DGeometryType::ANNULUS: {
          std::cerr << "Error, annulus rendering not coded up here yet!"
                    << std::endl;
          std::exit(EXIT_FAILURE);
          // break;
        }
        case RigidBody2DGeometryType::CIRCLE: {
          const CircleGeometry &circle{
              static_cast<const CircleGeometry &>(*geo.get())};
          glPushMatrix();
          glScaled(GLdouble(circle.r()), GLdouble(circle.r()), GLdouble(1.0));
          m_circle_renderer.renderCircle(color);
          glPopMatrix();
          break;
        }
        case RigidBody2DGeometryType::BOX: {
          const BoxGeometry &box{static_cast<const BoxGeometry &>(*geo.get())};
          renderBox(color, box.r());
          break;
        }
        }

        glPopMatrix();
      }
    }

    // Draw teleported versions of each body
    {
      const VectorXs &q{m_sim.discreteState().q()};
      assert(q.size() % 3 == 0);
      const unsigned nbodies{static_cast<unsigned>(q.size() / 3)};
      const std::vector<PlanarPortal> &planar_portals{
          m_sim.discreteState().planarPortals()};
      // For each planar portal
      // TODO: More efficient to invert these loops, only cache out pos, theta
      // once per body
      for (const PlanarPortal &planar_portal : planar_portals) {
        for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
          const Vector2s pos{q.segment<2>(3 * bdy_idx)};
          const scalar theta{q(3 * bdy_idx + 2)};
          // Compute the AABB for the current body
          Array2s min;
          Array2s max;
          m_sim.discreteState().bodyGeometry(bdy_idx)->computeAABB(pos, theta,
                                                                   min, max);
          assert((min < max).all());
          // If the AABB intersects a periodic boundary
          bool intersecting_index;
          if (planar_portal.aabbTouchesPortal(min, max, intersecting_index)) {
            Vector2s teleported_pos;
            planar_portal.teleportPoint(pos, intersecting_index,
                                        teleported_pos);
            glPushMatrix();
            glTranslated(GLdouble(teleported_pos.x()),
                         GLdouble(teleported_pos.y()), GLdouble(0.0));
            const scalar theta_degrees{180.0 * q(3 * bdy_idx + 2) /
                                       MathDefines::PI<scalar>()};
            glRotated(theta_degrees, 0.0, 0.0, 1.0);
            assert(bdy_idx < m_sim.discreteState().geometryIndices().size());

            const Vector3s color =
                m_sim.discreteSim().isKinematicallyScripted(bdy_idx)
                    ? Vector3s{0.5, 0.5, 0.5}
                    : m_template_colors[hashInt(
                                            m_sim.discreteState()
                                                .getUniqueBodyIndex(bdy_idx)) %
                                        m_template_colors.size()];

            const std::unique_ptr<RigidBody2DGeometry> &geo =
                m_sim.discreteState().bodyGeometry(bdy_idx);

            switch (geo->type()) {
            case RigidBody2DGeometryType::ANNULUS: {
              std::cerr << "Error, annulus not coded up here!" << std::endl;
              std::exit(EXIT_FAILURE);
              // break;
            }
            case RigidBody2DGeometryType::CIRCLE: {
              const CircleGeometry &circle{
                  static_cast<const CircleGeometry &>(*geo.get())};
              glPushMatrix();
              glScaled(GLdouble(circle.r()), GLdouble(circle.r()),
                       GLdouble(1.0));
              m_circle_renderer.renderCircleOutline(color);
              glPopMatrix();
              break;
            }
            case RigidBody2DGeometryType::BOX: {
              const BoxGeometry &box{
                  static_cast<const BoxGeometry &>(*geo.get())};
              renderTeleportedBox(color, box.r());
              break;
            }
            }

            glPopMatrix();
          }
        }
      }
    }

    // Draw each planar portal
    glPushAttrib(GL_COLOR);
    glPushAttrib(GL_LINE_WIDTH);
    glLineWidth(2.0);
    {
      // TODO: Create a set number of nice looking colors for the portal ahead
      // of time instead of regenerating them
      std::mt19937_64 mt{123456};
      std::uniform_int_distribution<int> color_gen{0, 255};
      const std::vector<PlanarPortal> &planar_portals{
          m_sim.discreteState().planarPortals()};
      for (const PlanarPortal &planar_portal : planar_portals) {
        const int r{color_gen(mt)};
        const int g{color_gen(mt)};
        const int b{color_gen(mt)};
        qglColor(QColor{r, g, b});
        paintPlanarPortal(planar_portal);
      }
    }
    glPopAttrib();
    glPopAttrib();

    // Draw each static plane
    glPushAttrib(GL_COLOR);
    glPushAttrib(GL_LINE_WIDTH);
    glLineWidth(2.0);
    qglColor(QColor{0, 0, 0});
    {
      const std::vector<RigidBody2DStaticPlane> &planes{
          m_sim.discreteState().planes()};
      for (const RigidBody2DStaticPlane &plane : planes) {
        paintInfiniteLine(plane.x(), plane.n());
      }
    }
    glPopAttrib();
    glPopAttrib();
  }

  if (m_display_mpm) {
    drawMaterialPoints(m_circle_renderer,
                       m_sim.continuumState().material_points);
  }

  if (m_display_grid) {
    drawGrid(m_sim.continuumState().physics_grid);
  }
}

static QString generateTimeString(const unsigned iteration,
                                  const Rational<std::intmax_t> &dt,
                                  const int display_precision,
                                  const scalar &end_time) {
  QString time_string{QObject::tr("  t: ")};
  time_string +=
      QString::number(iteration * scalar(dt), 'f', display_precision);
  if (end_time != SCALAR_INFINITY) {
    time_string += QString{QObject::tr(" / ")};
    time_string += QString::number(end_time);
  }
  return time_string;
}

void GLWidget::paintHUD() {
  static int text_width{0};

  glPushAttrib(GL_MATRIX_MODE);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();

  // Set an orthographic projection with height and width equal to window height
  // and width
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  GLint width;
  GLint height;
  getViewportDimensions(width, height);
  glOrtho(0, width, 0, height, -1, 1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Enable blending for transparent HUD elements
  glPushAttrib(GL_BLEND);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // Draw a semi-transparent overlay so text is visible regardless of background
  // color
  const Eigen::Matrix<GLdouble, 2, 1> overlay_start{0, height - 1 * 12 - 2};
  const Eigen::Matrix<GLdouble, 2, 1> overlay_extnt{text_width + 2 + 2, height};
  glColor4d(0.0, 0.0, 0.0, 0.5);
  glBegin(GL_QUADS);
  glVertex2d(GLdouble(overlay_start.x()), GLdouble(overlay_start.y()));
  glVertex2d(GLdouble(overlay_start.x() + overlay_extnt.x()),
             GLdouble(overlay_start.y()));
  glVertex2d(GLdouble(overlay_start.x() + overlay_extnt.x()),
             GLdouble(overlay_start.y() + overlay_extnt.y()));
  glVertex2d(GLdouble(overlay_start.x()),
             GLdouble(overlay_start.y() + overlay_extnt.y()));
  glEnd();

  glDisable(GL_BLEND);
  glPopAttrib();

  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glPopAttrib();

  // String to display in upper left corner
  const QString time_string{generateTimeString(
      m_sim.integratorState().overallIteration(),
      m_sim.integratorState().overallTimestep(), m_display_precision,
      scalar(m_sim.integratorState().endTime()))};
  //  const QString delta_H{ generateNumericString( " dH: ", m_delta_H0 ) };
  //  const QString delta_px{ generateNumericString( "dpx: ", m_delta_p0.x() )
  //  }; const QString delta_py{ generateNumericString( "dpy: ", m_delta_p0.y()
  //  ) }; const QString delta_L{ generateNumericString( " dL: ", m_delta_L0 )
  //  };
  {
    const QFontMetrics font_metrics{QFont{"Courier", 12}};
    text_width =
        std::max(text_width, font_metrics.boundingRect(time_string).width());
    //    text_width = std::max( text_width, font_metrics.boundingRect( delta_H
    //    ).width() ); text_width = std::max( text_width,
    //    font_metrics.boundingRect( delta_px ).width() ); text_width =
    //    std::max( text_width, font_metrics.boundingRect( delta_py ).width() );
    //    text_width = std::max( text_width, font_metrics.boundingRect( delta_L
    //    ).width() );
  }

  qglColor(QColor{255, 255, 255});
  const QFont font{"Courier", 12};
  renderText(2, font.pointSize(), time_string, font);
  //  renderText( 2, 2 * font.pointSize(), delta_H, font );
  //  renderText( 2, 3 * font.pointSize(), delta_px, font );
  //  renderText( 2, 4 * font.pointSize(), delta_py, font );
  //  renderText( 2, 5 * font.pointSize(), delta_L, font );

  assert(checkGLErrors());
}

void GLWidget::mousePressEvent(QMouseEvent *event) {
  if (m_lock_camera) {
    return;
  }

  bool repaint_needed{false};

  if (event->buttons() & Qt::LeftButton) {
    m_left_mouse_button_pressed = true;
    repaint_needed = true;
  }
  if (event->buttons() & Qt::RightButton) {
    m_right_mouse_button_pressed = true;
  }

  if (repaint_needed) {
    updateGL();
  }

  m_last_pos = event->pos();
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event) {
  if (m_lock_camera) {
    return;
  }

  bool repaint_needed{false};

  if (!(event->buttons() & Qt::LeftButton) && m_left_mouse_button_pressed) {
    m_left_mouse_button_pressed = false;
    repaint_needed = true;
  }
  if (!(event->buttons() & Qt::RightButton) && m_right_mouse_button_pressed) {
    m_right_mouse_button_pressed = false;
  }

  if (repaint_needed) {
    updateGL();
  }
}

void GLWidget::mouseMoveEvent(QMouseEvent *event) {
  if (m_lock_camera) {
    return;
  }

  const int dx{event->x() - m_last_pos.x()};
  const int dy{event->y() - m_last_pos.y()};
  m_last_pos = event->pos();

  bool repaint_needed{false};

  if (event->buttons() & Qt::LeftButton) {
    assert(m_left_mouse_button_pressed);
    m_camera_controller.translateView(dx, dy);
    repaint_needed = true;
  }

  if (event->buttons() & Qt::RightButton) {
    assert(m_right_mouse_button_pressed);
    m_camera_controller.zoomView(0.02 * m_camera_controller.scaleFactor() * dy);
    repaint_needed = true;
  }

  if (repaint_needed) {
    updateGL();
  }

  assert(checkGLErrors());
}

void GLWidget::wheelEvent(QWheelEvent *event) {
  if (m_lock_camera) {
    return;
  }

  m_camera_controller.zoomView(-0.002 * m_camera_controller.scaleFactor() *
                               event->delta());
  updateGL();
  assert(checkGLErrors());
}
