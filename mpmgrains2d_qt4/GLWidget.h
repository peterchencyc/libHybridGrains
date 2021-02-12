#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>
#include <QDir>

#include "mpmgrains2d/SimulationState.h"
#include "scisim/Math/Rational.h"

#include "TwoDimensionalDisplayController.h"
#include "GL2DCircleRenderer.h"

#include <random>

class GLWidget : public QGLWidget
{

  Q_OBJECT

public:

  explicit GLWidget( QWidget* parent = nullptr );
  ~GLWidget();

  QSize minimumSizeHint() const;
  QSize sizeHint() const;

  bool openScene( const QString& xml_scene_file_name, unsigned& fps, bool& render_at_fps, bool& lock_camera );

  // Methods to control the solver
  void stepSystem();
  void resetSystem();

  void renderAtFPS( const bool render_at_fps );

  void lockCamera( const bool lock_camera );

  void toggleHUD();
  void centerCamera( const bool update_gl );

  void saveScreenshot( const QString& file_name );
  #ifdef USE_HDF5
  void saveHDF5File( const QString& file_name );
  #endif

  void setMovieDir( const QString& dir_name );

  void setMovieFPS( const unsigned fps );

  void exportCameraSettings();

protected:

  void initializeGL();
  void resizeGL( int width, int height );
  void paintGL();

  void mousePressEvent( QMouseEvent* event );
  void mouseReleaseEvent( QMouseEvent* event );
  void mouseMoveEvent( QMouseEvent* event );
  void wheelEvent( QWheelEvent* event );

private:

  bool axesDrawingIsEnabled() const;
  void paintAxes() const;

  void drawPoints() const;

  void paintSystem() const;

  void paintHUD();

  TwoDimensionalDisplayController m_camera_controller;
  bool m_render_at_fps;
  bool m_lock_camera;
  QPoint m_last_pos;
  bool m_left_mouse_button_pressed;
  bool m_right_mouse_button_pressed;

  GL2DCircleRenderer m_circle_renderer;

  // Number of decimal places to display in time display
  const int m_display_precision;

  bool m_display_HUD;

  // Directory to save periodic screenshots of the simulation into
  QString m_movie_dir_name;
  QDir m_movie_dir;
  // Number of frames that have been saved in the movie directory
  unsigned m_output_frame;
  // Rate at which to output movie frames
  unsigned m_output_fps;
  // Number of timesteps between frame outputs
  unsigned m_steps_per_frame;

  // Current iteration of the solver
  unsigned m_iteration;
  // Current timestep
  Rational<std::intmax_t> m_dt;
  // End time of the simulation
  scalar m_end_time;
  
  SimulationState m_state;
};

#endif
