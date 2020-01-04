#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QDir>
#include <QGLWidget>

#include "hybridgrains2dnew/HybridGrains2DSim.h"

#include "DisplayController2D.h"
#include "GLCircleRenderer2D.h"

class GLWidget : public QGLWidget {

  Q_OBJECT

public:
  explicit GLWidget(QWidget *parent = nullptr);
  ~GLWidget();

  QSize minimumSizeHint() const;
  QSize sizeHint() const;

  bool openScene(const QString &xml_scene_file_name, const bool &render_on_load,
                 unsigned &fps, bool &render_at_fps, bool &lock_camera);

  // Methods to control the solver
  void stepSystem();
  void resetSystem();

  void renderAtFPS(const bool render_at_fps);

  void lockCamera(const bool lock_camera);

  void toggleHUD();

  void toggleGridDisplay();
  void toggleMPMDisplay();
  void toggleDiscreteDisplay();
  void toggleZoneIndicatorsDisplay();

  void centerCamera(const bool update_gl = true);

  void saveScreenshot(const QString &file_name);
#ifdef USE_HDF5
  void saveHDF5File(const QString &file_name);
#endif

  void setMovieDir(const QString &dir_name);

  void setMovieFPS(const unsigned fps);

  void exportCameraSettings();

  // void executeTest();

protected:
  void initializeGL();
  void resizeGL(int width, int height);
  void paintGL();

  void mousePressEvent(QMouseEvent *event);
  void mouseReleaseEvent(QMouseEvent *event);
  void mouseMoveEvent(QMouseEvent *event);
  void wheelEvent(QWheelEvent *event);

private:
  bool axesDrawingIsEnabled() const;
  void paintAxes() const;

  void paintSystem() const;

  void paintHUD();

  DisplayController2D m_camera_controller;
  bool m_render_at_fps;
  bool m_lock_camera;
  QPoint m_last_pos;
  bool m_left_mouse_button_pressed;
  bool m_right_mouse_button_pressed;

  GLCircleRenderer2D m_circle_renderer;

  // Colors to use for bodies in the scene
  std::vector<Vector3s> m_template_colors;

  // Number of decimal places to display in time display
  int m_display_precision;

  bool m_display_HUD;

  bool m_display_grid;
  bool m_display_mpm;
  bool m_display_discrete;
  bool m_display_zone_indicators;

  // Directory to save periodic screenshots of the simulation into
  QString m_movie_dir_name;
  QDir m_movie_dir;
  // Number of frames that have been saved in the movie directory
  unsigned m_output_frame;
  // Rate at which to output movie frames
  unsigned m_output_fps;
  // Number of timesteps between frame outputs
  unsigned m_steps_per_frame;

  // Initial state of the simulation
  HybridGrains2DSim m_sim0;
  // Current state of the simulation
  HybridGrains2DSim m_sim;
};

#endif
