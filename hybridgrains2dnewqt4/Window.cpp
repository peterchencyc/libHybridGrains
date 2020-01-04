#include "Window.h"

#include <cassert>

#include <QKeyEvent>
#include <QMenu>
#include <QMenuBar>

#include "ContentWidget.h"

Window::Window(const QString &scene_name, QWidget *parent)
    : QMainWindow(parent), m_content_widget(nullptr) {
  m_content_widget = new ContentWidget{scene_name, this};

  QMenu *file{menuBar()->addMenu(tr("File"))};

  QMenu *view{menuBar()->addMenu(tr("View"))};

  QAction *separator{new QAction(this)};
  separator->setSeparator(true);

  // File menu actions

  // Load the input xml file
  QAction *open_scene{new QAction{tr("Open..."), this}};
  open_scene->setShortcut(tr("Ctrl+o"));
  file->addAction(open_scene);
  connect(open_scene, SIGNAL(triggered()), m_content_widget, SLOT(openScene()));

  // Reload the current xml file
  QAction *reload_scene{new QAction{tr("Reload"), this}};
  reload_scene->setShortcut(tr("Ctrl+r"));
  file->addAction(reload_scene);
  connect(reload_scene, SIGNAL(triggered()), m_content_widget,
          SLOT(reloadScene()));

  // Add a separator
  file->addAction(separator);

  // Export an image of the scene
  QAction *export_image{new QAction{tr("Export Image..."), this}};
  export_image->setShortcut(tr("Ctrl+i"));
  file->addAction(export_image);
  connect(export_image, SIGNAL(triggered()), m_content_widget,
          SLOT(exportImage()));

  // Export a movie of the scene
  QAction *export_movie{new QAction{tr("Export Movie..."), this}};
  export_movie->setShortcut(tr("Ctrl+m"));
  file->addAction(export_movie);
  connect(export_movie, SIGNAL(triggered()), m_content_widget,
          SLOT(exportMovie()));

  // Add a separator
  QAction *separator2{new QAction{this}};
  separator2->setSeparator(true);
  file->addAction(separator2);

  // Export the current camera settings
  QAction *export_camera_settings{new QAction{tr("Export Camera..."), this}};
  file->addAction(export_camera_settings);
  connect(export_camera_settings, SIGNAL(triggered()), m_content_widget,
          SLOT(exportCameraSettings()));

  // View menu actions

  // Toggle the heads up display
  QAction *toggle_hud{new QAction(tr("Togge HUD"), this)};
  toggle_hud->setShortcut(tr("h"));
  view->addAction(toggle_hud);
  connect(toggle_hud, SIGNAL(triggered()), m_content_widget, SLOT(toggleHUD()));

  // Toggle the grid rendering
  QAction *toggle_grid_display{new QAction(tr("Togge Grid Display"), this)};
  // toggle_hud->setShortcut( tr( "h" ) );
  view->addAction(toggle_grid_display);
  connect(toggle_grid_display, SIGNAL(triggered()), m_content_widget,
          SLOT(toggleGridDisplay()));

  // Toggle the mpm sim rendering
  QAction *toggle_mpm_display{new QAction(tr("Toggle MPM Display"), this)};
  // toggle_hud->setShortcut( tr( "h" ) );
  view->addAction(toggle_mpm_display);
  connect(toggle_mpm_display, SIGNAL(triggered()), m_content_widget,
          SLOT(toggleMPMDisplay()));

  // Toggle the discrete sim rendering
  QAction *toggle_discrete_display{
      new QAction(tr("Toggle Discrete Display"), this)};
  //  //toggle_hud->setShortcut( tr( "h" ) );
  view->addAction(toggle_discrete_display);
  connect(toggle_discrete_display, SIGNAL(triggered()), m_content_widget,
          SLOT(toggleDiscreteDisplay()));

  // Toggle the zone indicator rendering
  QAction *toggle_zone_indicators_display{
      new QAction(tr("Toggle Zone Indicators Display"), this)};
  //  //toggle_hud->setShortcut( tr( "h" ) );
  view->addAction(toggle_zone_indicators_display);
  connect(toggle_zone_indicators_display, SIGNAL(triggered()), m_content_widget,
          SLOT(toggleZoneIndicatorsDisplay()));

  // Add a separator
  view->addAction(separator);

  // Center the camera
  QAction *center_camera{new QAction(tr("Center Camera"), this)};
  center_camera->setShortcut(tr("c"));
  view->addAction(center_camera);
  connect(center_camera, SIGNAL(triggered()), m_content_widget,
          SLOT(centerCamera()));

  setCentralWidget(m_content_widget);
}

void Window::keyPressEvent(QKeyEvent *event) {
  assert(event != nullptr);

  if (event->key() == Qt::Key_Space) {
    m_content_widget->toggleSimulationCheckbox();
  } else if (event->key() == Qt::Key_R) {
    m_content_widget->resetSystem();
  } else if (event->key() == Qt::Key_S) {
    m_content_widget->takeStep();
  }
  // else if( event->key() == Qt::Key_T )
  // {
  //   m_content_widget->executeTest();
  // }
}

void Window::closeEvent(QCloseEvent *event) {
  assert(m_content_widget != nullptr);
  m_content_widget->close();
}
