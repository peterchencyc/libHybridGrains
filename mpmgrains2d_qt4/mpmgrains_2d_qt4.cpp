#include <QApplication>
#include <QDesktopWidget>

#include "Window.h"
#include <iostream>

static void centerWindow( Window& window )
{
  const int screenWidth = QApplication::desktop()->screenGeometry().width();
  const int screenHeight = QApplication::desktop()->screenGeometry().height();

  const QSize windowSize = window.size();
  const int width = windowSize.width();
  const int height = windowSize.height();

  const int x = ( screenWidth - width ) / 2;
  const int y = ( screenHeight - height ) / 2;

  window.move( x, y );
}

int main( int argc, char** argv )
{
  QApplication app{ argc, argv };
  const QStringList arguments{ app.arguments() };
  if( arguments.count() > 2 )
  {
    std::cerr << "Error, must provide a valid configuration file name or no argument. Exiting." << std::endl;
    return EXIT_FAILURE;
  }
  Window window{ arguments.count() == 2 ? arguments[1] : "" };
  window.resize( window.sizeHint() );
  window.setWindowTitle( "MPM Grains 2D" );
  centerWindow( window );
  window.show();
  window.raise();
  return app.exec();
}
