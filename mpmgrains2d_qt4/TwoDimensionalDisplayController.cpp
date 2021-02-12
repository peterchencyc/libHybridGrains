#include "TwoDimensionalDisplayController.h"

#include <algorithm>

TwoDimensionalDisplayController::TwoDimensionalDisplayController()
: m_window_width( 512 )
, m_window_height( 512 )
, m_scale_factor( 1.0 )
, m_center_x( 0.0 )
, m_center_y( 0.0 )
{}

void TwoDimensionalDisplayController::reshape( const unsigned w, const unsigned h )
{
  // Record the new width and height
  m_window_width = w;
  m_window_height = h;
  // Reset the coordinate system before modifying
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  // Set the coordinate system to achieve the desired zoom level, center
  const GLdouble ratio = GLdouble( h ) / GLdouble( w );
  const GLdouble left   = m_center_x - m_scale_factor / ratio;
  const GLdouble right  = m_center_x + m_scale_factor / ratio;
  const GLdouble bottom = m_center_y - m_scale_factor;
  const GLdouble top    = m_center_y + m_scale_factor;
  const GLdouble nearVal = -1.0;
  const GLdouble farVal = 1.0;
  glOrtho( left, right, bottom, top, nearVal, farVal );
  // Set the viewport to be the entire window
  glViewport( 0, 0, GLsizei(w), GLsizei(h) );
}

void TwoDimensionalDisplayController::translateView( const GLdouble& dx, const GLdouble& dy )
{
  const GLdouble translate_x = 2.0 * m_scale_factor * dx / GLdouble( m_window_height );
  const GLdouble translate_y = 2.0 * m_scale_factor * dy / GLdouble( m_window_height );
  m_center_x -= translate_x;
  m_center_y += translate_y;
  reshape( m_window_width, m_window_height );
}

void TwoDimensionalDisplayController::zoomView( const GLdouble& delta )
{
  m_scale_factor = std::max( 0.00001, m_scale_factor + delta );
  reshape( m_window_width, m_window_height );
}

void TwoDimensionalDisplayController::setCenter( const GLdouble& x, const GLdouble& y )
{
  m_center_x = x;
  m_center_y = y;
}

void TwoDimensionalDisplayController::setScaleFactor( const GLdouble& scale )
{
  m_scale_factor = scale;
}

unsigned TwoDimensionalDisplayController::width() const
{
  return m_window_width;
}

unsigned TwoDimensionalDisplayController::height() const
{
  return m_window_height;
}

const GLdouble& TwoDimensionalDisplayController::scaleFactor() const
{
  return m_scale_factor;
}

const GLdouble& TwoDimensionalDisplayController::centerX() const
{
  return m_center_x;
}

const GLdouble& TwoDimensionalDisplayController::centerY() const
{
  return m_center_y;
}

void TwoDimensionalDisplayController::reset()
{
  m_scale_factor = 1.0;
  m_center_x = 0.0;
  m_center_y = 0.0;
}
