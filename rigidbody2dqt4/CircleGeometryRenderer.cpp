#include "CircleGeometryRenderer.h"

#include "rigidbody2d/CircleGeometry.h"
#include "GLCircleRenderer2D.h"

CircleGeometryRenderer::CircleGeometryRenderer( const CircleGeometry& geo, const GLCircleRenderer2D& circle_renderer )
: m_circle_geo( geo )
, m_circle_renderer( circle_renderer )
{}

void CircleGeometryRenderer::render( const Eigen::Matrix<GLdouble,3,1>& color )
{
  glPushMatrix();
  glScaled( GLdouble( m_circle_geo.r() ), GLdouble( m_circle_geo.r() ), GLdouble( 1.0 ) );
  m_circle_renderer.renderCircle( color );
  glPopMatrix();
}

void CircleGeometryRenderer::renderTeleported( const Eigen::Matrix<GLdouble,3,1>& color )
{
  glPushMatrix();
  glScaled( GLdouble( m_circle_geo.r() ), GLdouble( m_circle_geo.r() ), GLdouble( 1.0 ) );
  m_circle_renderer.renderCircleOutline( color );
  glPopMatrix();
}
