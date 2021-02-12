#include "AnnulusGeometryRenderer.h"

#include "rigidbody2d/AnnulusGeometry.h"
#include "GLCircleRenderer2D.h"

#include <iostream>

AnnulusGeometryRenderer::AnnulusGeometryRenderer( const AnnulusGeometry& geo, const GLCircleRenderer2D& circle_renderer )
: m_annulus_geo( geo )
, m_circle_renderer( circle_renderer )
{}

void AnnulusGeometryRenderer::render( const Eigen::Matrix<GLdouble,3,1>& color )
{
  glPushMatrix();
  glScaled( GLdouble( m_annulus_geo.r0() ), GLdouble( m_annulus_geo.r0() ), GLdouble( 1.0 ) );
  m_circle_renderer.renderCircleOutline( color );
  glPopMatrix();

  glPushMatrix();
  glScaled( GLdouble( m_annulus_geo.r1() ), GLdouble( m_annulus_geo.r1() ), GLdouble( 1.0 ) );
  m_circle_renderer.renderCircleOutline( color );
  glPopMatrix();
}

void AnnulusGeometryRenderer::renderTeleported( const Eigen::Matrix<GLdouble,3,1>& color )
{
  std::cout << "AnnulusGeometryRenderer::renderTeleported" << std::endl;
//   glPushMatrix();
//   glScaled( GLdouble( m_circle_geo.r() ), GLdouble( m_circle_geo.r() ), GLdouble( 1.0 ) );
//   m_circle_renderer.renderCircleOutline( color );
//   glPopMatrix();
}
