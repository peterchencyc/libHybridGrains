#ifndef ANNULUS_GEOMETRY_RENDERER_H
#define ANNULUS_GEOMETRY_RENDERER_H

#include "RigidBodyRenderer2D.h"

class AnnulusGeometry;
class GLCircleRenderer2D;

class AnnulusGeometryRenderer final : public RigidBodyRenderer2D
{

public:

  AnnulusGeometryRenderer( const AnnulusGeometry& geo, const GLCircleRenderer2D& circle_renderer );

  virtual ~AnnulusGeometryRenderer() override = default;

  virtual void render( const Eigen::Matrix<GLdouble,3,1>& color ) override;

  virtual void renderTeleported( const Eigen::Matrix<GLdouble,3,1>& color ) override;

private:

  const AnnulusGeometry& m_annulus_geo;
  const GLCircleRenderer2D& m_circle_renderer;

};

#endif
