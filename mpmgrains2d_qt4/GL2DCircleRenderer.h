#ifndef GL_2D_CIRCLE_RENDERER_H
#define GL_2D_CIRCLE_RENDERER_H

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <Eigen/Core>

// TODO: Write a version that uses display lists
class GL2DCircleRenderer final
{

public:

  // num_points_single_half: number of points in a single half of the circle
  explicit GL2DCircleRenderer( const unsigned num_pts = 16 );

  void renderCircle() const;

private:

  Eigen::Matrix<GLdouble,Eigen::Dynamic,1> m_x_crds;
  Eigen::Matrix<GLdouble,Eigen::Dynamic,1> m_y_crds;

};

#endif
