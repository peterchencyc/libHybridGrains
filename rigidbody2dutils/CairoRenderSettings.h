#ifndef CAIRO_RENDER_SETTINGS_H
#define CAIRO_RENDER_SETTINGS_H

#include "scisim/Math/MathDefines.h"

class CairoRenderSettings final {

public:
  CairoRenderSettings();

  const Vector2i &imageDimensions() const;
  const double &cameraScale() const;
  const Eigen::Vector2d &cameraCenter() const;
  const Eigen::Vector4d &backgroundColor() const;
  const Eigen::Vector4d &circleColor() const;
  const Eigen::Vector4d &fixedCircleColor() const;

  void setImageDimensions(const Vector2i &dims);
  void setCameraScale(const double &scale);
  void setCameraCenter(const Eigen::Vector2d &center);
  void setBackgroundColor(const Eigen::Vector4d &color);
  void setCircleColor(const Eigen::Vector4d &color);
  void setFixedCircleColor(const Eigen::Vector4d &color);

private:
  // Camera settings
  Vector2i m_image_dims;
  double m_cmra_scl;
  Eigen::Vector2d m_camera_center;
  Eigen::Vector4d m_bg_color;

  // Body render settings
  Eigen::Vector4d m_circle_color;
  Eigen::Vector4d m_fixed_circle_color;
};

#endif
