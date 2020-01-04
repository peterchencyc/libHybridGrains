#include "CairoRenderSettings.h"

CairoRenderSettings::CairoRenderSettings()
    : m_image_dims(1920, 1080), m_cmra_scl(10.0), m_camera_center(0.0, 0.0),
      m_bg_color(1.0, 1.0, 1.0, 1.0),
      m_circle_color(0.90196078431373, 0.16078431372549, 0.19607843137255, 1.0),
      m_fixed_circle_color(0.0, 0.0, 0.0, 1.0) {}

const Vector2i &CairoRenderSettings::imageDimensions() const {
  return m_image_dims;
}

const double &CairoRenderSettings::cameraScale() const { return m_cmra_scl; }

const Eigen::Vector2d &CairoRenderSettings::cameraCenter() const {
  return m_camera_center;
}

const Eigen::Vector4d &CairoRenderSettings::backgroundColor() const {
  return m_bg_color;
}

const Eigen::Vector4d &CairoRenderSettings::circleColor() const {
  return m_circle_color;
}

const Eigen::Vector4d &CairoRenderSettings::fixedCircleColor() const {
  return m_fixed_circle_color;
}

void CairoRenderSettings::setImageDimensions(const Vector2i &dims) {
  m_image_dims = dims;
}

void CairoRenderSettings::setCameraScale(const double &scale) {
  m_cmra_scl = scale;
}

void CairoRenderSettings::setCameraCenter(const Eigen::Vector2d &center) {
  m_camera_center = center;
}

void CairoRenderSettings::setBackgroundColor(const Eigen::Vector4d &color) {
  m_bg_color = color;
}

void CairoRenderSettings::setCircleColor(const Eigen::Vector4d &color) {
  m_circle_color = color;
}

void CairoRenderSettings::setFixedCircleColor(const Eigen::Vector4d &color) {
  m_fixed_circle_color = color;
}
