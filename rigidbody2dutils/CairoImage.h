#ifndef CAIRO_IMAGE_H
#define CAIRO_IMAGE_H

#include <cairo.h>

#include "scisim/Math/MathDefines.h"

class CairoImage final {

public:
  CairoImage(const Vector2i &image_dims, const Eigen::Vector4d &bg_color,
             const double &camera_scale, const Eigen::Vector2d &camera_center);

  CairoImage(const CairoImage &) = delete;
  CairoImage(CairoImage &&) = delete;
  CairoImage &operator=(const CairoImage &) = delete;
  CairoImage &operator=(CairoImage &&) = delete;

  ~CairoImage();

  void saveImage(const std::string &file_name);

  void drawFilledCircle(const Eigen::Vector2d &center, const double &radius,
                        const Eigen::Vector4d &color);

  void drawDrum(const Eigen::Vector2d &center, const double &drum_theta,
                const double &drum_radius, const Eigen::Vector4d &color);

private:
  void fillBackground(const Vector2i &image_dims,
                      const Eigen::Vector4d &bg_color);

  void setCamera(const Vector2i &image_dims, const double &camera_scale,
                 const Eigen::Vector2d &camera_center);

  cairo_surface_t *m_surface;
  cairo_t *m_cr;
};

#endif
