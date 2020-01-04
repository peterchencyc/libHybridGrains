#include "CairoImage.h"

CairoImage::CairoImage(const Vector2i &image_dims,
                       const Eigen::Vector4d &bg_color,
                       const double &camera_scale,
                       const Eigen::Vector2d &camera_center)
    : m_surface(cairo_image_surface_create(CAIRO_FORMAT_ARGB32, image_dims.x(),
                                           image_dims.y())),
      m_cr(cairo_create(m_surface)) {
  assert(m_surface != nullptr);
  assert(m_cr != nullptr);
  fillBackground(image_dims, bg_color);
  setCamera(image_dims, camera_scale, camera_center);
}

CairoImage::~CairoImage() {
  cairo_destroy(m_cr);
  cairo_surface_destroy(m_surface);
}

void CairoImage::saveImage(const std::string &file_name) {
  // TODO: Check the result of the save
  cairo_surface_write_to_png(m_surface, file_name.c_str());
}

void CairoImage::drawFilledCircle(const Eigen::Vector2d &center,
                                  const double &radius,
                                  const Eigen::Vector4d &color) {
  cairo_set_source_rgba(m_cr, color(0), color(1), color(2), color(3));
  cairo_arc(m_cr, center.x(), center.y(), radius, 0.0,
            2.0 * MathDefines::PI<double>());
  cairo_fill(m_cr);
}

void CairoImage::drawDrum(const Eigen::Vector2d &center,
                          const double &drum_theta, const double &drum_radius,
                          const Eigen::Vector4d &color) {
  const double body_line_width = 4.0 * 0.03;
  const double radius = drum_radius + 0.5 * body_line_width;

  cairo_set_line_width(m_cr, body_line_width);
  cairo_set_source_rgba(m_cr, color(0), color(1), color(2), color(3));

  cairo_arc(m_cr, center.x(), center.y(), radius, 0.0,
            2.0 * MathDefines::PI<double>());

  cairo_stroke(m_cr);

  cairo_save(m_cr);
  cairo_translate(m_cr, center.x(), center.y());
  cairo_rotate(m_cr, drum_theta);
  cairo_move_to(m_cr, -radius, 0);
  cairo_line_to(m_cr, radius, 0);
  cairo_stroke(m_cr);
  cairo_restore(m_cr);
}

void CairoImage::fillBackground(const Vector2i &image_dims,
                                const Eigen::Vector4d &bg_color) {
  cairo_set_source_rgba(m_cr, bg_color(0), bg_color(1), bg_color(2),
                        bg_color(3));
  cairo_rectangle(m_cr, 0, 0, image_dims.x(), image_dims.y());
  cairo_fill(m_cr);
}

void CairoImage::setCamera(const Vector2i &image_dims,
                           const double &camera_scale,
                           const Eigen::Vector2d &camera_center) {
  // Make the up direction a positive number
  cairo_scale(m_cr, 1.0, -1.0);
  // Make 0, 0 the center
  cairo_translate(m_cr, 0.5 * double(image_dims.x()),
                  -0.5 * double(image_dims.y()));
  // Rescale to requested zoom level
  cairo_scale(m_cr, double(image_dims.y()), double(image_dims.y()));
  // Map to the camera space
  cairo_scale(m_cr, 0.5 / camera_scale, 0.5 / camera_scale);
  // Move to requested center
  cairo_translate(m_cr, -camera_center.x(), -camera_center.y());
}
