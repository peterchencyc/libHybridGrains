'''Tools for rendering 2D rigid body simulations.'''

import sys
import math
import cairo
import numpy

class CairoImage( object ):
  '''Wrapper to automatically save an image generated with Cairo.'''
  def __init__( self, image_name, image_size ):
    if image_name.lower().endswith( '.pdf' ):
      self.ims = cairo.PDFSurface( image_name, image_size[0], image_size[1] )
    elif image_name.lower().endswith( '.png' ):
      self.ims = cairo.ImageSurface( cairo.FORMAT_ARGB32, image_size[0], image_size[1] )
    self.image_name = image_name

  def __enter__( self ):
    return cairo.Context( self.ims )

  def __exit__( self, type, value, traceback ):
    if self.image_name.lower().endswith( '.pdf' ):
      self.ims.show_page()
    elif self.image_name.lower().endswith( '.png' ):
      self.ims.write_to_png( self.image_name )

class Saved( object ):
  '''Wrapper to automatically save and restore cairo transforms.'''
  def __init__( self, cr ):
    self.cr = cr

  def __enter__( self ):
    self.cr.save()
    return self.cr

  def __exit__( self, type, value, traceback ):
    self.cr.restore()

def render_white_background( cr, img_size ):
  '''Draws a white rectangle covering the entire image.'''
  cr.set_source_rgba( 1, 1, 1, 1 )
  cr.rectangle( 0, 0, img_size[0], img_size[1] )
  cr.fill()

def render_time( cr, time ):
  '''Draws the simulation time.'''
  # Draw a semi-transparent background for text
  cr.set_source_rgba( 0, 0, 0, 0.5 )
  cr.rectangle( 30 - 4, 40 - 16 + 30, 220, 20 )
  cr.fill()
  # Draw the simulation time
  cr.set_source_rgb( 1, 1, 1 )
  cr.select_font_face( "Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL )
  cr.set_font_size( 20.0 )
  cr.move_to( 30, 40 + 30 )
  cr.show_text( "Sim Time:  " + "{0:0.3f}".format( time ).zfill(6) )

def render_git_hash( cr, git_hash ):
  '''Draws the git hash.'''
  # Draw a semi-transparent background for text
  cr.set_source_rgba( 0, 0, 0, 0.5 )
  cr.rectangle( 30 - 4, 40 - 16, 550, 20 )
  cr.fill()
  # Draw the simulation time
  cr.set_source_rgb( 1, 1, 1 )
  cr.select_font_face( "Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL )
  cr.set_font_size( 20.0 )
  cr.move_to( 30, 40 )
  cr.show_text( "Git: " + git_hash )

def set_camera( cr, cmra_psn, img_size, cmra_scl ):
  '''Configure Cairo with the given camera settings.'''
  # Make the up direction a positive number
  cr.scale( 1.0, -1.0 )
  # Make 0,0 the center
  cr.translate( 0.5 * float( img_size[0] ), -0.5 * float( img_size[1] ) )
  # Rescale to requested zoom level
  cr.scale( float( img_size[1] ), float( img_size[1] ) )
  # Map to the camera space
  cr.scale( 0.5 / cmra_scl, 0.5 / cmra_scl )
  # Move to requested center
  cr.translate( -cmra_psn[0], -cmra_psn[1] )


def draw_circle(cr, x, r):
  '''Draws a circle of given radius at given position'''
  cr.arc(x[0], x[1], r, 0, 2 * math.pi)

# copy from mpm2d_rendering
def drawLine( cr, x, n, max_length, line_width ):
  """ Draws a line in Cairo """
  with Saved(cr):
    # Translate to center of the line
    cr.translate( x[0], x[1] )
    # Translate in the negative normal direction by half the line width
    cr.translate( - 0.5 * line_width * n[0], - 0.5 * line_width * n[1] )
    # Rotate to the desired orientation
    theta = - math.atan2( n[0], n[1] )
    cr.rotate( theta )
    # Set the endpoints
    cr.move_to( - 0.5 * max_length, 0.0 )
    cr.line_to(   0.5 * max_length, 0.0 )
    # Draw the line
    cr.stroke()


def render_discrete_bodies(cr, img_size, cmra_psn, cmra_scl, discrete, colors):
  '''Draws circular rigid bodies as outlined half circles.'''
  # body_line_width = 2.0 * 0.009765625
  # body_line_width = 2.0 * 0.09765625
  body_line_width = 4.0 * 0.02
  render_half_cirlces = False

  box_verts = numpy.matrix([[-1.0, -1.0, 1.0,  1.0],
                            [-1.0,  1.0, 1.0, -1.0]])

  with Saved(cr):
    set_camera( cr, cmra_psn, img_size, cmra_scl )

    cr.set_line_width( body_line_width )

    # Draw each body
    assert len(colors) == discrete.nbodies
    for bdy_idx, (x, theta, _, _, _, fixed, bdy_geo) in enumerate( discrete.bodies() ):
      if not fixed:
        cr.set_source_rgba(colors[bdy_idx][0], colors[bdy_idx][1], colors[bdy_idx][2], 1)
      else:
        # cr.set_source_rgba(0.4, 0.4, 0.4, 1)
        cr.set_source_rgba(0.0, 0.0, 0.0, 1)
      # Circle
      if bdy_geo[0] == 0:
        if not render_half_cirlces:
          radius = bdy_geo[1]
          draw_circle(cr, x, radius)
          cr.fill()
        else:
          radius = bdy_geo[1] - 0.5 * body_line_width
          draw_circle( cr, x, radius )
          cr.stroke()
          with Saved(cr):
            cr.translate( x[0], x[1] )
            cr.rotate( theta )
            cr.move_to( -radius, 0 )
            cr.line_to( radius, 0 )
            cr.stroke()
      # Box
      elif bdy_geo[0] == 1:
        rotation_mat = numpy.matrix( [ [ math.cos(theta), -math.sin(theta) ], [ math.sin(theta), math.cos(theta) ] ] )
        # assert fixed == 0
        # Grab the box half-widths and adjust by the half the line width
        half_width = bdy_geo[1]
        # half_width[0] -= 0.5 * body_line_width
        # half_width[1] -= 0.5 * body_line_width
        width_transform = numpy.matrix( [ [ half_width[0], 0 ], [ 0, half_width[1] ] ] )
        # Rescale the mesh to have the proper radius
        current_mesh = width_transform * numpy.matrix( box_verts )
        # Rotate the mesh to have the proper orientation
        current_mesh = rotation_mat * current_mesh
        # Translate the mesh to the proper center of mass
        current_mesh = numpy.transpose(current_mesh)
        current_mesh += x
        current_mesh = numpy.transpose(current_mesh)

        # cr.move_to(current_mesh[0,0], current_mesh[1,0] + 20)
        # cr.line_to(current_mesh[0,1], current_mesh[1,1] + 20)
        # cr.stroke()
        # cr.line_to(current_mesh[0,1], current_mesh[1,1] + 20)
        # cr.line_to(current_mesh[0,2], current_mesh[1,2] + 20)
        # cr.stroke()
        # cr.line_to(current_mesh[0,2], current_mesh[1,2] + 20)
        # cr.line_to(current_mesh[0,3], current_mesh[1,3] + 20)
        # cr.stroke()
        # cr.line_to(current_mesh[0,3], current_mesh[1,3] + 20)
        # cr.line_to(current_mesh[0,0], current_mesh[1,0] + 20)
        # cr.stroke()

        cr.move_to(current_mesh[0,0], current_mesh[1,0])
        cr.line_to(current_mesh[0,1], current_mesh[1,1])
        cr.line_to(current_mesh[0,2], current_mesh[1,2])
        cr.line_to(current_mesh[0,3], current_mesh[1,3])
        cr.close_path()
        cr.fill()
      else:
        sys.exit('Invalid geometry type encountered')


def render_circles(cr, img_size, cmra_psn, cmra_scl, circle_centers, circle_colors, radius):
    '''Draws circles.'''
    assert circle_centers.shape[0] == 2
    ncircles = circle_centers.shape[1]
    assert len(circle_colors) == ncircles
    with Saved(cr):
      set_camera( cr, cmra_psn, img_size, cmra_scl )
      for idx in range(0, ncircles):
          x = circle_centers[:, idx]
          color = circle_colors[idx]
          cr.set_source_rgba(color[0], color[1], color[2], 1.0)
          # radius = 0.06
          draw_circle(cr, x, radius)
          cr.fill()

def render_static_drums(cr, img_size, cmra_psn, cmra_scl, discrete):
    '''Draws static drums.'''
    body_line_width = 4.0 * 0.03
    with Saved(cr):
        set_camera( cr, cmra_psn, img_size, cmra_scl )
        cr.set_line_width( body_line_width )
        for drum_center, drum_radius, drum_theta in discrete.staticDrums():
            cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
            radius = drum_radius + 0.5 * body_line_width
            draw_circle( cr, drum_center, radius )
            cr.stroke()
            with Saved(cr):
                cr.translate( drum_center[0], drum_center[1] )
                cr.rotate( drum_theta )
                cr.move_to( -radius, 0 )
                cr.line_to( radius, 0 )
                cr.stroke()

def render_static_planes(cr, img_size, cmra_psn, cmra_scl, discrete):
    '''Draws static planes.'''
    body_line_width = 4.0 * 0.03
    max_length = 120
    with Saved(cr):
        set_camera( cr, cmra_psn, img_size, cmra_scl )
        cr.set_line_width( body_line_width )
        for x, n in discrete.staticPlanes():
            cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
            drawLine(cr, x, n, max_length, body_line_width)
            # cr.set_source_rgba(0.0, 0.0, 0.0, 1.0)
            # radius = drum_radius + 0.5 * body_line_width
            # draw_circle( cr, drum_center, radius )
            # cr.stroke()
            # with Saved(cr):
            #     cr.translate( drum_center[0], drum_center[1] )
            #     cr.rotate( drum_theta )
            #     cr.move_to( -radius, 0 )
            #     cr.line_to( radius, 0 )
            #     cr.stroke()

def render_grid_points(cr, img_size, cmra_psn, cmra_scl, circle_centers, circle_colors):
    '''Draws circles.'''
    ncircles = len(circle_centers)
    assert len(circle_colors) == ncircles
    with Saved(cr):
      set_camera( cr, cmra_psn, img_size, cmra_scl )
      for idx in range(0, ncircles):
          x = circle_centers[idx]
          color = circle_colors[idx]
          cr.set_source_rgba(color[0], color[1], color[2], color[3])
          radius = 0.1
          draw_circle(cr, x, radius)
          cr.fill()
