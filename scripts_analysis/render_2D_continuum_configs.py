""" Generates an image of a 2D discrete rigid body SCISim simulation from an HDF5 configuration file. """

import os
import sys
import re
import argparse
import numpy
from rigidbody2d import rb2d_processing
from rigidbody2d import rb2d_rendering
from mpm2d import mpm2d_processing
from hybridgrains2d import hybrid2d_processing

# # Camera settings for the 2D Beverloo test, zoomed out release
# img_size = numpy.array([1920, 1080])
# cmra_psn = numpy.array([0.0, 7.0])
# cmra_scl = 25.0

# # Tmp
# img_size = numpy.array([1920, 1080])
# cmra_psn = numpy.array([0.0, 15.0])
# cmra_scl = 17.0

# Defense: column collapse test
img_size = numpy.array([1920, 1080])
cmra_psn = numpy.array([13.5, 4.0])
cmra_scl = 8.0
circle_rad = 0.05

# big hourglass
# img_size = numpy.array([1920, 1080])
# cmra_psn = numpy.array([0.0, 21.0])
# cmra_scl = 43.0
# render_headers = False
# render_background = True


parser = argparse.ArgumentParser( description='Renders images from 2D rigid body HDF5 files.' )
parser.add_argument( '-x', '--skip_existing', help='do not generate a plot if a file already exists', action='store_true' )
parser.add_argument( '-i', metavar='input_directory', type=str, nargs=1, help='input directory name', required=True )
parser.add_argument( '-o', metavar='output_directory', type=str, nargs=1, help='output directory name', required=True )
parser.add_argument( '-b', metavar='begin_file_num', type=int, nargs=1, help='number of the configuration file to begin processing from', required=False )
parser.add_argument( '-e', metavar='end_file_num', type=int, nargs=1, help='number of the configuration file to end processing at', required=False )
parser.add_argument( '-s', metavar='file_spacing', type=int, nargs=1, help='spacing between files to process', required=False )
args = parser.parse_args()

# Verify that the input directory exists
input_directory = args.i[0]
if not os.path.isdir( input_directory ):
  print 'Input directory', input_directory, 'does not exist'
  sys.exit(1)

# Verify that the output directory exists
output_directory = args.o[0]
if not os.path.isdir( output_directory ):
  print 'Output directory', output_directory, 'does not exist'
  sys.exit(1)

# Begin and end must be specified together
if ( args.b is None ) != ( args.e is None ):
  print 'Error, arguments b and e must be set together'
  sys.exit(1)

# Grab the start and end file numbers from the user provided arguments, if provided
range_set = args.b is not None and args.e is not None and args.s is not None
if range_set:
  start_file_num = args.b[0]
  if start_file_num < 0:
    print "Error, the begining file number must be a non-negative integer. Exiting."
    sys.exit(1)

  end_file_num = args.e[0]
  if start_file_num > end_file_num:
    print "Error, the begining file number must be less than the end file number. Exiting."
    sys.exit(1)

  file_skip = args.s[0]
  if file_skip < 1:
    print "Error, the spacing between files must be a positive integer. Exiting."
    sys.exit(1)

# Skip images that already exist
skip_existing = args.skip_existing

# Create a regular expression for matching config file names and for extracting the file number
# config_file_pattern = r'.+(\d+)\.h5'
config_file_pattern = r'config_(\d+)\.h5'
config_file_matcher = re.compile(config_file_pattern)


def draw_config(input_hdf5_name, output_image_name, image_size, camera_psn, camera_scl):
    '''Draws discrete bodies.'''
    # dstate, time, git_revision = rb2d_processing.load_discrete_state(input_hdf5_name)
    mstate, mtime, mhash = mpm2d_processing.load_configuration(input_hdf5_name)
    # hdstate, hmstate = hybrid2d_processing.load_hybrid_configuration(input_hdf5_name)

    # dcolors = []
    # for bdy_idx in range(0, dstate.nbodies):
    #     if dstate.fixed[bdy_idx]:
    #         dcolors.append([0., 0., 0.])
    #     elif hdstate.hybrid_weights[bdy_idx] != 1.0:
    #         if dstate.unique_ids[bdy_idx] % 2 == 0:
    #             dcolors.append(numpy.array([117, 13, 255]) / 255.)
    #         else:
    #             dcolors.append(numpy.array([255, 133, 28]) / 255.)
    #     else:
    #         if dstate.unique_ids[bdy_idx] % 2 == 0:
    #             dcolors.append(numpy.array([70., 91., 181.]) / 255.)
    #         else:
    #             dcolors.append(numpy.array([230., 41., 50.]) / 255.)

    circle_colors = []
    for bdy_idx in range(0, mstate.npoints):
        circle_colors.append([230. / 255., 41. / 255., 50. / 255., 1.0])

    # For plotting the grid
    # grid_point_colors = []
    # grid_points_to_plot = []
    # for node_idx in range(0, hmstate.nnodes):
    #     if hmstate.grid_hybrid_weights[node_idx] != 1.0:
    #         grid_point_colors.append([0., 0., 0., 1.0])
    #         grid_points_to_plot.append(hmstate.grid_node_positions[2 * node_idx: 2 * node_idx + 2])

    with rb2d_rendering.CairoImage(output_image_name, image_size) as cr:
        rb2d_rendering.render_white_background(cr, image_size)
        # rb2d_rendering.render_discrete_bodies(cr, image_size, camera_psn, camera_scl, dstate, dcolors)
        rb2d_rendering.render_circles(cr, image_size, camera_psn, camera_scl, mstate.q, circle_colors, circle_rad)
        # rb2d_rendering.render_grid_points(cr, image_size, camera_psn, camera_scl, grid_points_to_plot, grid_point_colors)
        rb2d_rendering.render_time(cr, mtime)
        rb2d_rendering.render_git_hash(cr, mhash)


for file_name in sorted(os.listdir(input_directory)):
    # Determine if this is a configuration file
    config_file_match = config_file_matcher.match(file_name)
    if config_file_match is None:
        continue

    # If a range was specified, check that this file is in the desired range
    if range_set:
        file_number = int(config_file_match.group(1))
        assert file_number >= 0
        # Skip this file if it is not in the valid file number range
        if file_number < start_file_num or file_number > end_file_num:
            continue
        if (file_number - start_file_num) % file_skip != 0:
            continue

    # Grab the base filename
    split_arg = os.path.splitext(file_name)
    assert split_arg[1] == '.h5'
    output_file_name = split_arg[0] + '.png'
    full_output_path = os.path.join(output_directory, output_file_name)
    # Get the full input file path
    input_file_name = os.path.join(input_directory, file_name)
    if skip_existing and os.path.exists(full_output_path):
        print output_file_name, "exists, skipping", input_file_name
        continue

    print "Processing", input_file_name
    draw_config(input_file_name, full_output_path, img_size, cmra_psn, cmra_scl)
    print "Generated", full_output_path, "from", input_file_name
