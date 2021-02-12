""" Generates an image of a 2D discrete rigid body SCISim simulation from an HDF5 configuration file. """

import os
import sys
import re
import argparse
import numpy
from rigidbody2d import rb2d_processing
from rigidbody2d import rb2d_rendering

# Zoomed out collapse: -d 1024 1024 -p 14 3.0 -z 15.0
# Zoomed in initial:   -d 1024 1024 -p 6.2273269459 3.0 -z 8.0

# Camera settings for the flow
# img_size = numpy.array( [ 1024, 1024 ] )
# cmra_psn = numpy.array( [ 6.2273269459, 3.0 ] )
# cmra_scl = 8.0

# Camera settings for column collapse
# img_size = numpy.array( [ 1024, 1024 ] )
# cmra_psn = numpy.array( [ 12.454653891710405 / 2.0, 12.454653891710405 / 4.0 ] )
# cmra_psn = numpy.array( [ 37.36396167513122 / 2.0, 12.454653891710405 / 4.0 ] )
# cmra_scl = 20.0

# Camera settings for the 2D Beverloo test, zoomed in settle
# img_size = numpy.array([1280, 720])
# cmra_psn = numpy.array([0.0, 31.0])
# cmra_scl = 12.0

# Camera settings for the 2D Beverloo test, zoomed out settle
# img_size = numpy.array([1920, 1080])
# cmra_psn = numpy.array([0.0, 40.0])
# cmra_scl = 25.0

# Camera settings for the 2D Beverloo test, zoomed out release
# img_size = numpy.array([1920, 1080])
# cmra_psn = numpy.array([0.0, 7.0])
# cmra_scl = 25.0
# # Tmp
# img_size = numpy.array([1920, 1080])
# cmra_psn = numpy.array([0.0, 15.0])
# cmra_scl = 17.0

# Defense: column collapse test
# img_size = numpy.array([1920, 1080])
# cmra_psn = numpy.array([13.5, 4.0])
# cmra_scl = 8.0
# render_headers = True

# Defense: Discrete drainage test
# img_size = numpy.array([1920, 1080])
# cmra_psn = numpy.array([0.0, 15.0])
# cmra_scl = 12.0
# render_headers = False
# render_background = True

# Defense: Bifurcation test
# img_size = numpy.array([1920, 1080])
# cmra_psn = numpy.array([0.0, 6.0])
# cmra_scl = 10.0
# render_headers = False
# render_background = True

# Thesis: Stake insertion
# img_size = numpy.array([1920, 1080])
# cmra_psn = numpy.array([6.2273269458552, 4.75])
# cmra_scl = 6.0
# render_headers = False
# render_background = True

# Defense: Small funnel drain
img_size = numpy.array([1920, 1080])
cmra_psn = numpy.array([0.0, 21.0])
cmra_scl = 23.0
render_headers = False
render_background = True

# Thesis: Drum
# img_size = numpy.array([1920, 1080])
# cmra_psn = numpy.array([0.0, 0.0])
# cmra_scl = 14.0
# render_headers = False
# render_background = True

# Defense: Box boucne test
# img_size = numpy.array([1920, 1080])
# cmra_psn = numpy.array([45.0, 25.0])
# cmra_scl = 30.0
# render_headers = False
# render_background = True

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
    dstate, time, git_revision = rb2d_processing.load_discrete_state(input_hdf5_name)
    dcolors = []
    for bdy_idx in range(0, dstate.nbodies):
        if dstate.fixed[bdy_idx]:
            dcolors.append([0., 0., 0.])
        else:
            # if dstate.unique_ids[bdy_idx] % 2 == 0:
            #     dcolors.append(numpy.array([70., 91., 181.]) / 255.)
            # else:
            #     dcolors.append(numpy.array([230., 41., 50.]) / 255.)
            # dcolors.append(numpy.array([70., 91., 181.]) / 255.) # Darker blue
            dcolors.append(0.9 * numpy.array([113., 172., 197.]) / 255.) # Blue
            # dcolors.append(numpy.array([147., 215., 250.]) / 255.)
            # dcolors.append(numpy.array([230., 41., 50.]) / 255.)
    with rb2d_rendering.CairoImage(output_image_name, image_size) as cr:
        if render_background:
            rb2d_rendering.render_white_background(cr, image_size)
        rb2d_rendering.render_static_drums(cr, image_size, camera_psn, camera_scl, dstate)
        rb2d_rendering.render_static_planes(cr, image_size, camera_psn, camera_scl, dstate)
        rb2d_rendering.render_discrete_bodies(cr, image_size, camera_psn, camera_scl, dstate, dcolors)
        if render_headers:
            rb2d_rendering.render_time(cr, time)
            rb2d_rendering.render_git_hash(cr, git_revision)


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
