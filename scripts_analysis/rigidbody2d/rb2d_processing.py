'''Tools for loading and processing discrete simulation data.'''

import sys
import h5py
import numpy
import numpy.linalg


class DiscreteState(object):
    '''A container for the state of a discrete rigid body simulation.'''
    def __init__(self, h5_file):
        if '/discrete/q' in h5_file:
            prefix = '/discrete/'
        else:
            prefix = ''

        try:
            self.q = h5_file[prefix + 'q'][:, 0]
            assert isinstance(self.q, numpy.ndarray)
            assert self.q.shape[0] % 3 == 0
            self.nbodies = self.q.shape[0] / 3
            self.v = h5_file[prefix + 'v'][:, 0]
            assert isinstance(self.v, numpy.ndarray)
            assert self.v.shape == self.q.shape
            self.m = h5_file[prefix + 'm'][:, 0]
            assert isinstance(self.m, numpy.ndarray)
            assert self.m.shape == self.q.shape
            self.fixed = h5_file[prefix + 'kinematically_scripted'][:, 0]
            assert isinstance(self.fixed, numpy.ndarray)
            assert self.fixed.shape == (self.nbodies,)
            assert numpy.all(self.fixed >= 0)
            assert numpy.all(self.fixed <= 1)
            if prefix + 'unique_index' in h5_file:
                self.unique_ids = h5_file[prefix + 'unique_index'][:, 0]
                assert isinstance(self.unique_ids, numpy.ndarray)
                assert self.unique_ids.shape == (self.nbodies,)
            if prefix + 'grid_node_masses' in h5_file:
                self.grid_node_masses = h5_file[prefix + 'grid_node_masses'][:, 0]
            # Copy the geometry into a flat array
            self.geometry = []
            if prefix + 'geometry/circles' in h5_file:
                circle_structs = h5_file[prefix + 'geometry/circles'][:]
            else:
                circle_structs = numpy.empty([0])
            for circle_idx in range(0, len(circle_structs)):
                self.geometry.append((0, circle_structs[circle_idx][0]))
            if prefix + 'geometry/boxes' in h5_file:
                box_structs = h5_file[prefix + 'geometry/boxes'][:]
            else:
                box_structs = numpy.empty([0])
            for box_idx in range(0, len(box_structs)):
                self.geometry.append((1, box_structs[box_idx][0].tolist()))
            # Compute each body's flat geometry index
            geo_indices = h5_file[prefix + 'geometry/geometry_indices'][:]
            assert isinstance(geo_indices, numpy.ndarray)
            assert geo_indices.shape[0] == 2
            assert geo_indices.shape[1] == self.nbodies

            # This code will need to get updated if we import annuli
            geo_indices[0, :] = geo_indices[0, :] - 1

            def flatidx(bdy_idx):
                '''Converts a geometry type and instance index pair into a flat geometry index.'''
                return geo_indices[0, bdy_idx] * len(circle_structs) + geo_indices[1, bdy_idx]
            self.geometry_indices = [flatidx(bdy_idx) for bdy_idx in range(0, self.nbodies)]

            # Load in the drum geometry
            if prefix + 'static_geometry/static_drums' in h5_file:
                self.drums = h5_file[prefix + 'static_geometry/static_drums'][:]
            else:
                self.drums = numpy.empty([0])
            self.ndrums = len(self.drums)

            # Load in the static plane geometry
            if prefix + 'static_geometry/static_planes' in h5_file:
                self.planes = h5_file[prefix + 'static_geometry/static_planes'][:]
            else:
                self.planes = numpy.empty([0])
            self.nplanes = len(self.planes)

        except KeyError as key_exception:
            sys.exit('HDF5 Key Error: ' + key_exception.message)

    def bodies(self):
        '''Iterates over the configuration, per body.'''
        bdy_idx = 0
        while bdy_idx < self.nbodies:
            x = self.q[3 * bdy_idx: 3 * bdy_idx + 2]
            theta = self.q[3 * bdy_idx + 2]
            vel = self.v[3 * bdy_idx: 3 * bdy_idx + 2]
            omega = self.v[3 * bdy_idx + 2]
            assert self.m[3 * bdy_idx] == self.m[3 * bdy_idx + 1]
            mass = self.m[3 * bdy_idx]
            geo_idx = self.geometry_indices[bdy_idx]
            fixed = self.fixed[bdy_idx]
            yield x, theta, vel, omega, mass, fixed, self.geometry[geo_idx]
            bdy_idx += 1

    def bodiesWithIndex(self):
        '''Iterates over the configuration, per body.'''
        bdy_idx = 0
        while bdy_idx < self.nbodies:
            x = self.q[3 * bdy_idx: 3 * bdy_idx + 2]
            theta = self.q[3 * bdy_idx + 2]
            vel = self.v[3 * bdy_idx: 3 * bdy_idx + 2]
            omega = self.v[3 * bdy_idx + 2]
            assert self.m[3 * bdy_idx] == self.m[3 * bdy_idx + 1]
            mass = self.m[3 * bdy_idx]
            geo_idx = self.geometry_indices[bdy_idx]
            fixed = self.fixed[bdy_idx]
            yield x, theta, vel, omega, mass, fixed, self.geometry[geo_idx], geo_idx
            bdy_idx += 1

    def simulatedBoundingBox(self):
        '''Computes a tight axis-aligned bounding box over all bodies.'''
        minimum = numpy.array([numpy.PINF, numpy.PINF])
        maximum = numpy.array([numpy.NINF, numpy.NINF])
        for bdy_idx in range(0, self.nbodies):
            if not self.fixed[bdy_idx]:
                geo = self.geometry[self.geometry_indices[bdy_idx]]
                # Circle
                if geo[0] == 0:
                    r = geo[1]
                    x = self.q[3 * bdy_idx: 3 * bdy_idx + 2]
                    minimum = numpy.minimum(minimum, x - r)
                    maximum = numpy.maximum(maximum, x + r)
                else:
                    sys.exit('Geometry type {} not supported in simulatedBoundingBox.'.format(geo[0]))
        return minimum, maximum

    def staticDrums(self):
        '''Iterates over the configuration, per body.'''
        idx = 0
        while idx < self.ndrums:
            x = self.drums[idx][0]
            r = self.drums[idx][1]
            theta = self.drums[idx][2]
            yield x, r, theta
            idx += 1
    
    def staticPlanes(self):
        '''Iterates over the configuration, per body.'''
        idx = 0
        while idx < self.nplanes:
            x = self.planes[idx][0]
            n = self.planes[idx][1]
            yield x, n
            idx += 1


def load_discrete_state(file_name):
    '''Loads the discrete state from a hybrid HDF5 config file.'''
    try:
        with h5py.File(file_name, 'r') as h5_file:
            discrete_state = DiscreteState(h5_file)
            if 'time' in h5_file:
                time = h5_file['time'][0, 0]
            else:
                time = h5_file['global_time'][0, 0]
            assert time >= 0.0
            git_revision = h5_file['git_hash'][0]
    except KeyError as key_exception:
        sys.exit('HDF5 Key Error: ' + key_exception.message)
    except IOError as io_exception:
        sys.exit('Failed to open HDF5 file: ' + io_exception.message)
    return discrete_state, time, git_revision

# def computeKineticStressTensor( nbodies, m, v, area ):
#   assert v.shape[0] == 3 * nbodies
#   assert m.shape == v.shape
#   # Compute the average velocity
#   v_average = numpy.array( [ 0.0, 0.0 ] )
#   for body_idx in range( 0, nbodies ):
#     v_average += v[ 3 * body_idx : 3 * body_idx + 2 ]
#   v_average /= float(nbodies)
#
#   # Compute the kinetic contribution to the stress tensor
#   stress_kinetic = numpy.zeros( ( 2, 2 ) )
#   for body_idx in range( 0, nbodies ):
#     assert m[ 3 * body_idx ] == m[ 3 * body_idx + 1 ]
#     v_delta = v[ 3 * body_idx : 3 * body_idx + 2 ] - v_average
#     stress_kinetic += m[ 3 * body_idx ] * numpy.outer( v_delta, v_delta )
#
#   stress_kinetic /= area
#   return stress_kinetic

# def computeHomogenizedVerticalVelocity( nbodies, q, v, r, num_bins ):
#   # Hard coded bounding box to sample within
#   # TODO: Make this a parameter
#   bbox_start = numpy.array( [ -6.0, -6.0 ], dtype=numpy.float64 )
#   bbox_end = numpy.array( [ 6.0, 6.0 ], dtype=numpy.float64 )
#
#   # Compute the average diameter of the balls
#   average_diameter = 2.0 * numpy.mean( r )
#   # Compute a bin size
#   bucket_height = 8 * average_diameter
#
#   # Spacing between buckets
#   center_delta = ( bbox_end[1] - bbox_start[1] ) / ( num_bins - 1 )
#
#   # Generate the bin bounds
#   # Portions of box that were below the bottom and so wrapped around to the top
#   bin_above_wraparound = numpy.empty( ( num_bins, 2 ), dtype=numpy.float64 )
#   bin_above_wraparound[:,0].fill( bbox_end[1] )
#   bin_above_wraparound[:,1].fill( numpy.inf )
#   # Portions of the box that were above the top and so wrapped around to the bottom
#   bin_below_wraparound = numpy.empty( ( num_bins, 2 ), dtype=numpy.float64 )
#   bin_below_wraparound[:,0].fill( -numpy.inf )
#   bin_below_wraparound[:,1].fill( bbox_start[1] )
#   #print bin_below_wraparound
#   # Bounds on the normal part of the bin
#   bin_bounds = numpy.empty( ( num_bins, 2 ), dtype=numpy.float64 )
#   # Center of the bin
#   bin_centers = numpy.empty( num_bins, dtype=numpy.float64 )
#   for bin_idx in range( 0, num_bins ):
#     center = bbox_start[1] + bin_idx * center_delta
#     bin_centers[bin_idx] = center
#     lower = center - 0.5 * bucket_height
#     upper = center + 0.5 * bucket_height
#     # if lower < bbox_start[1]:
#     #   delta_below = lower - bbox_start[1]
#     #   bin_above_wraparound[bin_idx,1] = bbox_end[1] + delta_below
#     bin_bounds[bin_idx,0] = max( bbox_start[1], lower )
#     bin_bounds[bin_idx,1] = min( bbox_end[1], upper )
#     # if upper > bbox_end[1]:
#     #   delta_above = upper - bbox_end[1]
#     #   bin_below_wraparound[bin_idx,0] = bbox_start[1] + delta_above
#
#   # bin_above_wraparound: If less than first column and greater than second column
#   # bin_below_wraparound: If less than first column and greater than second column
#
#   # Compute average horizontal velocity in each bin
#   ball_count = numpy.zeros( num_bins, dtype=numpy.int64 )
#   total_velocities = numpy.zeros( num_bins, dtype=numpy.float64 )
#   # For each ball
#   for ball_idx in range( 0, nbodies ):
#     # Grab the center of mass of the ball and the ball's velocity
#     x_ball = q[ 3 * ball_idx : 3 * ball_idx + 2 ]
#     v_ball = v[ 3 * ball_idx : 3 * ball_idx + 2 ]
#     # TODO: Acceleration possible here
#     # For each bin
#     for bin_idx in range( 0, num_bins ):
#       # If the bin overlaps the bottom and wraps around to the top
#       if( x_ball[1] <= bin_above_wraparound[bin_idx,0] and x_ball[1] >= bin_above_wraparound[bin_idx,1] ):
#         ball_count[bin_idx] += 1
#         total_velocities[bin_idx] += v_ball[0] + 8.0
#       # If the ball falls into the 'normal' part of the bin
#       if( x_ball[1] >= bin_bounds[bin_idx,0] and x_ball[1] <= bin_bounds[bin_idx,1] ):
#         ball_count[bin_idx] += 1
#         total_velocities[bin_idx] += v_ball[0]
#       if( x_ball[1] <= bin_below_wraparound[bin_idx,0] and x_ball[1] >= bin_below_wraparound[bin_idx,1] ):
#         ball_count[bin_idx] += 1
#         total_velocities[bin_idx] += v_ball[0] - 8.0
#   average_velocities = total_velocities / ball_count
#
#   return bin_centers, average_velocities

# def computeHomogenizedVerticalDensity( nbodies, q, r, num_bins ):
#   # Hard coded bounding box to sample within
#   # TODO: Make this a parameter
#   bbox_start = numpy.array( [ -6.0, -6.0 ], dtype=numpy.float64 )
#   bbox_end = numpy.array( [ 6.0, 6.0 ], dtype=numpy.float64 )
#
#   # Compute the average diameter of the balls
#   average_diameter = 2.0 * numpy.mean( r )
#   # Compute a bin size
#   bucket_height = 8.0 * average_diameter
#   assert numpy.all( 2.0 * r < bucket_height )
#
#   # Spacing between buckets
#   center_delta = ( bbox_end[1] - bbox_start[1] ) / ( num_bins - 1 )
#
#   # Bounds on the normal part of the bin
#   bin_bounds = numpy.empty( ( num_bins, 2 ), dtype=numpy.float64 )
#
#   # Center of the bin
#   bin_centers = numpy.empty( num_bins, dtype=numpy.float64 )
#
#   for bin_idx in range( 0, num_bins ):
#     center = bbox_start[1] + bin_idx * center_delta
#     bin_centers[bin_idx] = center
#     lower = center - 0.5 * bucket_height
#     upper = center + 0.5 * bucket_height
#     bin_bounds[bin_idx,0] = max( bbox_start[1], lower )
#     bin_bounds[bin_idx,1] = min( bbox_end[1], upper )
#
#   # Compute the total mass in each bin
#   total_mass = numpy.zeros( num_bins, dtype=numpy.float64 )
#   # For each bin
#   for bin_idx in range( 0, num_bins ):
#     # Create a bin polygon
#     current_bin = numpy.array( [ [ bbox_start[0], bin_bounds[bin_idx,1] ], [ bbox_start[0], bin_bounds[bin_idx,0] ],
#                                  [ bbox_end[0], bin_bounds[bin_idx,0] ], [ bbox_end[0], bin_bounds[bin_idx,1] ] ] )
#     # For each body
#     for bdy_idx in range( 0, nbodies ):
#       # Grab the center of mass
#       x_body = q[ 3 * bdy_idx : 3 * bdy_idx + 2 ]
#       # Grab the radius
#       r_body = r[ bdy_idx ]
#       assert r_body > 0.0
#       # Generate a circle mesh
#       current_mesh = rb2d_geometry_processing.generateCircle( r_body, x_body )
#       intersecting_part = rb2d_geometry_processing.computePolygonBoxOverlap( current_mesh, current_bin )
#       total_mass[bin_idx] += rb2d_geometry_processing.polygonArea( intersecting_part )
#
#   # Compute densities from the total mass
#   for bin_idx in range( 0, num_bins ):
#     bin_area = ( bbox_end[0] - bbox_start[0] ) * ( bin_bounds[bin_idx,1] - bin_bounds[bin_idx,0] )
#     total_mass[bin_idx] /= bin_area
#
#   return bin_centers, total_mass

# def writeMathematicaArray( output_file, array_name, array ):
#   assert array.ndim == 1
#   # TODO: Verify that the output_file is open
#   output_file.write( array_name + ' = {' )
#   for i in range( 0, array.shape[0] - 1 ):
#     output_file.write( str( array[i] ) + ',')
#   output_file.write( str( array[-1] ) )
#   output_file.write( '};\n' )
