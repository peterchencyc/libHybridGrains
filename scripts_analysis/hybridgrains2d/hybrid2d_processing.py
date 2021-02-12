'''Tools for loading 2D hybrid rigid body simulations computed with SCISim.'''

import sys
import math
import h5py
import numpy


class MPMState(object):
    '''A container for the state of a MPM simulation.'''
    def __init__(self, h5_file):
        try:
            self.npoints = h5_file['/continuum/npoints'][0, 0]
            assert self.npoints >= 0
            self.q = h5_file['/continuum/q'][:]
            assert isinstance(self.q, numpy.ndarray)
            assert self.q.shape == (2, self.npoints)
            self.v = h5_file['/continuum/v'][:]
            assert isinstance(self.v, numpy.ndarray)
            assert self.v.shape == self.q.shape
            self.r = h5_file['/continuum/r'][:, 0]
            assert isinstance(self.r, numpy.ndarray)
            assert self.r.shape == (self.npoints, )
            self.m = h5_file['/continuum/m'][:, 0]
            assert isinstance(self.m, numpy.ndarray)
            assert self.m.shape == (self.npoints, )
            self.sigma = h5_file['/continuum/sigma'][:]
            assert self.sigma.shape == (2, 2 * self.npoints)
            self.periodic = h5_file['/continuum/periodic'][0, 0] == 1
            self.lees_edwards_offset = h5_file['/continuum/lees_edwards_offset'][0, 0]
            self.simulation_region = h5_file['/continuum/simulation_region'][:, 0]
            self.region = h5_file['/continuum/region'][:, 0]
            self.cell_count = h5_file['/continuum/cell_count'][:, 0]
            assert isinstance(self.simulation_region, numpy.ndarray)
            assert self.simulation_region.shape == (4,)

            # For hybrid grains
            if '/continuum/hybrid_factors' in h5_file:
                self.hybrid_factors = h5_file['/continuum/hybrid_factors'][:, 0]
        except KeyError as key_exception:
            sys.exit('HDF5 Key Error: ' + key_exception.message)

    def cauchyStress(self, idx):
        '''Returns the Cauchy stress of a given point.'''
        return self.sigma[:, 2 * idx: 2 * idx + 2]

    def points(self):
        '''Iterates over the point locations, per material point.'''
        pnt_idx = 0
        while pnt_idx < self.npoints:
            yield self.q[:, pnt_idx]
            pnt_idx += 1

    def stresses(self):
        '''Iterates over the Cauchy stress, per material point.'''
        pnt_idx = 0
        while pnt_idx < self.npoints:
            yield self.sigma[:, 2 * pnt_idx: 2 * pnt_idx + 2]
            pnt_idx += 1

    def materialPoints(self):
        '''Iterates over the point positions, radii, and stresses.'''
        pnt_idx = 0
        while pnt_idx < self.npoints:
            yield self.q[:, pnt_idx], self.m[pnt_idx], self.r[pnt_idx], self.sigma[:, 2 * pnt_idx: 2 * pnt_idx + 2]
            pnt_idx += 1

    def totalMass(self):
        '''Sums the mass of all material points.'''
        return numpy.sum(self.m)

    def gridOrigin(self):
        '''Returns the lower left corner of the grid.'''
        return self.region[0:2]

    def cellCount(self):
        '''Returns the number of cells in the grid.'''
        return self.cell_count

    def cellSize(self):
        '''Returns the width of the grid cells.'''
        widths = (self.region[2:4] - self.region[0:2]) / self.cell_count.astype(self.region.dtype)
        assert abs(widths[1] - widths[0]) <= 1.0e-9
        return widths[0]

  # def totalMass( self ):
  #   """ Computes the total mass of the system """
  #   return numpy.sum( self.m )

  # def boundingBox( self ):
  #   """ Computes a tight axis-aligned bounding box on points with their radii. """
  #   minimum = numpy.array( [ numpy.PINF, numpy.PINF ] )
  #   maximum = numpy.array( [ numpy.NINF, numpy.NINF ] )
  #   for x, _, r, _ in self.points():
  #      minimum = numpy.minimum( minimum, x - r )
  #      maximum = numpy.maximum( maximum, x + r )
  #   return minimum, maximum

  # def computePointGrid( self, cell_size ):
  #   """ Computes a uniform grid of given cell width around the simulated geometry. """
  #   # Compute a bounding box around simulated geometry
  #   minimum, maximum = self.boundingBox()
  #   minimum -= 1.0e-6
  #   maximum += 1.0e-6
  #   # Compute the center of the bounding box
  #   grid_center = 0.5 * ( minimum + maximum )
  #   # Compute the number of cells
  #   cell_count = numpy.ceil( ( maximum - minimum ) / cell_size ).astype( int )
  #   # Compute the final grid width
  #   grid_width = cell_size * cell_count
  #   # Compute the lower left hand corner of the grid
  #   grid_start = grid_center - grid_width / 2.0
  #   return grid_start, cell_count


class HybridMPMState(object):
    '''A container for the state of a discrete rigid body simulation.'''
    def __init__(self, h5_file):
        try:
            self.grid_hybrid_weights = h5_file['/continuum/grid_node_hybrid_factors'][:, 0]
            self.grid_node_masses = h5_file['/continuum/grid_node_masses'][:, 0]
            assert len(self.grid_hybrid_weights) == len(self.grid_node_masses)
            self.grid_node_positions = h5_file['/continuum/grid_node_positions'][:, 0]
            assert len(self.grid_node_positions) == 2 * len(self.grid_node_masses)
            # assert isinstance(self.hybrid_weights, numpy.ndarray)
            self.nnodes = len(self.grid_hybrid_weights)
        except KeyError as key_exception:
            sys.exit('HDF5 Key Error: ' + key_exception.message)


class HybridDiscreteState(object):
    '''A container for the state of a discrete rigid body simulation.'''
    def __init__(self, h5_file):
        try:
            self.hybrid_weights = h5_file['/discrete/hybrid_factors'][:, 0]
            assert isinstance(self.hybrid_weights, numpy.ndarray)
            # self.q = h5_file['/discrete/q'][:,0]
            # assert isinstance( self.q, numpy.ndarray )
            # assert self.q.shape[0] % 3 == 0
            # self.nbodies = self.q.shape[0] / 3
            # self.v = h5_file['/discrete/v'][:,0]
            # assert isinstance( self.v, numpy.ndarray )
            # assert self.v.shape == self.q.shape
            # self.m = h5_file['/discrete/m'][:,0]
            # assert isinstance( self.m, numpy.ndarray )
            # assert self.m.shape == self.q.shape
            # self.fixed = h5_file['/discrete/kinematically_scripted'][:,0]
            # assert isinstance( self.fixed, numpy.ndarray )
            # assert self.fixed.shape == (self.nbodies,)
            # assert numpy.all( self.fixed >= 0 )
            # assert numpy.all( self.fixed <= 1 )
            # # Copy the geometry into a flat array
            # self.geometry = []
            # if '/discrete/geometry/circles' in h5_file:
            #     circle_structs = h5_file['/discrete/geometry/circles'][:]
            # else:
            #     circle_structs = numpy.empty([0])
            # for circle_idx in range(0,len(circle_structs)):
            #     self.geometry.append( ( 0, circle_structs[circle_idx][0] ) )
            # if '/discrete/geometry/boxes' in h5_file:
            #     box_structs = h5_file['/discrete/geometry/boxes'][:]
            # else:
            #     box_structs = numpy.empty([0])
            # for box_idx in range(0,len(box_structs)):
            #     self.geometry.append( ( 1, box_structs[box_idx][0].tolist() ) )
            # # Compute each body's flat geometry index
            # geo_indices = h5_file['/discrete/geometry/geometry_indices'][:]
            # assert isinstance( geo_indices, numpy.ndarray )
            # assert geo_indices.shape[0] == 2
            # assert geo_indices.shape[1] == self.nbodies
        except KeyError as key_exception:
            sys.exit('HDF5 Key Error: ' + key_exception.message)

    #     def flatidx( bdy_idx ):
    #       '''Converts a geometry type and instance index pair into a flat geometry index.'''
    #       return geo_indices[0,bdy_idx] * len(circle_structs) + geo_indices[1,bdy_idx]
    #     self.geometry_indices = [flatidx(bdy_idx) for bdy_idx in range(0,self.nbodies)]

    #     # For hybrid grains
    #     if '/discrete/hybrid_factors' in h5_file:
    #       self.hybrid_factors = h5_file['/discrete/hybrid_factors'][:,0]
    #     if '/discrete/always_fixed' in h5_file:
    #       self.always_fixed = h5_file['/discrete/always_fixed'][:,0]

    # def bodies( self ):
    #   '''Iterates over the configuration, per body.'''
    #   bdy_idx = 0
    #   while bdy_idx < self.nbodies:
    #     x = self.q[ 3 * bdy_idx : 3 * bdy_idx + 2 ]
    #     theta = self.q[ 3 * bdy_idx + 2 ]
    #     vel = self.v[ 3 * bdy_idx : 3 * bdy_idx + 2 ]
    #     omega = self.v[ 3 * bdy_idx + 2 ]
    #     assert self.m[ 3 * bdy_idx ] == self.m[ 3 * bdy_idx + 1 ]
    #     mass = self.m[ 3 * bdy_idx ]
    #     geo_idx = self.geometry_indices[bdy_idx]
    #     fixed = self.fixed[bdy_idx]
    #     yield x, theta, vel, omega, mass, fixed, geo_idx
    #     bdy_idx += 1

    # def centersOfMass( self ):
    #   '''Iterates over the centers of mass, per body.'''
    #   bdy_idx = 0
    #   while bdy_idx < self.nbodies:
    #     x = self.q[ 3 * bdy_idx : 3 * bdy_idx + 2 ]
    #     yield x
    #     bdy_idx += 1

    # def circleRadii( self ):
    #   '''If all bodies are circles, returns a list of their radii.'''
    #   radii = numpy.empty( (self.nbodies) )
    #   for bdy_idx in range(0,self.nbodies):
    #     geo_idx = self.geometry_indices[bdy_idx]
    #     geo = self.geometry[geo_idx]
    #     assert geo[0] == 0
    #     radii[bdy_idx] = geo[1]
    #   return radii

    # def totalKineticEnergy( self ):
    #   '''Computes the total kinetic energy of the system.'''
    #   kinetic_energy = 0.0
    #   for bdy_idx in range( 0, self.nbodies ):
    #     assert self.m[ 3 * bdy_idx ] == self.m[ 3 * bdy_idx + 1 ]
    #     total_mass = self.m[ 3 * bdy_idx ]
    #     assert total_mass > 0.0
    #     inertia = self.m[ 3 * bdy_idx + 2 ]
    #     assert inertia > 0.0
    #     cm_vel = self.v[ 3 * bdy_idx : 3 * bdy_idx + 2 ]
    #     omega = self.v[ 3 * bdy_idx + 2 ]
    #     # Translation component
    #     kinetic_energy += total_mass * numpy.dot( cm_vel, cm_vel )
    #     # Rotational component
    #     kinetic_energy += inertia * omega * omega
    #   return 0.5 * kinetic_energy

    # def rotationalKineticEnergy( self ):
    #   '''Computes the rotational component of the kinetic energy of the system.'''
    #   kinetic_energy = 0.0
    #   for bdy_idx in range( 0, self.nbodies ):
    #     inertia = self.m[ 3 * bdy_idx + 2 ]
    #     assert inertia > 0.0
    #     omega = self.v[ 3 * bdy_idx + 2 ]
    #     kinetic_energy += inertia * omega * omega
    #   return 0.5 * kinetic_energy

    # def translationalKineticEnergy( self ):
    #   '''Computes the translational component of the kinetic energy of the system.'''
    #   kinetic_energy = 0.0
    #   for bdy_idx in range( 0, self.nbodies ):
    #     assert self.m[ 3 * bdy_idx ] == self.m[ 3 * bdy_idx + 1 ]
    #     total_mass = self.m[ 3 * bdy_idx ]
    #     assert total_mass > 0.0
    #     cm_vel = self.v[ 3 * bdy_idx : 3 * bdy_idx + 2 ]
    #     kinetic_energy += total_mass * numpy.dot( cm_vel, cm_vel )
    #   return 0.5 * kinetic_energy

    # def isfixed( self, bdy_idx ):
    #   '''Returns true if a given body is fixed, false otherwise.'''
    #   assert bdy_idx >= 0
    #   assert bdy_idx < self.nbodies
    #   return self.fixed[ bdy_idx ] == 1

    # def is_always_fixed( self, bdy_idx ):
    #   '''Returns true if a given body is always fixed, false otherwise.'''
    #   assert bdy_idx >= 0
    #   assert bdy_idx < self.nbodies
    #   return self.always_fixed[ bdy_idx ] == 1

    # def x( self, bdy_idx ):
    #   '''Returns the center of mass of a body.'''
    #   assert bdy_idx >= 0
    #   assert bdy_idx < self.nbodies
    #   return self.q[ 3 * bdy_idx : 3 * bdy_idx + 2 ]

    # def bodyArea( self, geo_idx ):
    #   '''Returns the area of a body.'''
    #   assert geo_idx >= 0
    #   assert geo_idx < len(self.geometry)
    #   geo = self.geometry[geo_idx]
    #   if geo[0] == 0:
    #     return math.pi * geo[1] * geo[1]
    #   else:
    #     sys.exit('Unsupported body type in bodyArea')

    # def centerOfMassBoundingBox( self ):
    #   '''Computes a tight axis-aligned bounding box over all centers of mass.'''
    #   minimum = numpy.array( [ numpy.PINF, numpy.PINF ] )
    #   maximum = numpy.array( [ numpy.NINF, numpy.NINF ] )
    #   for bdy_idx in range(0,self.nbodies):
    #     cm = self.x( bdy_idx )
    #     minimum = numpy.minimum( minimum, cm )
    #     maximum = numpy.maximum( maximum, cm )
    #   return minimum, maximum


# class DiscreteForces( object ):
#   '''A container for the forces of a discrete rigid body simulation.'''
#   def __init__( self, h5_file ):
#     try:
#       self.collision_count = h5_file['collision_count'][0,0]
#       if self.collision_count == 0:
#         self.collision_points = numpy.empty( ( 2, 0 ), numpy.float64 )
#         self.collision_indices = numpy.empty( ( 2, 0 ), numpy.int32 )
#         self.collision_forces = numpy.empty( ( 2, 0 ), numpy.float64 )
#         self.collision_normals = numpy.empty( ( 2, 0 ), numpy.float64 )
#       else:
#         self.collision_points = h5_file['collision_points'][:]
#         self.collision_indices = h5_file['collision_indices'][:]
#         self.collision_forces = h5_file['collision_forces'][:]
#         self.collision_normals = h5_file['collision_normals'][:]
#       assert isinstance( self.collision_points, numpy.ndarray )
#       assert self.collision_points.ndim == 2
#       assert self.collision_points.shape[0] == 2
#       assert self.collision_points.shape[1] == self.collision_count
#       assert isinstance( self.collision_indices, numpy.ndarray )
#       assert self.collision_indices.shape == self.collision_points.shape
#       assert isinstance( self.collision_forces, numpy.ndarray )
#       assert self.collision_forces.shape == self.collision_points.shape
#       assert isinstance( self.collision_normals, numpy.ndarray )
#       assert self.collision_normals.shape == self.collision_points.shape
#     except KeyError as key_exception:
#       sys.exit('HDF5 Key Error: ' + key_exception.message)

#   def points( self ):
#     '''Iterates over the collision points, per collision.'''
#     clsn_idx = 0
#     while clsn_idx < self.collision_count:
#       yield self.collision_points[:,clsn_idx]
#       clsn_idx += 1

#   def impulses( self ):
#     '''Iterates over the impulses, per collision.'''
#     clsn_idx = 0
#     while clsn_idx < self.collision_count:
#       yield self.collision_forces[:,clsn_idx]
#       clsn_idx += 1

#   def collisions( self ):
#     '''Iterates over the collision state, per collision.'''
#     clsn_idx = 0
#     while clsn_idx < self.collision_count:
#       yield self.collision_indices[:,clsn_idx], self.collision_points[:,clsn_idx], self.collision_forces[:,clsn_idx], self.collision_normals[:,clsn_idx]
#       clsn_idx += 1


# class GridContainer( object ):
#   '''A generic uniform grid.'''
#   def __init__( self, origin, cell_size, cell_count ):
#     self.origin = origin
#     assert self.origin.shape == (2,)
#     self.cell_size = cell_size
#     assert self.cell_size > 0.0
#     self.cell_count = cell_count
#     assert self.cell_count.shape == (2,)
#     self.cell_contents = {}

#   def addEmptyListContents( self, content_name ):
#     '''Creates a new property that lives on the grid.'''
#     self.cell_contents[content_name] = [ [] for _ in range(0,self.cell_count[0] * self.cell_count[1]) ]

#   def addDefaultContents( self, content_name, default_value ):
#     '''Creates a new property that lives on the grid.'''
#     self.cell_contents[content_name] = [ default_value for _ in range(0,self.cell_count[0] * self.cell_count[1]) ]

#   def hasContents( self, content_name ):
#     '''True if the gird has the given content.'''
#     return content_name in self.cell_contents

#   def contents( self, content_name, i, j ):
#     '''Returns the contents of given name in the i, jth cell.'''
#     assert i < self.cell_count[0]
#     assert j < self.cell_count[1]
#     return self.cell_contents[content_name][i * self.cell_count[1] + j]

#   def setContents( self, content_name, i, j, value ):
#     '''Sets the contents of given name in the i, jth cell.'''
#     assert i < self.cell_count[0]
#     assert j < self.cell_count[1]
#     self.cell_contents[content_name][i * self.cell_count[1] + j] = value

#   def containingCell( self, x ):
#     '''Determines the cell that contains the point x.'''
#     return numpy.floor( ( x - self.origin ) / self.cell_size ).astype(int)

#   def cellIsValid( self, i, j ):
#     '''True if the given cell indices are valid.'''
#     if i < 0 or j < 0 or i >= self.cell_count[0] or j >= self.cell_count[1]:
#       return False
#     else:
#       return True

#   def cellExtent( self, i, j ):
#     '''Returns the lower left and upper right corner of a cell.'''
#     lower_left = self.origin + numpy.array([i,j]) * self.cell_size
#     upper_right = self.origin + numpy.array([i + 1,j + 1]) * self.cell_size
#     return lower_left, upper_right

#   def cellArea( self ):
#     '''Returns the area of a cell.'''
#     return self.cell_size * self.cell_size

#   def cells( self ):
#     '''Iterates over the cells.'''
#     for idx0 in range(0,self.cell_count[0]):
#       for idx1 in range(0,self.cell_count[1]):
#         yield idx0, idx1


# class HybridData( object ):
#   '''A container for the hybrid data.'''
#   def __init__( self, h5_file ):
#     try:
#       self.coupling_stiffness = h5_file['/hybrid/coupling_stiffness']
#       self.simulation_type = h5_file['/hybrid/sim_type']
#       # self.level_set = h5_file['/hybrid/level_set'][:]
#     except KeyError as key_exception:
#       sys.exit('HDF5 Key Error: ' + key_exception.message)


def load_hybrid_configuration(file_name):
    '''Loads the discrete and continuum configuration from a hybrid HDF5 config file.'''
    try:
        with h5py.File(file_name, 'r') as h5_file:
            # mpm_state = MPMState( h5_file )
            hybrid_discrete_state = HybridDiscreteState(h5_file)
            hybrid_mpm_state = HybridMPMState(h5_file)
            # time = h5_file['global_time'][0,0]
            # assert time >= 0.0
            # git_revision = h5_file['git_hash'][0]
            # hybrid_data = HybridData( h5_file )
    except KeyError as key_exception:
        sys.exit('HDF5 Key Error in load_hybrid_configuration: ' + key_exception.message)
    except IOError as io_exception:
        sys.exit('Failed to open HDF5 file in load_hybrid_configuration: ' + io_exception.message)
    # return discrete_state, mpm_state, hybrid_data, time, git_revision
    return hybrid_discrete_state, hybrid_mpm_state


# def load_lm_impact( file_name ):
#   '''Loads the velocities from a hybrid HDF5 config file.'''
#   try:
#     with h5py.File( file_name, 'r' ) as h5_file:
#       mpm_post_v = h5_file['LMCoupling/mpm_post_v'][:]
#       mpm_coupled_v = h5_file['LMCoupling/mpm_coupled_v'][:]
#       dem_post_v = h5_file['LMCoupling/dem_post_v'][:]
#       dem_coupled_v = h5_file['LMCoupling/dem_coupled_v'][:]
#       time = h5_file['global_time'][0,0]
#       assert time >= 0.0
#       git_revision = h5_file['git_hash'][0]
#   except KeyError as key_exception:
#     # sys.exit('HDF5 Key Error in load_lm_impact: ' + key_exception.message)
#     with h5py.File( file_name, 'r' ) as h5_file:
#       mpm_post_v = h5_file['/continuum/v'][:]
#       mpm_post_v = [flatten for inner in mpm_post_v for flatten in inner]
#       mpm_coupled_v = h5_file['/continuum/v'][:]
#       mpm_coupled_v = [flatten for inner in mpm_coupled_v for flatten in inner]
#       dem_post_v = h5_file['/discrete/v'][:,0]
#       dem_coupled_v = h5_file['/discrete/v'][:,0]
#       time = h5_file['global_time'][0,0]
#       git_revision = h5_file['git_hash'][0]
#   except IOError as io_exception:
#     # sys.exit('Failed to open HDF5 file in load_lm_impact: ' + io_exception.message)
#     with h5py.File( file_name, 'r' ) as h5_file:
#       mpm_post_v = h5_file['/continuum/v'][:]
#       mpm_post_v = [flatten for inner in mpm_post_v for flatten in inner]
#       mpm_coupled_v = h5_file['/continuum/v'][:]
#       mpm_coupled_v = [flatten for inner in mpm_coupled_v for flatten in inner]
#       dem_post_v = h5_file['/discrete/v'][:,0]
#       dem_coupled_v = h5_file['/discrete/v'][:,0]
#       time = h5_file['global_time'][0,0]
#       git_revision = h5_file['git_hash'][0]
#   return mpm_post_v, mpm_coupled_v, dem_post_v, dem_coupled_v, time, git_revision

# def load_discrete_state( file_name ):
#   '''Loads the discrete state from a hybrid HDF5 config file.'''
#   try:
#     with h5py.File( file_name, 'r' ) as h5_file:
#       discrete_state = DiscreteState( h5_file )
#       if 'global_time' in h5_file:
#         time = h5_file['global_time'][0,0]
#       else:
#         time = h5_file['time'][0,0]
#       # time = h5_file['time'][0,0]
#       assert time >= 0.0
#       git_revision = h5_file['git_hash'][0]
#   except KeyError as key_exception:
#     sys.exit('HDF5 Key Error: ' + key_exception.message)
#   except IOError as io_exception:
#     sys.exit('Failed to open HDF5 file: ' + io_exception.message)
#   return discrete_state, time, git_revision


# def load_mpm_state(file_name):
#     '''Loads the MPM state from a hybrid HDF5 config file.'''
#     try:
#         with h5py.File(file_name, 'r') as h5_file:
#             mpm_state = MPMState(h5_file)
#             time = h5_file['global_time'][0, 0]
#             assert time >= 0.0
#             git_revision = h5_file['git_hash'][0]
#     except KeyError as key_exception:
#         sys.exit('HDF5 Key Error: ' + key_exception.message)
#     except IOError as io_exception:
#         sys.exit('Failed to open HDF5 file: ' + io_exception.message)
#     return mpm_state, time, git_revision


# def load_discrete_forces( file_name ):
#   '''Loads the discrete forces from a hybrid HDF5 force file.'''
#   try:
#     with h5py.File( file_name, 'r' ) as h5_file:
#       discrete_forces = DiscreteForces( h5_file )
#       # time = h5_file['global_time'][0,0]
#       time = h5_file['global_time'][0,0]
#       assert time >= 0.0
#       # time_step = h5_file['discrete_timestep'][0,0]
#       time_step = h5_file['discrete_timestep'][0,0]
#       assert time_step > 0.0
#       git_revision = h5_file['git_hash'][0]
#   except KeyError as key_exception:
#     sys.exit('HDF5 Key Error for {} in load_discrete_forces: {}'.format(file_name,key_exception.message))
#   except IOError as io_exception:
#     sys.exit('Failed to open {} in load_discrete_forces: {}'.format(file_name,io_exception.message))
#   return discrete_forces, time, time_step, git_revision
