""" Tools for loading 2D MPM simulations """

import h5py
import numpy
import sys


class MPMState(object):
    """ A container for the state of a MPM simulation """
    def __init__(self, h5_file):
        if '/continuum/npoints' in h5_file:
            self.npoints = h5_file['/continuum/npoints'][0, 0]
            assert self.npoints >= 0
            self.q = h5_file['/continuum/q'][:]
            assert isinstance(self.q, numpy.ndarray)
            assert self.q.shape == (2, self.npoints)
            # self.v = h5_file['/continuum/v'][:]
            # assert isinstance( self.v, numpy.ndarray )
            # assert self.v.shape == self.q.shape
            # self.r = h5_file['/continuum/r'][:,0]
            # assert isinstance( self.r, numpy.ndarray )
            # assert self.r.shape == ( self.npoints, )
            # self.m = h5_file['/continuum/m'][:,0]
            # assert isinstance( self.m, numpy.ndarray )
            # assert self.m.shape == ( self.npoints, )
            # self.sigma = h5_file['/continuum/sigma'][:]
            # assert self.sigma.shape == ( 2, 2 * self.npoints )
            # self.periodic = h5_file['/continuum/periodic'][0,0] == 1
            # self.lees_edwards_offset = h5_file['/continuum/lees_edwards_offset'][0,0]
            # self.simulation_region = h5_file['/continuum/simulation_region'][:,0]
            # assert isinstance( self.simulation_region, numpy.ndarray )
            # assert self.simulation_region.shape == (4,)
        elif '/npoints' in h5_file:
            self.npoints = h5_file['/npoints'][0, 0]
            assert self.npoints >= 0
            self.q = h5_file['/q'][:]
            assert isinstance(self.q, numpy.ndarray)
            assert self.q.shape == (2, self.npoints)
        else:
            sys.exit('Error, HDF5 file does not appear to have continuum MPM data.')

#   def points( self ):
#     """ Iterates over the point positions, radii, and stresses. """
#     pnt_idx = 0
#     while pnt_idx < self.npoints:
#       yield self.q[:,pnt_idx], self.r[ pnt_idx ], self.sigma[ :, 2 * pnt_idx : 2 * pnt_idx + 2 ]
#       pnt_idx += 1


def load_configuration(file_name):
    """ Loads the continuum configuration from a MPM HDF5 config file """
    with h5py.File(file_name, 'r') as h5_file:
        mpm_state = MPMState(h5_file)
        if '/global_time' in h5_file:
            time = h5_file['/global_time'][0, 0]
        elif '/time' in h5_file:
            time = h5_file['/time'][0, 0]
        else:
            sys.exit('Error, failed to locate time when parsing mpm2d file.')
        assert time >= 0.0
        git_hash = h5_file['/git_hash'][0]
    return mpm_state, time, git_hash
