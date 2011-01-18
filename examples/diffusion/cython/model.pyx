#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

"""

import numpy
cimport numpy

from pydas cimport DASSL

cdef class DiffusionModel(DASSL):
    """
    An extension type for solving the diffusion equation in a 1D thin membrane. 
    The attribute `N` is the number of grid points to use to discretize the
    spatial direction.    
    """
    
    cdef public int N
    
    def __init__(self, N=10):
         self.N = N
    
    def residual(self, double t, numpy.ndarray[numpy.float64_t, ndim=1] y, numpy.ndarray[numpy.float64_t, ndim=1] dydt):
        
        cdef Py_ssize_t i
        cdef double dx
        cdef numpy.ndarray[numpy.float64_t, ndim=1] delta
        
        dx = 1.0 / (self.N - 1)
        
        delta = numpy.zeros(y.shape[0], numpy.float64)
        # Internal nodes
        for i in range(1, self.N-1):
            delta[i] = (y[i+1] - 2 * y[i] + y[i-1]) / (dx * dx) - dydt[i]
        # Left boundary (x = 0)
        i = 0
        delta[i] = y[i] - 1.0 
        # Right boundary (x = 1)
        i = self.N - 1
        delta[i] = y[i]
        
        return delta, 0
