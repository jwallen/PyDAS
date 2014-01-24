import numpy as np
cimport numpy as np

cimport cython

cdef class DASSL:

    cdef public int maxOrder
    cdef public object tstop
    cdef public double initialStep
    cdef public double maximumStep
    cdef public object bandwidths
    cdef public bint nonnegative
    
    cdef public double t
    cdef public np.ndarray y
    cdef public np.ndarray dydt
    
    cdef np.ndarray info
    cdef np.ndarray atol
    cdef np.ndarray rtol
    cdef np.ndarray rwork
    cdef np.ndarray iwork
    cdef np.ndarray rpar
    cdef np.ndarray ipar
    cdef int idid
    
    cpdef initialize(self, double t0, np.ndarray y0, np.ndarray dydt0=?, atol=?, rtol=?)
    
    cpdef advance(self, double tout)
    
    cpdef step(self, double tout)
    
    cdef solve(self, double tout)