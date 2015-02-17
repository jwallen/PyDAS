################################################################################
#
#   PyDAS - A Python wrapper for several differential algebraic system solvers
#
#   Copyright (c) 2010-2014 by Joshua W. Allen (joshua.w.allen@gmail.com) 
#                           extended by Connie W. Gao (connieg@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This `Cython <http://www.cython.org/>`_ module exposes the DASSL differential 
algebraic system solver to Python and provides a Python extension type, the 
:class:`DASSL` base class, for working with DASSL.

To use DASSL, write a Python class or Cython extension type that derives from
the :class:`DASSL` class and implement the :meth:`residual`
method, which accepts :math:`t`, :math:`\\mathbf{y}`, and 
:math:`\\mathbf{\\frac{d \\mathbf{y}}{dt}}` as arguments and returns the
corresponding value of :math:`\\mathbf{g} \\left(t, \\mathbf{y}, \\frac{d \\mathbf{y}}{dt} \\right)`.
Run by calling the :meth:`initialize` method to set the initial conditions
and solver tolerances, then by using the :meth:`advance` or :meth:`step`
methods to move forward or backward in the independent variable. If you 
know the form of the Jacobian, you can implement the :meth:`jacobian`
method to provide it to the solver.

You can implement your derived class in pure Python, but for a significant
speed boost consider using Cython to compile the module. You can see the
proper Cython syntax for the residual and jacobian methods by looking at the
corresponding methods in the :class:`DASSL` base class.
"""

import numpy as np
cimport numpy as np

cimport cython

################################################################################

# Expose the (double-precision) DASSL function
cdef extern from "dassl.h":
    int ddassl_(
        void* res,      # The residual function that defines the ODE/DAE system
        int* neq,       # The number of equations to be solved
        double* t,      # The current value of the independent variable
        double* y,      # The current values of the dependent variables
        double* yprime, # The current values of the first derivatives of the dependent variables
        double* tout,   # The value of the independent variable at which a solution is desired
        int* info,      # Parameters controlling how the integration is performed
        double* rtol,   # The relative error tolerance(s), either as a scalar or a vector
        double* atol,   # The absolute error tolerance(s), either as a scalar or a vector
        int* idid,      # Report of the solver actions, used to control subsequent calls
        double* rwork,  # Work space for double-precision values
        int* lrw,       # The length of the double-precision workspace
        int* iwork,     # Work space for integer values
        int* liw,       # The length of the integer workspace
        double* rpar,   # Double-precision parameters to pass to the residual and Jacobian functions
        int* ipar,      # Integer parameters to pass to the residual and Jacobian functions
        void* jac       # The Jacobian function
    )

################################################################################

class DASSLError(Exception):
    """
    An exception class for exceptions relating to use of DASSL.
    """
    
    def __init__(self, msg):
        self.msg = msg
    
    def __str__(self):
        return self.msg

################################################################################

cdef class DASSL:
    """
    A base class for using the DASSL differential algebraic system solver by
    L. R. Petzold. DASSL can be used to solve systems of the form
    
    .. math:: \\mathbf{g} \\left(t, \\mathbf{y}, \\frac{d \\mathbf{y}}{dt} \\right) = \\mathbf{0}
    
    where :math:`t` is the independent variable and :math:`\\mathbf{y}` is the
    vector of dependent variables. DASSL uses the backward differentiation
    formulas of orders one through five, and so is suitable for stiff problems.
    In particular, DASSL is much more robust than VODE, the solver used by
    SciPy.
    
    The DASSL solver options can be set by modifying the following attributes:
    
    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `atol`              ``object``          The absolute tolerance, either a scalar or a vector
    `rtol`              ``object``          The relative tolerance, either a scalar or a vector
    `maxOrder`          ``double``          The maximum order of the backward differentiation formulas, from ``1`` to ``5`` (default is ``5``)
    `initialStep`       ``double``          An initial step size to use, or ``0`` to let DASSL choose
    `maximumStep`       ``double``          An maximum allowed step size to use, or ``0``
    `tstop`             ``object``          A value of the independent variable to avoid integrating past, e.g. if system is undefined there
    `bandwidths`        ``object``          If the Jacobian matrix is banded, contains the lower and upper half-bandwidths
    `nonnegative`       ``bool``            :data:`True` to force nonnegativity constraint on dependent variable, :data:`False` if not
    =================== =================== ====================================
    
    These options are passed to DASSL when the :meth:`initialize()` method is
    called.
    """
    

    
    def __init__(self, maxOrder=5, initialStep=0, maximumStep=0, tstop=None, bandwidths=None, nonnegative=False, **kwargs):
        self.maxOrder = maxOrder
        self.initialStep = initialStep
        self.maximumStep = maximumStep
        self.tstop = tstop
        self.bandwidths = bandwidths
        self.nonnegative = nonnegative
    
    cpdef initialize(self, double t0, np.ndarray y0, np.ndarray dydt0=None, np.ndarray senpar=np.zeros(1, np.int32), atol=1e-16, rtol=1e-8):
        """
        Initialize the DASSL solver by setting the initial values of the
        independent variable `t0`, dependent variables `y0`, and first
        derivatives `dydt0`. If provided, the derivatives must be consistent
        with the other initial conditions; if not provided, DASSL will attempt
        to estimate a consistent set of initial values for the derivatives.
        You can also set the absolute and relative tolerances `atol` and `rtol`,
        respectively, either as single values for all dependent variables or 
        individual values for each dependent variable. Here, `senpar` is a dummy 
        variable and is unused since DASSL does not have sensitivity analysis support.
        """
        
        cdef double rwork[3]
        cdef int iwork[3]
        cdef int i
        cdef int neq, lrw, liw
        
        # Determine the number of equations
        neq = len(y0)
        if dydt0 is not None and len(dydt0) != neq:
            raise DASSLError('Expected %i values of dydt0, got %i.' % (neq, len(dydt0)))
        
        # Initialize all DASSL options to default values (i.e. all zeros)
        # Note that only the first 11 elements are used in DASSL
        self.info = np.zeros(15, np.int32)
        
        # Tell DASSL we are initializing the problem
        self.info[0] = 0
        
        # Tell DASSL about format of absolute and relative tolerances
        try:
            atol = float(atol)
            rtol = float(rtol)
            # If we're here then atol and rtol are both scalars
            self.info[1] = 0
        except TypeError:
            # If we're here then atol and rtol should both be arrays
            self.info[1] = 1
        self.atol = np.array(atol, np.float64)
        self.rtol = np.array(rtol, np.float64)
            
        # Set this option in advance() or step()
        self.info[2] = 0
        
        # Set tstop if specified
        if self.tstop:
            self.info[3] = 1
            rwork[0] = float(self.tstop)
        else:
            self.info[3] = 0
            rwork[0] = 0
            
        # Set whether or not jacobian function is provided
        # This is determined by whether or not you have implemented a
        # jacobian method
        if hasattr(self, 'jacobian'):
            self.info[4] = 1
        else:
            self.info[4] = 0
            
        # For banded matrices, set half-bandwidths
        if self.bandwidths:
            self.info[5] = 1
            iwork[0] = int(self.bandwidths[0])
            iwork[1] = int(self.bandwidths[1])
        else:
            self.info[5] = 0
            iwork[0] = 0
            iwork[1] = 0
        
        # Set maximum step size
        if self.maximumStep:
            self.info[6] = 1
            rwork[1] = float(self.maximumStep)
        else:
            self.info[6] = 0
            rwork[1] = 0
        
        # Set initial step size
        if self.initialStep:
            self.info[7] = 1
            rwork[2] = float(self.initialStep)
        else:
            self.info[7] = 0
            rwork[2] = 0
        
        # Set maximum order of backward differentiation formulas
        if self.maxOrder == 5 or not self.maxOrder:
            self.maxOrder = 5
            self.info[8] = 0
            iwork[2] = 0
        elif self.maxOrder >= 1 and self.maxOrder <= 5:
            self.info[8] = 1
            iwork[2] = int(self.maxOrder)
        else:
            raise DASSLError('The maxOrder attribute should have a value between 1 and 5.')
        
        # Turn on/off nonnegativity constraint
        if self.nonnegative:
            self.info[9] = 1
        else:
            self.info[9] = 0
                
        # If no initial derivatives provided, tell DASSL to attempt to guess them
        if dydt0 is not None:
            self.info[10] = 1
        else:
            self.info[10] = 0
        
        # Allocate rwork array
        if self.info[4] == 1 and self.info[5] == 1:
            lrw = 40 + (self.maxOrder + 4) * neq + (2*self.bandwidths[0]+self.bandwidths[1]+1) * neq
        elif self.info[4] == 0 and self.info[5] == 1:
            lrw = 40 + (self.maxOrder + 4) * neq + (2*self.bandwidths[0]+self.bandwidths[1]+1) * neq + 2*(neq/(self.bandwidths[0]+self.bandwidths[1]+1)+1)
        else:
            lrw = 40 + (self.maxOrder + 4) * neq + neq * neq
        self.rwork = np.zeros(lrw, np.float64)
        for i in range(3): self.rwork[i] = rwork[i]
        
        # Declare iwork array
        liw = 20 + neq
        self.iwork = np.zeros(liw, np.int32)
        for i in range(3): self.iwork[i] = iwork[i]
        
        # We don't use the rpar or ipar arrays, so set them to dummy values
        self.rpar = np.zeros(1, np.float64)
        self.ipar = np.zeros(1, np.int32)
        
        # Set initial conditions
        self.t = t0
        self.y = np.zeros(neq, np.float64)
        for i in range(neq):
            self.y[i] = y0[i]
        self.dydt = np.zeros(neq, np.float64)
        if dydt0 is not None:
            for i in range(neq):
                self.dydt[i] = dydt0[i]
    
    cpdef advance(self, double tout):
        """
        Simulate from the current value of the independent variable to a 
        specified value `tout`, taking as many steps as necessary. The resulting
        values of :math:`t`, :math:`\\mathbf{y}`, and 
        :math:`\\frac{d \\mathbf{y}}{dt}` can then be accessed via the `t`, `y`,
        and `dydt` attributes.
        """
        
        cdef int idid
        
        # Tell DASSL to return solution at tout
        self.info[2] = 0 
        # Call DASSL
        idid = self.solve(tout)
        return idid
        
    cpdef step(self, double tout):
        """
        Perform one simulation step from the current value of the independent 
        variable toward (but not past) a specified value `tout`. The resulting
        values of :math:`t`, :math:`\\mathbf{y}`, and 
        :math:`\\frac{d \\mathbf{y}}{dt}` can then be accessed via the `t`, `y`,
        and `dydt` attributes.
        """
        
        cdef int idid
        
        # Tell DASSL to only take one simulation step towards tout
        self.info[2] = 1
        # Call DASSL
        idid = self.solve(tout)
        return idid
        
    cdef solve(self, double tout):
        """
        Invoke DASSL with the given state of the object.
        """
        
        # Set the global DASSL object to this object (so we can get back to
        # this object's residual and jacobian methods
        global dasslObject
        dasslObject = self
        
        cdef int neq = self.y.shape[0]
        cdef int lrw = self.rwork.shape[0]
        cdef int liw = self.iwork.shape[0]
        cdef bint first = True
        
        cdef void* res = <void*> residual
        cdef void* jac = <void*> 0
        if self.info[4] == 1:
            jac = <void*> jacobian
        
        # Call DASSL
        while first or self.idid == -1:
            ddassl_(
                res,
                &(neq),
                &(self.t),
                <np.float64_t*> self.y.data,
                <np.float64_t*> self.dydt.data,
                &(tout),
                <int*> self.info.data,
                <np.float64_t*> self.rtol.data,
                <np.float64_t*> self.atol.data,
                &(self.idid),
                <np.float64_t*> self.rwork.data,
                &(lrw),
                <int*> self.iwork.data,
                &(liw),
                <np.float64_t*> self.rpar.data,
                <int*> self.ipar.data,
                jac
            )
            first = False
            if self.idid == -1:
                print('Attempting another 500 steps...')
                self.info[0] = 1
        
        # DASSL wrote onto the self.idid parameter automatically
        # Let's return it to the user
        return self.idid
        
    @cython.boundscheck(False)
    def residual(self, double t, np.ndarray[np.float64_t,ndim=1] y, np.ndarray[np.float64_t,ndim=1] dydt, **kwargs):
        """
        Evaluate the residual function for this model, given the current value
        of the independent variable `t`, dependent variables `y`, and first
        derivatives `dydt`. Return a numpy array with the values of the residual
        function and an integer with status information (0 if okay, -2 to
        terminate).
        """
        print('DASSLError: You must implement the residual() method in your derived class.')
        return np.zeros(y.shape[0], np.float64), -2
    
    #@cython.boundscheck(False)
    #def jacobian(self, double t, np.ndarray[np.float64_t,ndim=1] y, np.ndarray[np.float64_t,ndim=1] dydt, double cj):
    #   """
    #   Evaluate the Jacobian matrix for this model, given the current value
    #   of the independent variable `t`, dependent variables `y`, and first
    #   derivatives `dydt`. Return a numpy array with the values of the 
    #   Jacobian matrix.
    #   """
    #   print('DASSLError: You must implement the jacobian() method in your derived class.')
    #   return np.zeros((y.shape[0],y.shape[0]), np.float64)
        
################################################################################

# A module-level variable that contains the currently active DASSL object
# The residual and jacobian functions use this to call the object's residual
# and jacobian functions
cdef DASSL dasslObject

@cython.boundscheck(False)
cdef void residual(double* t, double* y, double* yprime, double* delta, int* ires, double* rpar, int* ipar):
    cdef np.ndarray[np.float64_t,ndim=1] res
    cdef int i
    res, ires[0] = dasslObject.residual(dasslObject.t, dasslObject.y, dasslObject.dydt)
    for i in range(res.shape[0]):
        delta[i] = res[i]

@cython.boundscheck(False)
cdef void jacobian(double* t, double* y, double* yprime, double* pd, double* cj, double* rpar, int* ipar):
    cdef np.ndarray[np.float64_t,ndim=2] jac
    cdef int i, j
    jac = dasslObject.jacobian(dasslObject.t, dasslObject.y, dasslObject.dydt, cj[0])
    N = jac.shape[0]
    for i in range(N):
        for j in range(N):
            pd[j*N+i] = jac[i,j]

