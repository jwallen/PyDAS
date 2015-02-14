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
This `Cython <http://www.cython.org/>`_ module exposes the DASPK3.1 differential 
algebraic system solver to Python and provides a Python extension type, the 
:class:`DASPK` base class, for working with DASPK.

To use DASPK, write a Python class or Cython extension type that derives from
the :class:`DASPK` class and implement the :meth:`residual`
method, which accepts :math:`t`, :math:`\\mathbf{y}`, and 
:math:`\\mathbf{\\frac{d \\mathbf{y}}{dt}}` as arguments and returns the
corresponding value of :math:`\\mathbf{g} \\left(t, \\mathbf{y}, \\frac{d \\mathbf{y}}{dt} \\right), \\mathbf{senpar}`.
Run by calling the :meth:`initialize` method to set the initial conditions
and solver tolerances, then by using the :meth:`advance` or :meth:`step`
methods to move forward or backward in the independent variable. If you 
know the form of the Jacobian, you can implement the :meth:`jacobian`
method to provide it to the solver.

You can implement your derived class in pure Python, but for a significant
speed boost consider using Cython to compile the module. You can see the
proper Cython syntax for the residual and jacobian methods by looking at the
corresponding methods in the :class:`DASPK` base class.
"""

import numpy as np
cimport numpy as np

cimport cython

################################################################################

# Expose the (double-precision) DASPK function
cdef extern from "daspk.h":
    int ddaspk_(
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
        void* jac,       # The Jacobian function
        void* psol,      # Psol function
        double* senpar,  # Sensitivity parameters
        void* g_res     # G_res function
    )

################################################################################

class DASPKError(Exception):
    """
    An exception class for exceptions relating to use of DASSL.
    """
    
    def __init__(self, msg):
        self.msg = msg
    
    def __str__(self):
        return self.msg

################################################################################
cdef class DASPK:
    """
    A base class for using the DASPK differential algebraic system solver by
    L. R. Petzold. DASPK can be used to solve systems of the form
    
    .. math:: \\mathbf{g} \\left(t, \\mathbf{y}, \\frac{d \\mathbf{y}}{dt} \\right) = \\mathbf{0}
    
    where :math:`t` is the independent variable and :math:`\\mathbf{y}` is the
    vector of dependent variables. DASPK uses the backward differentiation
    formulas of orders one through five, and so is suitable for stiff problems.
    In particular, DASPK is much more robust than VODE, the solver used by
    SciPy.
    
    The DASPK solver options can be set by modifying the following attributes:
    
    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `atol`              ``object``          The absolute tolerance, either a scalar or a vector
    `rtol`              ``object``          The relative tolerance, either a scalar or a vector
    `maxOrder`          ``int``             The maximum order of the backward differentiation formulas, from ``1`` to ``5`` (default is ``5``)
    `initialStep`       ``double``          An initial step size to use, or ``0`` to let DASSL choose
    `maximumStep`       ``double``          An maximum allowed step size to use, or ``0``
    `tstop`             ``object``          A value of the independent variable to avoid integrating past, e.g. if system is undefined there
    `bandwidths`        ``object``          If the Jacobian matrix is banded, contains the lower and upper half-bandwidths
    `nonnegative`       ``bool``            :data:`True` to force nonnegativity constraint on dependent variable, :data:`False` if not
    `sensitivity`       ``bool``            :data:`True` if sensitivities are to be calculated, :data:`False` if not
    `sensmethod`        ``int``             ``0`` to use default simultaneous corrector, ``1`` for staggered corrector, ``2`` for staggered direct method
    =================== =================== ====================================
    
    These options are passed to DASSL when the :meth:`initialize()` method is
    called.
    """
    
    def __init__(self, maxOrder=5, initialStep=0, maximumStep=0, tstop=None, bandwidths=None, nonnegative=False, sensitivity=False, sensmethod=0):
        self.maxOrder = maxOrder
        self.initialStep = initialStep
        self.maximumStep = maximumStep
        self.tstop = tstop
        self.bandwidths = bandwidths
        self.nonnegative = nonnegative
        self.sensitivity = sensitivity
        self.sensmethod = sensmethod
    
    cpdef initialize(self, double t0, np.ndarray y0, np.ndarray dydt0=None, np.ndarray senpar=None, atol=1e-16, rtol=1e-8):
        """
        Initialize the DASPK solver by setting the initial values of the
        independent variable `t0`, dependent variables `y0`, and first
        derivatives `dydt0`. If provided, the derivatives must be consistent
        with the other initial conditions; if not provided, DASPK will attempt
        to estimate a consistent set of initial values for the derivatives.
        You can also set the absolute and relative tolerances `atol` and `rtol`,
        respectively, either as single values for all dependent variables or 
        individual values for each dependent variable.
        """
        
        cdef double rwork[3]
        cdef int iwork[3]
        cdef int i
        cdef int neq, lrw, liw, npar, ny, itmp1, itmp2, itmp3
        
        # Determine the number of equations
        neq = len(y0)
        if dydt0 is not None and len(dydt0) != neq:
            raise DASPKError('Expected %i values of dydt0, got %i.' % (neq, len(dydt0)))
        
        # Initialize all DASSL options to default values (i.e. all zeros)
        # Note that only the first 11 elements are used in DASSL
        self.info = np.zeros(30, np.int32)
        
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
        # Set this for direct methods only. Setting bandwidths improves performance of solver
        if self.bandwidths:
            self.info[5] = 1
            iwork[0] = int(self.bandwidths[0])  # lower bandwidth 
            iwork[1] = int(self.bandwidths[1])  # upper bandwidth
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
            raise DASPKError('The maxOrder attribute should have a value between 1 and 5.')
        
        # Turn on/off nonnegativity constraint
        # Invoking nonnegativity should only be done in extreme circumstances
        if self.nonnegative:
            self.info[9] = 1
            # nonnegative is a dummy variable, 
            # iwork must be further specified for inequality constraints if this is specified
        else:
            self.info[9] = 0
                
        # Are the initial T, Y, Yprime consisten?
        self.info[10] = 0    # yes they are
        # self.info[10] = 1   # no they are not consistent
        
        # Use direct method in DDASPK
        # self.info[11] = 0
        # For iterative Krylov method
        # self.info[11] = 1  
        
        # Related to Krylov method
        # self.info[12] = 0
        
        # Proceed with integration after initial condition calculation is done, used with info(11) > 0
        # self.info[13] = 0 # yes
        # self.info[13] = 1 # no
        
        # Related to Krylov method
        # self.info[14] = 0
        
        # Option to exclude algebraic variables from error test
        # self.info[15]
    
        # Used when initial condition not specified info(11) > 0
        # self.info[16]
        
        # Error control option for quadrature/output variables. This option is only used when info(28) > 0
        # self.info[17]
        
        # Sensitivity analysis toggling
        npar = 0  # number of parameter is zero unless sensitivity analysis is involved
        if self.sensitivity:
            # number of parameters involved in system to be solved, 
            # including parameters that appear only in the initial conditions
            self.senpar = senpar
            npar = len(senpar)
            self.info[18] = npar
            
            # Method used to calculate sensitivities
            # Finite difference method is used for sensitivity analysis     
            # Second order centered FDM is used
            #self.info[19] = 0   
            # self.info[19] = 1 # First order forward FDM is used
            self.info[19] = 2 # Provide a user-supplied RES, where RES should also compute the residuals of the sensitivity equations
            # IRES is used to determine whether to compute the sensitivities or not in the RES routine
            # self.info[19] = 3, 4, or 5 related to ADIFOR-generated routines
            
            # Perturbation factor used when finite differences are used for sensitivity
            self.info[20] = 0  # use default sensitivity
            # self.info[20] = 1 # alter the value for the perturbation used
                       
            # Number of parameters that are involved in RES
            self.info[21] = npar 
            
            # Error control for the sensitivity variables
            self.info[22] = 0 # include error testsl
            # self.info[22] = 1 # don't include error test
            
            # Sensitivity of a derived quantity
            # In addition to computing the sensitivity of the solution, you may want to compute the sensitivity
            # of a derived quantity with respect to the parameters
            self.info[23] = 0 # Do not include sensitivity for derived quantity
            # self.info[23] = 1 # Do sensitivity on derived quantity
            
            # Method option for sensitivity analysis
            self.info[24] = self.sensmethod
            
            # For sensitivity calculations using the adjoint method.  If not using, set both to 0.
            self.info[28] = 0
            self.info[29] = 0
        
        else:
            self.info[18] = 0
        
        self.info[25] = 0  # serial computing, change this if doing parallel computing
        self.info[26] = 0  # for parallel computing only if changed
        self.info[27] = 0  # option for efficient computation of quadrature variables                  
                           
        
        # Allocate rwork array
        ny = neq / (npar + 1) # total number of state variables
        lrw = 50 + max(self.maxOrder + 4, 7) * neq
#        if info[11] == 1: # Krylov method, not applicable here
#            lrw += (MAXL+3+MIN0(1,MAXL-KMP))*NY + (MAXL+3)*MAXL + 1 + LENWP
        
        # initialize itmp values
        itmp1 = 0
        itmp2 = 0
        itmp3 = 0
        
        if self.info[15] == 1:
            lrw += neq
        
        # Non-krylov method in use, ie. standard direct method
        if self.info[11] == 0:
            # dense matrix
            if self.info[5] == 0:
                lrw += ny * ny
                if self.info[4] == 3:
                    itmp1 = ny * (3 * ny)
            if self.info[5] == 1: # banded matrix
                lrw += (2*self.bandwidths[0]+self.bandwidths[1]+1) * ny
                if self.info[4] == 0:
                    itmp1 = 2*(ny/(self.bandwidths[0]+self.bandwidths[1]+1)+1)
                if self.info[4] == 3:
                    itmp1 = (self.bandwidths[0] + self.bandwidths[1] + 1) * (3 * ny)
        # Krylov method, we don't actually use this
        if self.info[11] == 1:
            if self.info[4] == 2:
                itmp1 = ny
        
        # index-2 two step process
        if self.info[10] == 4:
                if self.info[19] < 2 and self.info[18] == 0:
                    itmp2 = 2 * ny
                if self.info[19] < 2 and self.info[18] > 0:
                    itmp2 = 2* neq + 4 * ny + 2 * self.info[23]
                if self.info[19] > 2 and self.info[18] == 0:
                    itmp2 = ny
                if self.info[19] > 2 and self.info[18] > 0:
                    itmp2 = 2 * ny + self.info[21]
                    
        # sensitivity analysis is on
        if self.info[18] > 0:
            if self.info[19] < 2:
                itmp3 = 4 * ny + 2 * self.info[23]
            if self.info[19] == 3:
                itmp3 = self.info[18] * (2*ny + max(ny, self.info[23]) + self.info[21])
            if self.info[19] == 4:
                itmp3 = self.info[21]
#            if info(19) == 5:   # provide an adifor generated routine for g_res
#                itmp3 = info(18)*(ny+1) + rwork(17)*(2*ny*ny)
        
        lrw += max(itmp1, itmp2, itmp3)
            
        self.rwork = np.zeros(lrw, np.float64)
        for i in range(3): self.rwork[i] = rwork[i]
        
        # Declare iwork array
        liw = 40 + ny
        if self.info[9] == 1 or self.info[9] == 3:
            liw += ny
        if self.info[10] == 1 or self.info[10] == 3:
            liw += ny
        if self.info[15] == 1:
            liw += ny
        if self.info[10] == 4 or self.info[10] == 5:
            liw += 2*ny
        if self.info[11] == 0 and self.info[4] == 2:
            liw += ny
        if self.info[18] > 0 and self.info[19] == 4:
            liw += ny + 1 + self.rwork[16]*(2*ny*ny) + self.info[21]
        if self.info[11] == 0 and self.info[4] == 2:
            liw += 3*ny + self.info[21]
        if self.info[19] == 5:
            liw += 3*ny + self.info[21]
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
        
        # Tell DASPK to return solution at tout
        self.info[2] = 0 
        # Call DASPK
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
        Invoke DASPK with the given state of the object.
        """
        
        # Set the global DASSL object to this object (so we can get back to
        # this object's residual and jacobian methods
        global daspkObject
        daspkObject = self
        
        cdef int neq = self.y.shape[0]
        cdef int lrw = self.rwork.shape[0]
        cdef int liw = self.iwork.shape[0]
        cdef bint first = True
        
        cdef void* res = <void*> residual
        cdef void* jac = <void*> 0        
        cdef void* psol = <void*> 0
        cdef void* g_res = <void*> 0
        if self.info[4] == 1:
            jac = <void*> jacobian
        
        # Call DASPK
        while first or self.idid == -1:
            ddaspk_(
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
                jac,
                psol,
                <np.float64_t*> self.senpar.data,
                g_res,
            )
            first = False
            if self.idid == -1:
                print('Attempting another 500 steps...')
                self.info[0] = 1
        
        # DASPK wrote onto the self.idid parameter automatically
        # Let's return it to the user
        return self.idid
        
    @cython.boundscheck(False)
    def residual(self, double t, np.ndarray[np.float64_t,ndim=1] y, np.ndarray[np.float64_t,ndim=1] dydt, np.ndarray[np.float64_t,ndim=1] senpar):
        """
        Evaluate the residual function for this model, given the current value
        of the independent variable `t`, dependent variables `y`, and first
        derivatives `dydt`. Return a numpy array with the values of the residual
        function and an integer with status information (0 if okay, -2 to
        terminate).
        """
        print('DASPKError: You must implement the residual() method in your derived class.')
        return np.zeros(y.shape[0], np.float64), -2
    
    #@cython.boundscheck(False)
    #def jacobian(self, double t, np.ndarray[np.float64_t,ndim=1] y, np.ndarray[np.float64_t,ndim=1] dydt, double cj, np.ndarray[np.float64_t,ndim=1] senpar):
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
cdef DASPK daspkObject

@cython.boundscheck(False)
cdef void residual(double* t, double* y, double* yprime, double* cj, double* delta, int* ires, double* rpar, int* ipar, double *senpar):
    cdef np.ndarray[np.float64_t,ndim=1] res
    cdef int i
    res, ires[0] = daspkObject.residual(daspkObject.t, daspkObject.y, daspkObject.dydt, daspkObject.senpar)
    for i in range(res.shape[0]):
        delta[i] = res[i]

@cython.boundscheck(False)
cdef void jacobian(double* t, double* y, double* yprime, double* pd, double* cj, double* rpar, int* ipar, double* senpar, int* ijac):
    cdef np.ndarray[np.float64_t,ndim=2] jac
    cdef int i, j
    jac = daspkObject.jacobian(daspkObject.t, daspkObject.y, daspkObject.dydt, cj[0], daspkObject.senpar)
    N = jac.shape[0]
    for i in range(N):
        for j in range(N):
            pd[j*N+i] = jac[i,j]