#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   PyDAS - A Python wrapper for several differential algebraic system solvers
#
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
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

import numpy

import pydas._dassl as _dassl

################################################################################

class DASSLError(Exception):
    """
    An exception class for exceptions relating to use of DASSL. Pass a string
    describing the circumstances of the exceptional behavior.
    """
    pass

################################################################################

class DASSL:
    """
    A Pythonic interface to the DASSL differential algebraic equation solver.
    """

    def __init__(self, maxOrder=5, initialStep=0, maximumStep=0, tstop=None, bandwidths=None, nonnegative=False):
        self.maxOrder = maxOrder
        self.initialStep = initialStep
        self.maximumStep = maximumStep
        self.tstop = tstop
        self.bandwidths = bandwidths
        self.nonnegative = nonnegative
    
    def initialize(self, t0, y0, dydt0=None, atol=1e-16, rtol=1e-8):
        """
        Initialize the DASSL solver by setting the initial values of the
        independent variable `t0`, dependent variables `y0`, and first
        derivatives `dydt0`. If provided, the derivatives must be consistent
        with the other initial conditions; if not provided, DASSL will attempt
        to estimate a consistent set of initial values for the derivatives.
        You can also set the absolute and relative tolerances `atol` and `rtol`,
        respectively, either as single values for all dependent variables or 
        individual values for each dependent variable.
        """
        
        rwork = numpy.zeros(3, numpy.float64)
        iwork = numpy.zeros(3, numpy.int32)
        
        # Determine the number of equations
        neq = len(y0)
        if dydt0 is not None and len(dydt0) != neq:
            raise DASSLError('Expected %i values of dydt0, got %i.' % (neq, len(dydt0)))
        
        # Initialize all DASSL options to default values (i.e. all zeros)
        # Note that only the first 11 elements are used in DASSL
        self.info = numpy.zeros(15, numpy.int32)
        
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
        self.atol = numpy.array(atol, numpy.float64)
        self.rtol = numpy.array(rtol, numpy.float64)
            
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
        self.info[4] = 1 if hasattr(self, 'jacobian') else 0
            
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
            lrw = 40 + (self.maxOrder + 4) * neq + (2*self.bandwidths[0]+self.bandwidths[0]+1) * neq
        elif self.info[4] == 0 and self.info[5] == 1:
            lrw = 40 + (self.maxOrder + 4) * neq + (2*self.bandwidths[0]+self.bandwidths[0]+1) * neq + 2*(neq/(2*self.bandwidths[0]+self.bandwidths[0]+1)+1)
        else:
            lrw = 40 + (self.maxOrder + 4) * neq + neq * neq
        self.rwork = numpy.zeros(lrw, numpy.float64)
        for i in range(3): self.rwork[i] = rwork[i]
        
        # Declare iwork array
        liw = 20 + neq
        self.iwork = numpy.zeros(liw, numpy.int32)
        for i in range(3): self.iwork[i] = iwork[i]
        
        # We don't use the rpar or ipar arrays, so set them to dummy values
        self.rpar = numpy.zeros(1, numpy.float64)
        self.ipar = numpy.zeros(1, numpy.int32)
        # However, we do use them to get neq to the res and jac functions
        self.ipar[0] = neq
        
        # Set initial conditions
        self.t = t0
        self.y = numpy.zeros(neq, numpy.float64)
        for i in range(neq):
            self.y[i] = y0[i]
        self.dydt = numpy.zeros(neq, numpy.float64)
        if dydt0 is not None:
            for i in range(neq):
                self.dydt[i] = dydt0[i]
        
        self.idid = 0
        
    def advance(self, tout):
        """
        Simulate from the current value of the independent variable to a 
        specified value `tout`, taking as many steps as necessary. The resulting
        values of :math:`t`, :math:`\\mathbf{y}`, and 
        :math:`\\frac{d \\mathbf{y}}{dt}` can then be accessed via the `t`, `y`,
        and `dydt` attributes.
        """
        # Tell DASSL to return solution at tout
        self.info[2] = 0 
        # Call DASSL
        return self.solve(tout)
        
    def step(self, tout):
        """
        Perform one simulation step from the current value of the independent 
        variable toward (but not past) a specified value `tout`. The resulting
        values of :math:`t`, :math:`\\mathbf{y}`, and 
        :math:`\\frac{d \\mathbf{y}}{dt}` can then be accessed via the `t`, `y`,
        and `dydt` attributes.
        """
        # Tell DASSL to only take one simulation step towards tout
        self.info[2] = 1
        # Call DASSL
        return self.solve(tout)
        
    def solve(self, tout):
        """
        Invoke DASSL with the given state of the object.
        """
        res = self.__residual
        jac = self.__jacobian if hasattr(self, 'jacobian') else dummy_jacobian
        self.t = numpy.array(self.t)
        
        first = True
        while first or self.idid == -1:
            self.idid = _dassl.ddassl(
                res,                # The residual function
                self.t,             # The current value of the independent variable
                self.y,             # The current values of the dependent variables
                self.dydt,          # The current values of the first derivatives of the independent variable
                tout,               # The value of the independent variable we are integrating toward
                self.info,          # Configuration information for DASSL
                self.rtol,          # The relative tolerance
                self.atol,          # The absolute tolerance
                self.rwork,         # Allocated double-precision array for DASSL to use for storing floating-point values
                self.iwork,         # Allocated integer array for DASSL to use for storing integer values
                self.rpar,
                self.ipar,
                jac,
            )
            
            first = False
            if self.idid == -1:
                print('Attempting another 500 steps...')
                self.info[0] = 1
        
        # DASSL wrote onto the self.idid parameter automatically
        # Let's return it to the user
        return self.idid
    
    def __residual(self, t, y, dydt):
        """
        Evaluate the residual function for this model, given the current value
        of the independent variable `t`, dependent variables `y`, and first
        derivatives `dydt`. Return a numpy array with the values of the residual
        function and an integer with status information (0 if okay, -2 to
        terminate).
        """
        try:
            return self.residual(t, y, dydt)
        except AttributeError:
            raise DASSLError('You must implement the residual() method in your derived class.')
    
    def __jacobian(self, t, y, dydt, cj):
       """
       Evaluate the Jacobian matrix for this model, given the current value
       of the independent variable `t`, dependent variables `y`, and first
       derivatives `dydt`. Return a numpy array with the values of the 
       Jacobian matrix.
       """
       try:
           return self.jacobian(t, y, dydt, cj)
       except AttributeError:
           raise DASSLError('You must implement the jacobian() method in your derived class.')

def dummy_jacobian(t, y, dydt, cj):
    return numpy.zeros(y.shape[0], numpy.float64)
