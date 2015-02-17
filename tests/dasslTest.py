#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This file contains unit tests for PyDAS.
"""

import unittest

from pydas.dassl import DASSL
import math
import numpy

################################################################################


class SimpleModel(DASSL):
    """
    A model of first-order irreversible reactions in series
    
        A -> B -> C
    
    occurring in a batch reactor. In such a system the concentration
    of the intermediate B has a maximum that depends on the relative
    rate constants `k1` and `k2`. These are stored as data members of
    the class so that they are available to the residual function.
    """
    
    def __init__(self, k1=0.0, k2=0.0):
        self.k1 = k1
        self.k2 = k2
    
    def residual(self, t, y, dydt):
        delta = numpy.zeros(y.shape[0], numpy.float64)
        delta[0] = -self.k1 * y[0] - dydt[0]
        delta[1] =  self.k1 * y[0] - self.k2 * y[1] - dydt[1]
        delta[2] =  self.k2 * y[1] - dydt[2]
        return delta, 0
    
    def jacobian(self, t, y, dydt, cj):
        pd = -cj * numpy.identity(y.shape[0], numpy.float64)
        pd[0,0] += -self.k1
        pd[1,0] += self.k1
        pd[1,1] += -self.k2
        pd[2,1] += self.k2
        return pd

#################################################################################

class DASSLCheck(unittest.TestCase):
    """
    Contains unit tests of the DASSL wrapper.
    """
    
    def testSimple(self):
        """
        This test solves a simple set of ODEs numerically using DASSL and
        compares the result to a known analytical solution. The physical
        system this problem represents is two first-order chemical reactions
        in series in a homogeneous batch reactor.
        
        The governing equations are

            dA/dt = -k1*A
            dB/dt = k1*A - k2*B
            dC/dt = k2*B

        with initial conditions

            A(0) = 1.0      B(0) = 0.0      C(0) = 0.0

        and k1 = 1.0 and k2 = 0.25. The analytical solution is
        
            A(t) = exp(-k1*t)
            B(t) = k1/(k2-k1) * (exp(-k1*t) - exp(-k2*t))
            C(t) = (1 - k2/(k2-k1) * exp(-k1*t) - k1/(k2-k1) * exp(-k2*t)
            
        Note that C(t) = A(0) - A(t) - B(t) since the governing equations are
        conservative.
        """
        
        k1 = 1.0; k2 = 0.25
        model = SimpleModel(k1, k2)
        t0 = 0.0; y0 = numpy.array([1.0, 0.0, 0.0], numpy.float64)
        dydt0 = model.residual(t0, y0, numpy.zeros(3, numpy.float64))[0]
        model.initialize(t0, y0, dydt0, atol=1e-16, rtol=1e-8)
        
        tmax = 100; iter = 0; maxiter = 1000
        while iter < 1000 and model.t < 16:
            model.step(tmax)
            t = model.t
            A, B, C = model.y
            
            Atrue = math.exp(-k1*t)
            Btrue = k1 / (k2 - k1) * (math.exp(-k1*t) - math.exp(-k2*t))
            Ctrue = 1.0 - Atrue - Btrue
            if Atrue > 1e-8:
                self.assertAlmostEqual(A / Atrue, 1.0, 6)
            if Btrue > 1e-8:
                self.assertAlmostEqual(B / Btrue, 1.0, 6, 'At t = %g: B = %g, but Btrue = %g' % (t, B, Btrue))
            if Ctrue > 1e-8:
                self.assertAlmostEqual(C / Ctrue, 1.0, 6, 'At t = %g: C = %g, but Ctrue = %g' % (t, C, Ctrue))
            
################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
