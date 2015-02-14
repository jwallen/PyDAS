#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This file contains unit tests for PyDASPK.
"""

import unittest

from pydas.daspk import DASPK
import math
import numpy
import matplotlib.pyplot as plt
################################################################################

class DASPKSimpleModel(DASPK):
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
    
    def residual(self, t, y, dydt, senpar):
        delta = numpy.zeros(y.shape[0], numpy.float64)
        delta[0] = -self.k1 * y[0] - dydt[0]
        delta[1] =  self.k1 * y[0] - self.k2 * y[1] - dydt[1]
        delta[2] =  self.k2 * y[1] - dydt[2]
        return delta, 0
    
    def jacobian(self, t, y, dydt, cj, senpar):
        pd = -cj * numpy.identity(y.shape[0], numpy.float64)
        pd[0,0] += -self.k1
        pd[1,0] += self.k1
        pd[1,1] += -self.k2
        pd[2,1] += self.k2
        return pd

################################################################################

class DASPKSensitivityModel(DASPK):
    """
    A model of first-order irreversible reactions in series
    
        A -> B -> C
    
    occurring in a batch reactor. In such a system the concentration
    of the intermediate B has a maximum that depends on the relative
    rate constants `k1` and `k2`. These are stored as data members of
    the class so that they are available to the residual function.
    """
    
    def __init__(self, sensitivity, sensmethod):
        self.sensitivity = sensitivity
        self.sensmethod = sensmethod
    
    def residual(self, t, y, dydt, senpar):
        delta = numpy.zeros(9, numpy.float64)
        delta[0] = -senpar[0]* y[0] 
        delta[1] =  senpar[0]* y[0] - senpar[1] * y[1]
        delta[2] =  senpar[1] * y[1]
        numReactions = 2
        numState = 3
        for j in range(numReactions):
            for i in range(numState):
                for k in range(numState):
                    delta[(j+1)*numState + i] += self.jacobian(t,y,dydt,0,senpar)[i,k]*y[(j+1)*numState + k]                     
                delta[(j+1)*numState + i] += self.dgdk(y)[i,j]
                
                
        for i in range(9):
            delta[i] -= dydt[i]
        return delta, 1
    
    def jacobian(self, t, y, dydt, cj, senpar):
        pd = -cj * numpy.identity(3, numpy.float64)
        pd[0,0] += -senpar[0]
        pd[1,0] += senpar[0]
        pd[1,1] += -senpar[1]
        pd[2,1] += senpar[1]
        return pd
    
    def dgdk(self,y):
        deriv = numpy.zeros((3,2),numpy.float64)
        deriv[0,0] = -y[0]
        deriv[1,0] = y[0]
        deriv[1,1] = -y[1]
        deriv[2,1] = y[1]
        return deriv

################################################################################

class DASPKCheck(unittest.TestCase):
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
        model = DASPKSimpleModel(k1, k2)
        t0 = 0.0; y0 = numpy.array([1.0, 0.0, 0.0], numpy.float64)
        dydt0 = model.residual(t0, y0, numpy.zeros(3, numpy.float64),numpy.zeros(3, numpy.float64))[0]
        senpar = numpy.array([0.0,0.0,0.0], numpy.float64)        
        
        model.initialize(t0, y0, dydt0, senpar, atol=1e-16, rtol=1e-8)
        
        tmax = 100; iter = 0; maxiter = 1000
        tvec = []
        Avec = []
        Bvec = []
        Cvec = []
        while iter < 1000 and model.t < 16:
            model.step(tmax)
            t = model.t
            tvec.append(t)
            A, B, C = model.y
            Avec.append(A)
            Bvec.append(B)
            Cvec.append(C)
            Atrue = math.exp(-k1*t)
            Btrue = k1 / (k2 - k1) * (math.exp(-k1*t) - math.exp(-k2*t))
            Ctrue = 1.0 - Atrue - Btrue
            if Atrue > 1e-8:
                self.assertAlmostEqual(A / Atrue, 1.0, 6)
            if Btrue > 1e-8:
                self.assertAlmostEqual(B / Btrue, 1.0, 6, 'At t = %g: B = %g, but Btrue = %g' % (t, B, Btrue))
            if Ctrue > 1e-8:
                self.assertAlmostEqual(C / Ctrue, 1.0, 6, 'At t = %g: C = %g, but Ctrue = %g' % (t, C, Ctrue))
        
        plt.plot(tvec,Avec,tvec,Bvec,tvec,Cvec)
        plt.show()

    def testSensitivity(self):
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

        This tests whether finite difference sensitivities work in daspik
        """

        k1 = 1.0; k2 = 0.25
        model = DASPKSensitivityModel(True, 1)
        t0 = 0.0
        n0 = numpy.array([1.0, 0.0, 0.0], numpy.float64)   #state variables 
        senpar = numpy.array([k1,k2], numpy.float64)        
        neq = len(n0)*(len(senpar)+1)
        y0 = numpy.zeros(neq, numpy.float64)
        atol = numpy.ones(neq,numpy.float64)*1e-10
        rtol = numpy.ones(neq,numpy.float64)*1e-8
        for i in range(len(n0)):
            y0[i] = n0[i]
            atol[i] = 1e-16
            rtol[i] = 1e-8
        dydt0 = numpy.zeros(neq, numpy.float64)
        res = model.residual(t0, y0, numpy.zeros(neq, numpy.float64),senpar)[0]
        dydt0[:len(res)] = res
        model.initialize(t0, y0, dydt0, senpar, atol=atol, rtol=rtol)

        tmax = 100; iter = 0; maxiter = 1000
        tvec = []
        dAdk1vec = []
        dBdk1vec = []
        dCdk1vec = []        
        dAdk2vec = []
        dBdk2vec = []
        dCdk2vec = []
        
        Avec = []
        Bvec = []
        Cvec = []
        while iter < 1000 and model.t < 16:
            model.step(tmax)
            t = model.t
            A, B, C = model.y[:3]
            tvec.append(model.t)
            dAdk1vec.append(model.y[3]*k1/A)
            dBdk1vec.append(model.y[4]*k1/B)
            dCdk1vec.append(model.y[5]*k1/C)
            dAdk2vec.append(model.y[6]*k2/A)
            dBdk2vec.append(model.y[7]*k2/B)
            dCdk2vec.append(model.y[8]*k2/C)
            
            Avec.append(A)
            Bvec.append(B)
            Cvec.append(C)
            Atrue = math.exp(-k1*t)
            Btrue = k1 / (k2 - k1) * (math.exp(-k1*t) - math.exp(-k2*t))

            Ctrue = 1.0 - Atrue - Btrue
            if Atrue > 1e-8:
                self.assertAlmostEqual(A / Atrue, 1.0, 6)
            if Btrue > 1e-8:
                self.assertAlmostEqual(B / Btrue, 1.0, 6, 'At t = %g: B = %g, but Btrue = %g' % (t, B, Btrue))
            if Ctrue > 1e-8:
                self.assertAlmostEqual(C / Ctrue, 1.0, 6, 'At t = %g: C = %g, but Ctrue = %g' % (t, C, Ctrue))
                
        print "At 16 seconds after simulation..."
        print "dA/dk1 = %f" % dAdk1vec[-1]        
        print "dB/dk1 = %f" % dBdk1vec[-1]        
        print "dC/dk1 = %f" % dCdk1vec[-1]            
        print "dA/dk2 = %f" % dAdk2vec[-1]
        print "dB/dk2 = %f" % dBdk2vec[-1]
        print "dC/dk2 = %f" % dCdk2vec[-1]
        plt.figure(1)
        plt.plot(tvec, Avec, label='A')
        plt.plot(tvec, Bvec, label='B')
        plt.plot(tvec, Cvec, label='C')
        plt.figure(2)
        plt.plot(tvec,dAdk1vec,label='dA/dk1')
        plt.plot(tvec,dAdk2vec,label='dA/dk2')
        plt.plot(tvec,dBdk1vec,label='dB/dk1')
        plt.plot(tvec,dBdk2vec,label='dB/dk2')
        plt.plot(tvec,dCdk1vec,label='dC/dk1')
        plt.plot(tvec,dCdk2vec,label='dC/dk2')
        
        plt.legend()
        plt.show()

################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
