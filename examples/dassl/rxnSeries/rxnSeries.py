#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This file uses DASSL to solve a model of a set of irreversible first-order 
reactions in series occurring in a well-mixed batch reactor:

    A -> B -> C

The governing equations are

    d[A]/dt = -k1*[A]
    d[B]/dt = k1*[A] - k2*[B]
    d[C]/dt = k2*[B]

with initial conditions

    [A](0) = 1.0    [B](0) = 0.0    [C] = 0.0

and k1 = 1.0 and k2 = 0.25.
"""

import sys
sys.path.append('.')

from pydas import DASSL
import numpy

class Model(DASSL):
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
	
	#def jacobian(self, t, y, dydt, cj):
		#pd = -cj * numpy.identity(y.shape[0], numpy.float64)
		#pd[0,0] += -self.k1
		#pd[1,0] += self.k1
		#pd[1,1] += -self.k2
		#pd[2,1] += self.k2
		#return pd

################################################################################

if __name__ == '__main__':
	
	# Initialize solution vectors
	t = []
	y = []
	
	# Set maximum simulation time and maximum number of simulation steps to allow
	tmax = 16
	maxiter = 1000
	
	# Initialize the model
	model = Model(k1=1.0, k2=0.25)
	t0 = 0.0; y0 = numpy.array([1.0, 0.0, 0.0], numpy.float64)
	# Since the model is explicit -- dy/dt = f(t, y) -- it is easy to
	# evaluate dydt0; doing this is recommended because it lets DASSL
	# choose a proper initial time step
	dydt0 = - model.residual(t0, y0, numpy.zeros(3, numpy.float64))[0]
	model.initialize(t0, y0, dydt0)
	
	# Generate the solution by stepping until tmax is reached
	# This will give you the solution at time points automatically selected
	# by the solver
	iter = 0
	while iter < maxiter and model.t < tmax:
		model.step(tmax)
		t.append(model.t)
		# You must make a copy of y because it is overwritten by DASSL at
		# each call to step()
		y.append(model.y.copy())
	
	# Convert the solution vectors to numpy arrays
	t = numpy.array(t, numpy.float64)
	y = numpy.array(y, numpy.float64)
	
	# If matplotlib is installed, show a plot of the results
	# Otherwise do nothing
	try:
		import pylab
		pylab.plot(t, y)
		pylab.legend(['A', 'B', 'C'], loc=1)
		pylab.xlabel('Time')
		pylab.ylabel('Concentration')
		pylab.show()
	except ImportError:
		pass
	
	