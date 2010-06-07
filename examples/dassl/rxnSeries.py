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

import PyDAS.dassl as dassl
import numpy

class Model(dassl.DASSL):
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
	model.initialize(0, [1.0, 0.0, 0.0])
	# Take a small initial step, short enough that it is before the time
	# scale of interest
	model.step(1e-10)
	
	# Generate the solution by stepping until tmax is reached
	iter = 0
	while iter < maxiter and model.t < tmax:
		model.step(tmax)
		t.append(model.t)
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
	
	