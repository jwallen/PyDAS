#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This file investigates the methane oxidation kinetic mechanism found in

    A. A. Boni and R. C. Penner. "Sensitivity Analysis of a Mechanism for 
    Methane Oxidation Kinetics." Comb. Sci. Tech. 15, p. 99-106 (1977).

The mechanism can be integrated in time with various initial conditions with
or without sensitivity analysis.
"""

import numpy
from solver import *

if __name__ == '__main__':
	
	Nspec = 15
	Nrxn = 46
	
	sensitivity = False
	
	tmin = 1e-7; tmax = 1e-2
	maxiter = 1000
	
	# Set the rate coefficients at 2000 K
	# The first 23 are the forward rates, the second 23 the reverse rates
	# The data is taken from Table 1 of Boni and Penner
	# Units are cm^3/molecule*s (bimolecular) or cm^6/molecule^2*s (termolecular)
	K = numpy.array(
	    [    7.2e-17,      4.3e-11,      3.3e-11,      3.6e-12,     1.66e-10,
	        3.32e-14,     2.63e-11,     1.84e-10,      8.7e-12,     6.38e-16,
	        1.66e-10,     1.66e-10,     3.32e-10,      7.0e-14,      8.9e-13,
	         3.0e-11,      1.2e-11,     5.34e-12,      5.8e-33,     9.65e-32,
	        2.58e-10,     5.32e-33,     1.57e-11,
	        3.05e-31,     1.74e-13,     1.34e-12,      1.1e-13,      8.8e-17,
	         7.2e-20,      1.5e-14,     1.42e-14,      7.0e-15,      5.4e-32,
	         1.9e-19,      2.1e-20,      4.6e-19,      1.0e-35,      3.8e-13,
	         3.0e-12,      1.1e-11,      2.3e-11,     5.14e-21,     8.56e-20,
	         1.0e-15,      3.6e-14,      2.0e-12   ], 
	    numpy.float64
	)
	
	# Stoichiometry matrix (rows are species, columns are reactions)
	S = numpy.zeros((Nspec,Nrxn), numpy.int)
	S[:,0:23] = numpy.array([ #        |                   |                   |                   |
          [ -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ], # CH4
          [  1,  1,  1,  1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ], # CH3
          [  1,  0, -1,  0,  1,  0,  0,  0, -1,  1,  0,  0, -1,  1,  1,  1,  1, -1, -1, -1, -1, -1,  0 ], # H
          [  0, -1,  0,  1,  0,  1,  1, -1,  0,  0,  1, -1,  0,  0, -1, -1,  1,  1, -1, -1,  2,  0, -2 ], # OH
          [  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  1,  1,  0,  0,  1 ], # H2O
          [  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0 ], # H2
          [  0,  0,  0, -1, -1,  0, -1,  0,  0,  0, -1,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  1 ], # O
          [  0,  0,  0,  0,  1,  1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ], # CH2O
          [  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0, -1,  0 ], # O2
          [  0,  0,  0,  0,  0,  0,  1,  1,  1,  1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0 ], # CHO
          [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0 ], # CO
          [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0 ], # CO2
          [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  1,  0 ], # HO2
          [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ], # Ar
          [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ], # M
        ])
	S[:,Nrxn/2:Nrxn] = -S[:,0:Nrxn/2]
	
	# Reactant order matrix
	Sr = -S.copy()
	for i in range(Nspec):
		for j in range(Nrxn):
			if Sr[i,j] < 0: Sr[i,j] = 0
	Sr[14,0] = 1;   Sr[14,0+23] = 1
	Sr[14,9] = 1;   Sr[14,9+23] = 1
	Sr[14,13] = 1;  Sr[14,13+23] = 1
	Sr[13,18] = 1;  Sr[13,18+23] = 1
	Sr[4,19] = 1;   Sr[4,19+23] = 1
	Sr[14,21] = 1;  Sr[14,21+23] = 1
	
	# Set the desired initial conditions
	# Units are molecules/cm^3
	t0 = 0
	C0 = numpy.zeros(Nspec, numpy.float64)
	C0[0] = 2.265e18
	C0[8] = 4.530e18
	C0[13] = 1.790e19
	C0[5] = 1
	C0[11] = 1
	C0[14] = 2.469e19
	
	if sensitivity:
		y0 = numpy.zeros(Nspec+Nspec*Nrxn,numpy.float64)
		y0[0:Nspec] = C0
	else:
		y0 = C0
	
	# Create solver object
	model = MethaneCombustion(K, S, Sr, Nspec, Nrxn, sensitivity)
	dydt0 = - model.residual(t0, y0, numpy.zeros(Nspec+Nspec*Nrxn, numpy.float64))[0]
	model.initialize(t0, y0, dydt0, atol=1e-16, rtol=1e-12)
	
	# Initialize solution vectors
	t = []
	C = []
	Z = []
	
	# Generate the solution by stepping until tmax is reached
	# This will give you the solution at time points automatically selected
	# by the solver
	model.advance(tmin)
	iter = 0
	while iter < maxiter and model.t < tmax:
		model.step(tmax)
		t.append(model.t)
		# You must make a copy of y because it is overwritten by DASSL at
		# each call to step()
		C.append(model.y[0:Nspec].copy())
		if sensitivity:
			Z.append(numpy.reshape(model.y[Nspec:Nspec+Nspec*Nrxn].copy(), Nspec, Nrxn))
	# Convert the solution vectors to numpy arrays
	t = numpy.array(t, numpy.float64)
	C = numpy.array(C, numpy.float64)
	if sensitivity:
		Z = numpy.array(Z, numpy.float64)
	
	# If matplotlib is installed, show a plot of the results
	# Otherwise do nothing
	try:
		import pylab
		pylab.loglog(t, C)
		names = ['CH4', 'CH3', 'H', 'OH', 'H2O', 'H2', 'O', 'CH2O', 'O2', 'CHO', 'CO', 'CO2', 'HO2', 'Ar', 'M']
		pylab.xlabel('Time (s)')
		pylab.ylabel('Concentration (molecules/cm^3)')
		pylab.xlim(tmin, tmax)
		pylab.ylim(1e10, 1e20)
		for i, name in enumerate(names):
			pylab.annotate(name, (t[-1]*0.999,C[-1,i]), xytext=(20,-5), textcoords='offset points', arrowprops=dict(arrowstyle="-"))
			pylab.annotate(name, (t[0]*1.001,C[0,i]), xytext=(20,-5), textcoords='offset points', arrowprops=dict(arrowstyle="-"))
		pylab.show()
	except ImportError:
		pass
    