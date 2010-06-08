#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This file investigates the methane oxidation kinetic mechanism found in

    A. A. Boni and R. C. Penner. "Sensitivity Analysis of a Mechanism for 
    Methane Oxidation Kinetics." Comb. Sci. Tech. 15, p. 99-106 (1977).

The mechanism can be integrated in time with various initial conditions with
or without sensitivity analysis.
"""

import sys
sys.path.append('.')

from PyDAS import dassl
import numpy
import cython

################################################################################

class MethaneCombustion(dassl.DASSL):
	"""
	"""
	
	def __init__(self, K, S, Sr, Nspec, Nrxn, sensitivity):
		self.K = K
		self.S = S
		self.Sr = Sr
		self.Nspec = Nspec
		self.Nrxn = Nrxn
		self.sensitivity = sensitivity
	
	def residual(self, t, y, dydt):
		
		cython.declare(i=cython.int, j=cython.int, delta=numpy.ndarray)
		cython.declare(Nspec=cython.int, Nrxn=cython.int)
		
		delta = numpy.zeros(y.shape[0], numpy.float64)
		
		Nspec = self.Nspec; Nrxn = self.Nrxn
		
		# Separate y into concentration and sensitivity coefficient
		# components
		C = y[0:Nspec]
		if self.sensitivity:
			Z = numpy.reshape(y[Nspec:Nspec+Nspec*Nrxn], Nspec, Nrxn)
		
		# Calculate the species fluxes
		rxnRates = numpy.zeros(Nrxn, numpy.float64)
		for j in range(Nrxn):
			rxnRates[j] = self.K[j] * numpy.prod(C ** self.Sr[:,j])
		delta[0:Nspec] = numpy.dot(self.S, rxnRates)
		# Set Ar and N2 concentrations to be constant
		delta[13] = 0.0
		delta[14] = 0.0
		for i in range(Nspec):
			delta[i] -= dydt[i]
		
		return delta, 0
	
	#def jacobian(self, t, y, dydt, cj):
		
		#pd = -cj * numpy.identity(y.shape[0], numpy.float64)
		
		#Nspec = self.Nspec; Nrxn = self.Nrxn
		
		## Separate y into concentration and sensitivity coefficient
		## components
		#C = y[0:Nspec]
		#if self.sensitivity:
			#Z = numpy.reshape(y[Nspec:Nspec+Nspec*Nrxn], Nspec, Nrxn)
		
		#for i in range(self.Nspec):
			#for l in range(self.Nspec):
				#sk = numpy.zeros(self.Nspec, numpy.int)
				#sk[l] = 1
				#if C[l] > 0:
					#for j in range(self.Nrxn):
						#pd[i,l] += self.S[i,j] * self.K[j] * self.Sr[l,j] * numpy.prod(C**(self.Sr[:,j] - sk))
		#return pd
	
