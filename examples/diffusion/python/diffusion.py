#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

"""

import numpy
from model import DiffusionModel

if __name__ == '__main__':
    
    # The times at which to obtain the solution
    tlist = numpy.array([10**i for i in range(-6, 1)], numpy.float64)
    # The number of grid points to use to discretize the spatial dimension
    N = 501
    
    # Set initial conditions
    t0 = 0.0
    y0 = numpy.zeros((N), numpy.float64)
    y0[0] = 1
    
    # Initialize the model
    model = DiffusionModel(N=N)
    dydt0 = - model.residual(t0, y0, numpy.zeros((N), numpy.float64))[0]
    model.initialize(t0, y0, dydt0)
    
    # Integrate to get the solution at each time point
    t = []; y = []
    for t1 in tlist:
        model.advance(t1)
        t.append(model.t)
        # You must make a copy of y because it is overwritten by DASSL at
        # each call to advance()
        y.append(model.y.copy())

    # Convert the solution vectors to numpy arrays
    t = numpy.array(t, numpy.float64)
    y = numpy.array(y, numpy.float64)

    x = numpy.linspace(0, 1, N)
    
    # Generate plot of solutions
    import pylab
    pylab.figure(figsize=(7,5))
    pylab.plot(x, y.T)
    pylab.xlim((0,1))
    pylab.ylim((0,1))
    pylab.xlabel('Membrane position $x$')
    pylab.ylabel('Concentration $\\theta$')
    pylab.legend([
        't = $\\mathdefault{10^{-6}}$', 
        't = $\\mathdefault{10^{-5}}$', 
        't = $\\mathdefault{10^{-4}}$', 
        't = $\\mathdefault{10^{-3}}$', 
        't = $\\mathdefault{10^{-2}}$', 
        't = $\\mathdefault{10^{-1}}$', 
        't = $\\mathdefault{10^{0}}$'], loc=1)
    pylab.show()
    