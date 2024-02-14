# -*- coding: utf-8 -*-
"""debug_point_kinetics.py

function debugs point kinetics solver

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""


import numpy as np
from scipy.interpolate import interp1d
from kinetics.functions.inversepointkinetics import inversepke
from kinetics.containers.inputs import pointkineticsInputsContainer
from kinetics.containers.outputs import pointkineticsOutputsContainer
from kinetics.functions.plotters import lineplot


# test recovery of results
res = pointkineticsOutputsContainer()
res.recover("pke.h5")


powerfunc = interp1d(res.timepoints, res.power)
timepoints = np.linspace(0.0, 100.0, 1000)


beta = 0.00689 * np.array([0.033, 0.219, 0.196, 0.395, 0.115, 0.042])
lamda = np.array([0.0124, 0.0305, 0.1110, 0.3011, 1.1400, 3.0100])
promptL = 6E-05 
volume = 1.0 #m3
nubar = 2.434
Q = 200.0
v = 1.86E+04 #m/s


#construct inputs container
inputs = \
    pointkineticsInputsContainer(beta=beta, lamda=lamda, promptL=promptL,\
        volume=volume, nubar=nubar, Q=Q, v=v, power=powerfunc,
        timepoints=timepoints, typ="invpke")
        

#initialize inverse point kinetics solver
invpke = inversepke(inputs)


#solver inverse kinetics problem
invpke.solve()




