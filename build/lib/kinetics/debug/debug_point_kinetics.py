# -*- coding: utf-8 -*-
"""debug_point_kinetics.py

function debugs point kinetics solver

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""

import numpy as np
from kinetics.functions.pointkinetics import pointkinetics
from kinetics.functions.control import generalControlRule


#define reactivity scenario
timepoints = np.linspace(0.0, 100.0, 50)
rhoext = generalControlRule(['linear'], [[0.0, 0.0]], [105.0])


#define point kinetics instance
beta = 0.0075 * np.array([0.033, 0.219, 0.196, 0.395, 0.115, 0.042])
lamda = np.array([0.0124, 0.0305, 0.1110, 0.3011, 1.1400, 3.0100])
promptL = 6E-05
P0 = 1.0
volume = 1.0
nubar = 2.434
Q = 200.0


#initalize point kinetics solver
pke = \
    pointkinetics(beta=beta, lamda=lamda, promptL=promptL, P0=P0, volume=volume,
                  nubar=nubar, Q=Q, rhoext=rhoext)

#execute solver
pke.solve()

#export results to hdf5 file


#test recovery of results