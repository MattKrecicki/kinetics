# -*- coding: utf-8 -*-
"""

function debugs point kinetics input container

Created on Sun Dec 31 13:06:47 2023

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""


import numpy as np
from kinetics.containers.inputs import pointkineticsContainer
from kinetics.functions.control import generalControlRule


# define reactivity scenario
timepoints = np.linspace(0.0, 100.0, 500)

rhoext = generalControlRule(['linear', 'linear'],
                            [[0.0, 0.0], [0.0, 0.003]],
                            [0.5, 20.0])


# define point kinetics parameters
beta = 0.00689 * np.array([0.033, 0.219, 0.196, 0.395, 0.115, 0.042])
lamda = np.array([0.0124, 0.0305, 0.1110, 0.3011, 1.1400, 3.0100])
promptLc = 38.6e-6
promptLr = 91.2e-6

#define transfer probablities
frc = 0.4877
fcr = 0.4880

P0 = 1e+3 #initial reactor power
volumec = 1.0 #m3
volumer = 1.0 #m3
nubar = 2.434 #n/fission
Q = 200.0 #MeV/fission
vc = 1.86E+04 #m/s
vr = 1E+04 #m/s


inputs = \
    pointkineticsContainer(beta=beta, lamda=lamda, promptLc=promptLc, \
        promptLr=promptLr, frc=frc, fcr=fcr, P0=P0, volumec=volumec, 
        volumer=volumer, nubar=nubar, Q=Q, vc=vc, vr=vr, rhoext=rhoext,
        timepoints=timepoints, typ="pke2region")
        

inputs.whatis("P0")