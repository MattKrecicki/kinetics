# -*- coding: utf-8 -*-
"""debug_point_kinetics.py

function debugs point kinetics solver

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""


import numpy as np
from scipy.interpolate import interp1d
from kinetics.functions.inversepointkinetics import srcinversepke
from kinetics.containers.inputs import pointkineticsInputsContainer
from kinetics.containers.outputs import pointkineticsOutputsContainer
from kinetics.functions.plotters import lineplot


# test recovery of results
res = pointkineticsOutputsContainer()
res.recover("srcpke.h5")

#generate power function
powerfunc = interp1d(res.timepoints, res.power)
timepoints = np.linspace(0.0, 100.0, 300)
rho0 = -0.05

#kinetic parameters
beta = 0.00689 * np.array([0.033, 0.219, 0.196, 0.395, 0.115, 0.042])
lamda = np.array([0.0124, 0.0305, 0.1110, 0.3011, 1.1400, 3.0100])
promptL = 6E-05 
volume = 1.0 #m3
nubar = 2.434
Q = 200.0
v = 1.86E+04 #m/s

#source parameters
S0 = 3.7e+10    #n/s, 1 Ci source
epsilon = 0.3   #source efficiency


#construct inputs container
inputs = \
    pointkineticsInputsContainer(beta=beta, lamda=lamda, promptL=promptL,\
        rho0=rho0, volume=volume, nubar=nubar, Q=Q, v=v, S0=S0, epsilon=epsilon,
        power=powerfunc, timepoints=timepoints, typ="srcinvpke")
        

#initialize inverse point kinetics solver
srcinvpke = srcinversepke(inputs)


#solver inverse kinetics problem
srcinvpke.solve()


#plot computed excess reactivity to actual excess reactivity
lineplot([res.timepoints, srcinvpke.outputs.timepoints],
         [res.rho, srcinvpke.outputs.rhototal],
         xlabel="time, seconds", ylabel="reactivity, dk/k",
         label=["ref.", "computed"], colors=["black", "red"],
         linestyles=["-", "--"], markers=["None", "None"], grid=True)

err = srcinvpke.outputs.rhototal - res.rho

lineplot([res.timepoints],
         [err],
         xlabel="time, seconds", ylabel="reactivity error, dk/k",
         label=["error"], colors=["blue"],
         linestyles=["None"], markers=["s"], grid=True)



