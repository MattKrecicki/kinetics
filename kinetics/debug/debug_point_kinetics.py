# -*- coding: utf-8 -*-
"""debug_point_kinetics.py

function debugs point kinetics solver

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""


import numpy as np
from kinetics.functions.pointkinetics import pke
from kinetics.functions.control import generalControlRule
from kinetics.containers.inputs import pointkineticsInputsContainer
from kinetics.containers.outputs import pointkineticsOutputsContainer
from kinetics.functions.plotters import lineplot


# define reactivity scenario
timepoints = np.linspace(0.0, 100.0, 500)

rhoext = generalControlRule(['linear', 'linear'],
                            [[0.0, 0.0], [0.0, 0.003445]],
                            [0.5, 10.0])

# define point kinetics parameters
beta = 0.00689 * np.array([0.033, 0.219, 0.196, 0.395, 0.115, 0.042])
lamda = np.array([0.0124, 0.0305, 0.1110, 0.3011, 1.1400, 3.0100])
promptL = 6E-05 
P0 = 1e+3 #watts
volume = 1.0 #m3
nubar = 2.434
Q = 200.0
v = 1.86E+04 #m/s

#construct inputs container
inputs = \
    pointkineticsInputsContainer(beta=beta, lamda=lamda, promptL=promptL, P0=P0,\
        volume=volume, nubar=nubar, Q=Q, v=v, rhoext=rhoext,
        timepoints=timepoints, typ="pke")

# initalize point kinetics solver
pkesolver = pke(inputs)

# execute solver
pkesolver.solve()

#plot change in neutron population
lineplot([pkesolver.inputs.timepoints],
         [pkesolver.outputs.power], markers=["None"],
         linestyles=["--"], grid=True, xlabel="Time, seconds",
         ylabel="Power, Watts")

#plot change in delayed neutron prescusors
lineplot([pkesolver.inputs.timepoints]*6,
         [pkesolver.outputs.dnt[0,:], pkesolver.outputs.dnt[1,:],
          pkesolver.outputs.dnt[2,:], pkesolver.outputs.dnt[3,:],
          pkesolver.outputs.dnt[4,:], pkesolver.outputs.dnt[5,:]],
         label=["Group 1", "Group 2", "Group 3",
                "Group 4", "Group 5", "Group 6"],
         markers=["None"]*6, loc="upper right",
         linestyles=["--"]*6, grid=True, xlabel="Time, seconds",
         ylabel="Delay neutron precusor conc., a.u.")

#plot total excess reactivity as a function of time
lineplot([pkesolver.inputs.timepoints], [pkesolver.outputs.rho],
         markers=["None"], linestyles=["--"], grid=True, xlabel="Time, seconds",
         ylabel="Excess reactivity, dk/k")


# export results to hdf5 file
pkesolver.outputs.export("pke.h5")

# test recovery of results
res = pointkineticsOutputsContainer()
res.recover("pke.h5")