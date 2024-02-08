# -*- coding: utf-8 -*-
"""debug_point_kinetics.py

function debugs point kinetics solver

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""


import numpy as np
from kinetics.functions.tworegionkinetics import srcpke2region
from kinetics.functions.control import controlrule
from kinetics.containers.inputs import pointkineticsInputsContainer
from kinetics.containers.outputs import pointkineticsOutputsContainer
from kinetics.functions.plotters import lineplot

# define reactivity scenario
timepoints = np.linspace(0.0, 100.0, 500)
rhoext = controlrule(['linear', 'linear', 'linear'],
                     [[0.0, 0.0], [0.0, 0.02], [0.0, 0.04]],
                     [5.0, 50.0, 105.0])

# define point kinetics parameters
beta = 0.00689 * np.array([0.033, 0.219, 0.196, 0.395, 0.115, 0.042])
lamda = np.array([0.0124, 0.0305, 0.1110, 0.3011, 1.1400, 3.0100])
promptLc = 38.6e-6
promptLr = 91.2e-6

# define transfer probablities
frc = 0.4877
fcr = 0.4880

rhoi = -0.05  #inital excess reactivity, dk/k
S0 = 3.7e+10  #external source strength n/s, 1 Ci source
epsilon = 0.3 #source efficiency
volumec = 1.0 #m3
volumer = 1.0 #m3
nubar = 2.434 #n/fission
Q = 200.0     #MeV/fission
vc = 1.86E+04 #m/s
vr = 1E+04    #m/s

# construct two-region kinetics inputs container
inputs = \
    pointkineticsInputsContainer(beta=beta, lamda=lamda, promptLc=promptLc, \
        promptLr=promptLr, frc=frc, fcr=fcr, S0=S0, epsilon=epsilon, rhoi=rhoi,
        volumec=volumec, volumer=volumer, nubar=nubar, Q=Q, vc=vc, vr=vr,
        rhoext=rhoext, timepoints=timepoints, typ="srcpke2region")

# initalize point kinetics solver
pkesolver = srcpke2region(inputs)

# execute solver
pkesolver.solve(rtol=1e-8)

# plot change in reactor power
lineplot([pkesolver.inputs.timepoints], [pkesolver.outputs.power], markers=["None"],
         linestyles=["--"], grid=True, xlabel="Time, seconds",
         ylabel="Power, Watts")

# plot change in neutron populations
lineplot([pkesolver.inputs.timepoints]*2,
         [pkesolver.outputs.ntc, pkesolver.outputs.ntr], markers=["None"]*2,
         linestyles=["--"]*2, label=["core", "reflector"], grid=True,
         xlabel="Time, seconds", ylabel="Neutron population, a.u.")

# plot change in delayed neutron prescusors
lineplot([pkesolver.inputs.timepoints]*6,
         [pkesolver.outputs.dnt[0,:], pkesolver.outputs.dnt[1,:],
          pkesolver.outputs.dnt[2,:], pkesolver.outputs.dnt[3,:],
          pkesolver.outputs.dnt[4,:], pkesolver.outputs.dnt[5,:]],
         label=["Group 1", "Group 2", "Group 3",
                "Group 4", "Group 5", "Group 6"],
         markers=["None"]*6, loc="upper left",
         linestyles=["--"]*6, grid=True, xlabel="Time, seconds",
         ylabel="Delay neutron precusor conc., a.u.")

# plot total excess reactivity as a function of time
lineplot([pkesolver.inputs.timepoints], [pkesolver.outputs.rho], markers=["None"],
         linestyles=["--"], grid=True, xlabel="Time, seconds",
         ylabel="Excess reactivity, dk/k")

# export results to hdf5 file
pkesolver.outputs.export("srcpke2region.h5")

# test recovery of results
res = pointkineticsOutputsContainer()
res.recover("srcpke2region.h5")
