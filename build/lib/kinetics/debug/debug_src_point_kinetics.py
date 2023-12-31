# -*- coding: utf-8 -*-
"""debug_point_kinetics.py

function debugs point kinetics solver

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""


import numpy as np
from kinetics.functions.pointkinetics import srcpke
from kinetics.functions.control import generalControlRule
from kinetics.containers.outputs import pointkineticscontainer
from kinetics.functions.plotters import lineplot


# define reactivity scenario
timepoints = np.linspace(0.0, 100.0, 500)

rhoext = generalControlRule(['linear', 'linear'],
                            [[0.0, 0.0], [0.0, 0.0]],
                            [0.5, 20.0])


# define point kinetics parameters
beta = 0.00689 * np.array([0.033, 0.219, 0.196, 0.395, 0.115, 0.042])
lamda = np.array([0.0124, 0.0305, 0.1110, 0.3011, 1.1400, 3.0100])
promptL = 6E-05 #seconds
S0 = 3.7e+10 #n/s, 1 Ci source
epsilon = 0.3 #source efficiency
rhoi = -0.05     #dk/k
volume = 1.0 #m3
nubar = 2.434 #n/fission
Q = 200.0 #MeV/fission
v = 1.86E+04 #m/s

#
nest = nubar * S0 * epsilon * (1 / (1+rhoi))

# initalize point kinetics solver
pkesolver = srcpke(beta=beta, lamda=lamda, promptL=promptL, S0=S0, rhoi=rhoi,
                   epsilon=epsilon, volume=volume, nubar=nubar, Q=Q, v=v,
                   rhoext=rhoext, timepoints=timepoints)

# execute solver
pkesolver.solve(rtol=1e-10)

#plot change in neutron population
lineplot([pkesolver.timepoints], [pkesolver.solution.power], markers=["None"],
         linestyles=["--"], grid=True, xlabel="Time, seconds",
         ylabel="Power, Watts")

#plot change in delayed neutron prescusors
lineplot([pkesolver.timepoints]*6,
         [pkesolver.solution.dnt[0,:], pkesolver.solution.dnt[1,:],
          pkesolver.solution.dnt[2,:], pkesolver.solution.dnt[3,:],
          pkesolver.solution.dnt[4,:], pkesolver.solution.dnt[5,:]],
         label=["Group 1", "Group 2", "Group 3",
                "Group 4", "Group 5", "Group 6"],
         markers=["None"]*6, loc="upper right",
         linestyles=["--"]*6, grid=True, xlabel="Time, seconds",
         ylabel="Delay neutron precusor conc., a.u.")

#plot total excess reactivity as a function of time
lineplot([pkesolver.timepoints], [pkesolver.solution.rho], markers=["None"],
         linestyles=["--"], grid=True, xlabel="Time, seconds",
         ylabel="Excess reactivity, dk/k")


# export results to hdf5 file
pkesolver.solution.export("srcpke.h5")

# test recovery of results
res = pointkineticscontainer()
res.recover("srcpke.h5")