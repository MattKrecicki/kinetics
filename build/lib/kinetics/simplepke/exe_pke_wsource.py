# -*- coding: utf-8 -*-
"""exe_pke

Executes the PKE solver for constant reactivity insertion

"""

import numpy as np
from pke_w_source import PKEwSource
from kinetics.functions.plotters import lineplot


# Define kinetic parameters
# -----------------------------------------------------------------------------
beta = 0.00689 * np.array([0.033, 0.219, 0.196, 0.395, 0.115, 0.042])
lamda = np.array([0.0124, 0.0305, 0.1110, 0.3011, 1.1400, 3.0100])
promptL = 6E-05


#define external source
# -----------------------------------------------------------------------------
S0 = 37000000000          #1 Ci source
epsilon = 0.2             #source is only 30% efficient in causing fissions
Q = 200                   #MeV per fission


# Generate reactivity scenario
# -----------------------------------------------------------------------------
rhoi = -0.05  # initial negative reactivity
simulationTime = 800                                   # seconds
npoints = 500                                          # number of time points
timepoints = np.linspace(0, simulationTime, npoints)   # absolute time vector


# define external reactivity function
def externalrho(t):
    
    if t > 0.5 and t < 100.0:
        rhoext = 0.00049*t
    else:
        rhoext = 0.00049*100.0
        
    return rhoext 


# Execute the PKE solver
# -----------------------------------------------------------------------------
pke = PKEwSource(rhoi=rhoi, rho=externalrho, beta=beta, lamda=lamda, promptL=promptL,
                 S0=S0, Q=Q, epsilon=epsilon, timepoints=timepoints, rtol=1E-15)
pke.solve()


# Plot results
# -----------------------------------------------------------------------------

lineplot([timepoints]*1, [pke.power],
         xlabel="Time, sec", ylabel="Reactor power, Watts",
         markers=["None"]*1, linestyles=["--"]*1, colors=["black"],
         grid=True)


lineplot([timepoints[1:]]*1, [pke.totalrho[1:]],
         xlabel="Time, sec", ylabel="Total reactivity, dk/k",
         markers=["None"]*1, linestyles=["--"]*1, label=[r"$\rho$"],
         grid=True)


lineplot([timepoints]*6, [pke.dg1, pke.dg2, pke.dg3, pke.dg4, pke.dg5, pke.dg6],
         xlabel="Time, sec", ylabel="Delayed neutron precusor conc., a.u.",
         markers=["None"]*6, linestyles=["--"]*6, grid=True,
         label=["Group 1", "Group 2", "Group 3", "Group 4", "Group 5",
                 "Group 6"])


