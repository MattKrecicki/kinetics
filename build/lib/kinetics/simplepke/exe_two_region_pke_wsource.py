# -*- coding: utf-8 -*-
"""exe_pke

Executes the two-region PKE solver w/ external source

"""

import numpy as np
from two_region_pke_w_source import twoRegionPKEwSource
from kinetics.functions.plotters import lineplot


# Define kinetic parameters
# -----------------------------------------------------------------------------

#define core kinetic parameters
beta = 0.00689 * np.array([0.033, 0.219, 0.196, 0.395, 0.115, 0.042])
lamda = np.array([0.0124, 0.0305, 0.1110, 0.3011, 1.1400, 3.0100])
promptLc = 38.6e-6
promptLr = 91.2e-6

#define transfer probablities
f = 0.2380
frc = 0.4877
fcr = 0.4880

#define external source
# -----------------------------------------------------------------------------
S0 = 37000000000          #1 Ci source, n/s
epsilon = 0.2             #source is only 30% efficient in causing fissions
Q = 200.0                 #MeV per fission


# Generate reactivity scenario
# -----------------------------------------------------------------------------
rhoi = -0.05 # initial negative reactivity
simulationTime = 9000                                   # seconds
npoints = 500                                          # number of time points
timepoints = np.linspace(0, simulationTime, npoints)   # absolute time vector


# define external reactivity function
def externalrho(t):
    
    rhoadd = 0.04995
    dt = 1.0
    tstart = 100.0
    rate = rhoadd / dt
    
    if t > tstart and t < dt+tstart:
        rhoext = rate * (t-tstart)
        
    elif t >= dt+tstart:
        rhoext = rhoadd
        
    else:
        rhoext = 0.0
    
            
    return rhoext 


# Execute the PKE solver
# -----------------------------------------------------------------------------
pke = twoRegionPKEwSource(rhoi=rhoi, rho=externalrho, beta=beta, lamda=lamda,
                          promptLc=promptLc, promptLr=promptLr, f=f, frc=frc,
                          fcr=fcr, S0=S0, Q=Q, epsilon=epsilon,
                          timepoints=timepoints)
pke.solve(rtol=1e-5)


# Plot results
# -----------------------------------------------------------------------------

lineplot([timepoints]*2, [pke.Nc, pke.Nr],
         xlabel="Time, sec", ylabel="Neutron Population, a.u.",
         markers=["None"]*2, linestyles=["--"]*2, label=["Core", "Reflector"],
         grid=True)

lineplot([timepoints]*1, [pke.totalrho],
         xlabel="Time, sec", ylabel="Total excess reactivity, dkk",
         markers=["None"]*1, linestyles=["--"]*1, label=["Total", "Source"],
         grid=True)

lineplot([timepoints]*6, [pke.dg1, pke.dg2, pke.dg3, pke.dg4, pke.dg5, pke.dg6],
         xlabel="Time, sec", ylabel="Delayed neutron precusor conc., a.u.",
         markers=["None"]*6, linestyles=["--"]*6, grid=True,
         label=["Group 1", "Group 2", "Group 3", "Group 4", "Group 5",
                 "Group 6"])


