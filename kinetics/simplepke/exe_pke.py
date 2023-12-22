# -*- coding: utf-8 -*-
"""exe_pke

Executes the PKE solver for constant reactivity insertion

"""

import numpy as np

from solvepke import SolvePKEs
from kinetics.functions.plotters import lineplot

# Define kinetic parameters
# -----------------------------------------------------------------------------
beta = 0.0075 * np.array([0.033, 0.219, 0.196, 0.395, 0.115, 0.042])
lamda = np.array([0.0124, 0.0305, 0.1110, 0.3011, 1.1400, 3.0100])
promptL = 6E-05
P0 = 1  # relative power


# Generate external reactivity scenario
# -----------------------------------------------------------------------------
# constant reactivity insertion
simulationTime = 10  # seconds
npoints = 50  # number of time points
timepoints = np.linspace(0, simulationTime, npoints)  # absolute time vector
rho = 0.1 * beta.sum()

# Execute the PKE solver
# -----------------------------------------------------------------------------
xt = SolvePKEs(rho, beta, lamda, promptL, P0, timepoints, rtol=1E-10)

# Plot results
# -----------------------------------------------------------------------------

lineplot([timepoints], [xt[:,0]], xlabel="Time, sec", ylabel="Normalized Flux",
         markers=["s"], linestyles=["--"], grid=True)

#plot1D(timepoints, xt[:, 0], xlabel="Time, sec", ylabel="Normalized Flux",
#       fontsize=8, marker="*r", markerfill=False, markersize=6)
