# -*- coding: utf-8 -*-
"""debug_point_kinetics.py

function debugs point kinetics solver

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""


import numpy as np
from scipy.interpolate import interp1d
from kinetics.functions.inversepointkinetics import inversepke, inversepke_old
from kinetics.containers.inputs import pointkineticsInputsContainer
from kinetics.containers.outputs import pointkineticsOutputsContainer
from kinetics.functions.plotters import lineplot


# test recovery of results
res = pointkineticsOutputsContainer()
res.recover("pke.h5")


#generate power function
powerfunc = interp1d(res.timepoints, res.power, fill_value="extrapolate")
timepoints = np.linspace(0.0, 40.0, 1000000)

#kinetic parameters
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

print("solving")
#solver inverse kinetics problem
invpke.solve(rtol=1e-3, factor=1e-2)
print("solved")

intv = 100

#plot computed excess reactivity to actual excess reactivity
lineplot([res.timepoints[::intv], invpke.outputs.timepoints[::intv]],
         [res.rho[::intv], invpke.outputs.rho[::intv]],
         xlabel="time, seconds", ylabel="reactivity, dk/k",
         label=["ref.", "computed"], colors=["black", "red"],
         linestyles=["-", "--"], markers=["None", "None"], grid=True)

lineplot([res.timepoints[::intv]],
         [(invpke.outputs.rho[::intv] - res.rho[::intv])*1e+5],
         xlabel="time, seconds", ylabel="reactivity error, pcm",
         label=["error"], colors=["blue"],
         linestyles=["--"], markers=["None"], grid=True)



Cpke = res.dnt
Cinv = invpke.outputs.dnt

Cdiff = 100*(Cpke - Cinv)/Cpke

lineplot([res.timepoints[::intv]]*6,
         [Cinv[0,::intv], Cinv[1,::intv], Cinv[2,::intv], Cinv[3,::intv],
          Cinv[4,::intv], Cinv[5,::intv]],
         xlabel="time, seconds", ylabel="Rrecusor conc, a.u.",
         label=["Grp-1", "Grp-2", "Grp-3", "Grp-4", "Grp-5", "Grp-6"],
         linestyles=["--"]*6, markers=["None"]*6, grid=True)

lineplot([res.timepoints[::intv]]*6,
         [Cdiff[0,::intv], Cdiff[1,::intv], Cdiff[2,::intv], Cdiff[3,::intv],
          Cdiff[4,::intv], Cdiff[5,::intv]],
         xlabel="time, seconds", ylabel="Rrecusor conc. rel. difference, %",
         label=["Grp-1", "Grp-2", "Grp-3", "Grp-4", "Grp-5", "Grp-6"],
         linestyles=["--"]*6, markers=["None"]*6, grid=True)

