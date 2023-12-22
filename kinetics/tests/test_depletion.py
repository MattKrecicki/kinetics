# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 09:33:12 2021

@author: matt krecicki
"""

import h5py
import numpy as np
from ntptransient.depletion import exampleSolver
import matplotlib.pyplot as plt


with h5py.File("simpleXenon.h5", "r") as f:
    ref = f["results"]["Nt"][()]


n0 = [0.0,    0.0,    0.0,    0.0,      6.43230E-04*1E+24, 0.00257292*1E+24]
power = 330E+6
volume = 332097.8
Tc = 2700.0
Pc = 6.89E+6
timepoints = np.linspace(0, 30.0*60, 10) #30 min burn


sol = exampleSolver(n0, power, volume, Tc, Pc, timepoints, rtol=1E-10)


plt.scatter(timepoints/60, sol[:,0]/1E+24, label="I-135, NTP-TRANSIENT")
plt.scatter(timepoints/60, sol[:,1]/1E+24, label="Xe-135, NTP-TRANSIENT")
plt.plot(timepoints/60, ref[0,:], label="I-135, pyIsoDep")
plt.plot(timepoints/60, ref[1,:], label="Xe-135, pyIsoDep")
plt.xlabel("Time, minutes")
plt.ylabel("Concentration, atoms/b-cm")
plt.legend()
plt.show()
plt.close()


plt.scatter(timepoints/60, sol[:,2]/1E+24, label="Pm-149, NTP-TRANSIENT", c="C2")
plt.scatter(timepoints/60, sol[:,3]/1E+24, label="Sm-149, NTP-TRANSIENT", c="C3")
plt.plot(timepoints/60, ref[2,:], label="Pm-149, pyIsoDep", c="C2")
plt.plot(timepoints/60, ref[3,:], label="Sm-149, pyIsoDep", c="C3")
plt.xlabel("Time, minutes")
plt.ylabel("Concentration, atoms/b-cm")
plt.legend()
plt.show()
plt.close()