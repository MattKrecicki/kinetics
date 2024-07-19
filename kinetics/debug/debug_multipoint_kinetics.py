# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 17:17:59 2024

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""

import numpy as np
from kinetics.containers.inputs import multiPointKineticsInputsContainer as \
    inputsContainer
from kinetics.functions.multipointkinetics import avery


# ----- input kinetic parameters from Valocchi et al., 2020 "Reduced order
# models in reactor kinetics: A comparison between point kinetics and
# multipoint kinetics

inputs = inputsContainer(typ="avery")

#add region 1
inputs.add(Id="1",
           typ="avery",
           Q=200.0,
           v=1.86E+04, #placeholder 
           volume=1.0, #placeholder
           lamdaki = np.array([0.40529]),
           Bki = np.array([0.003327]),
           # coupling to:  1      2       3      4
           Bjk = np.array([332.7, 1097.1, 496.7, 298.4])*1e-5,
           # coupling to:  1      2       3      4
           Kjk = np.array([0.89357, 0.28779, 0.15006, 0.05032]),
           # coupling to:  1      2       3      4
           Ljk = np.array([0.2465, 0.3491, 0.3935, 0.4495])*1e-6)


#add region 2
inputs.add(Id="2",
           typ="avery",
           Q=200.0,
           v=1.86E+04, #placeholder
           volume=1.0, #placeholder
           lamdaki = np.array([0.40529]),
           Bki = np.array([0.002899]),
           # coupling to:  1      2      3      4
           Bjk = np.array([119.3, 289.9, 170.0, 194.4])*1e-5,
           # coupling to:  1        2        3        4
           Kjk = np.array([0.03211, 0.12127, 0.07364, 0.02286]),
           # coupling to:  1       2       3       4
           Ljk = np.array([0.2938, 0.1370, 0.2646, 0.6829])*1e-6)


#add region 3
inputs.add(Id="3",
           typ="avery",
           Q=200.0,
           v=1.86E+04, #placeholder
           volume=1.0, #placeholder
           lamdaki = np.array([0.40529]),
           Bki = np.array([0.002899]),
           # coupling to:  1      2       3      4
           Bjk = np.array([339.8, 1486.1, 776.2, 714.7])*1e-5,
           # coupling to:  1        2        3        4
           Kjk = np.array([0.12139, 0.28695, 0.38502, 0.19824]),
           # coupling to:  1       2       3       4
           Ljk = np.array([2.8120, 2.5776, 2.4007, 4.2755])*1e-6)


#add region 4
inputs.add(Id="4",
           typ="avery",
           Q=200.0,
           v=1.86E+04, #placeholder
           volume=1.0, #placeholder
           lamdaki = np.array([0.40529]),
           Bki = np.array([0.002899]),
           # coupling to:  1      2       3      4
           Bjk = np.array([307.7, 1390.7, 738.7, 809.4])*1e-5,
           # coupling to:  1        2        3        4
           Kjk = np.array([0.16864, 0.43211, 0.58580, 0.75710]),
           # coupling to:  1        2        3        4
           Ljk = np.array([17.8861, 18.8878, 19.1796, 25.4637])*1e-6)


#run final error checking to make sure each region is coupled correctly
inputs.validate()


# ----- solve for initial conditions

mpk = avery(kineticData=inputs, P0=100.0)



#mpk.solve(rtol=1e-10, k0=1.0, flux0=None, maxitr=500)
