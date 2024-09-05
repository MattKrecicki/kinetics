# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 17:17:59 2024

@author: matt krecicki
@email: matthewkrecicki@gmail.com

References:
    1. Verification of iterative matrix solutions for multipoint kinetics equations
    2. Investigating decoupling effects in CFV type of reactors using Avery’s coupled reactors theory

Need to move the dependencies to the region object and out of the container object

container object should have the interpolation functions within it and the ablity to generate the matrix


"""


import numpy as np
from kinetics.containers.inputs import regionKineticsData
from kinetics.containers.inputs import multiPointKineticsInputsContainer as \
    inputsContainer
from kinetics.functions.multipointkinetics import avery as transientsolver

# -----------------------------------------------------------------------------
# ----- Setup Multipoint Kinetics Data Container for Transient Analysis -------
# -----------------------------------------------------------------------------

mpkdata = \
    inputsContainer(typ="avery", nregions=3, ndelayed=1,
                    dependencies=["hrod"], order=["1", "2", "3"],
                    detectors=["N", "S", "E", "W"])


# when defining regions the Id will act as the fission site and the defined 
# value will act as the fission sink. Therefore the Kjk defined in the first 
# region is defined as k11, k21, k31


# ----- initialize region 1 data ----------------------------------------------
r1 = \
    regionKineticsData(Id="1", typ="avery", x=0.0, y=-1.0, z=0.0, 
                       dependencies=["hrod"], volume=1.0)


#add region 1's 0cm rod height function

r1.add(#average energy released per fission
       Q=200.0,
       #one-group neutron velocity
       v=1.86E+04,        
       
       #detail how couling is generated so matrix can be constructed correctly
       coupling = [    "1",     "2",      "3"],
       
       # prompt fission coupling coefficients (sink < src)
       # coupling to:   1 < 1    2 < 1    3 < 1
       Kjk =  np.array([0.92804, 0.03754, 0.03756]),  
       Kjkd = np.array([0.92804, 0.03754, 0.03756]),
       
       # mean generation times (sink < src)
       # coupling to:   1 < 1   2 < 1   3 < 1
       Ljk =  np.array([0.3613, 0.6126, 0.6154])*1e-6,
       Ljkd = np.array([0.3613, 0.6126, 0.6154])*1e-6,
       
       #delayed neutron group data, only a single group included
       # coupling to:           1 < 1     2 < 1     3 < 1
       Bjk =      {1: np.array([0.003760, 0.003639, 0.003657])},
       lamdajk =  {1: np.array([0.40529,  0.40529,  0.40529])},
       
       #specific detector response functions
       detectors={"N": {"kp": 0.1, "lp": 10e-6, "kd": 0.1, "ld": 10e-6},
                  "S": {"kp": 1.0, "lp": 10e-6, "kd": 1.0, "ld": 10e-6},
                  "E": {"kp": 0.2, "lp": 10e-6, "kd": 0.2, "ld": 10e-6},
                  "W": {"kp": 0.2, "lp": 10e-6, "kd": 0.2, "ld": 10e-6}},
       
       #specify dependency variables
       hrod=0.0)


#add region 1's 35cm rod height function
r1.add(#average energy released per fission
       Q=200.0,
       #one-group neutron velocity
       v=1.86E+04,        
       
       #detail how couling is generated so matrix can be constructed correctly
       coupling = [    "1",     "2",      "3"],
       
       # prompt fission coupling coefficients (sink < src)
       # coupling to:   1 < 1    2 < 1    3 < 1
       Kjk =  np.array([0.95226, 0.03820, 0.03616]),  
       Kjkd = np.array([0.95226, 0.03820, 0.03616]),
       
       # mean generation times (sink < src)
       # coupling to:   1 < 1   2 < 1   3 < 1
       Ljk =  np.array([0.3690, 0.6115, 0.6162])*1e-6,
       Ljkd = np.array([0.3690, 0.6115, 0.6162])*1e-6,
       
       #delayed neutron group data, only a single group included
       # coupling to:           1 < 1     2 < 1     3 < 1
       Bjk =      {1: np.array([0.003766, 0.003633, 0.003651])},
       lamdajk =  {1: np.array([0.40529,  0.40529,  0.40529])},
       
       #specific detector response functions
       detectors={"N": {"kp": 0.1, "lp": 10e-6, "kd": 0.1, "ld": 10e-6},
                  "S": {"kp": 1.0, "lp": 10e-6, "kd": 1.0, "ld": 10e-6},
                  "E": {"kp": 0.2, "lp": 10e-6, "kd": 0.2, "ld": 10e-6},
                  "W": {"kp": 0.2, "lp": 10e-6, "kd": 0.2, "ld": 10e-6}},
       
       #specify dependency variables
       hrod=35.0)

#add 1st region to data container
mpkdata.add(r1)

# ----- initialize region 2 data ----------------------------------------------
r2 = \
    regionKineticsData(Id="2", typ="avery", x=1.0, y=0.0, z=0.0, 
                       dependencies=["hrod"], volume=1.0)


#add region 2's 0cm rod height function
r2.add(#average energy released per fission
       Q=200.0,
       #one-group neutron velocity
       v=1.86E+04,        
       
       #detail how couling is generated so matrix can be constructed correctly
       coupling = [    "1",     "2",      "3"],
       
       # prompt fission coupling coefficients (sink < src)
       # coupling to:   1 < 2    2 < 2    3 < 2
       Kjk =  np.array([0.03756, 0.92804, 0.03754]),  
       Kjkd = np.array([0.03756, 0.92804, 0.03754]),
       
       # mean generation times (sink < src)
       # coupling to:   1 < 2   2 < 2   3 < 2
       Ljk =  np.array([0.6154, 0.3613, 0.6126])*1e-6,
       Ljkd = np.array([0.6154, 0.3613, 0.6126])*1e-6,
       
       #delayed neutron group data, only a single group included
       # coupling to:           1 < 2     2 < 2     3 < 2
       Bjk =      {1: np.array([0.003657, 0.003760, 0.003639])},
       lamdajk =  {1: np.array([0.40529,  0.40529,  0.40529])},
       
       #specific detector response functions
       detectors={"N": {"kp": 0.1, "lp": 10e-6, "kd": 0.1, "ld": 10e-6},
                  "S": {"kp": 1.0, "lp": 10e-6, "kd": 1.0, "ld": 10e-6},
                  "E": {"kp": 0.2, "lp": 10e-6, "kd": 0.2, "ld": 10e-6},
                  "W": {"kp": 0.2, "lp": 10e-6, "kd": 0.2, "ld": 10e-6}},
       
       #specify dependency variables
       hrod=0.0)


#add region 2's 35cm rod height function
r2.add(#average energy released per fission
       Q=200.0,
       #one-group neutron velocity
       v=1.86E+04,        
       
       #detail how couling is generated so matrix can be constructed correctly
       coupling = [    "1",     "2",      "3"],
       
       # prompt fission coupling coefficients (sink < src)
       # coupling to:   1 < 2    2 < 2    3 < 2
       Kjk =  np.array([0.03766, 0.92029, 0.03488]),  
       Kjkd = np.array([0.03766, 0.92029, 0.03488]),
       
       # mean generation times (sink < src)
       # coupling to:   1 < 2   2 < 2   3 < 2
       Ljk =  np.array([0.6161, 0.3588, 0.6178])*1e-6,
       Ljkd = np.array([0.6161, 0.3588, 0.6178])*1e-6,
       
       #delayed neutron group data, only a single group included
       # coupling to:           1 < 2     2 < 2     3 < 2
       Bjk =      {1: np.array([0.003651, 0.003756, 0.003632])},
       lamdajk =  {1: np.array([0.40529,  0.40529,  0.40529])},
       
       #specific detector response functions
       detectors={"N": {"kp": 0.1, "lp": 10e-6, "kd": 0.1, "ld": 10e-6},
                  "S": {"kp": 1.0, "lp": 10e-6, "kd": 1.0, "ld": 10e-6},
                  "E": {"kp": 0.2, "lp": 10e-6, "kd": 0.2, "ld": 10e-6},
                  "W": {"kp": 0.2, "lp": 10e-6, "kd": 0.2, "ld": 10e-6}},
       
       #specify dependency variables
       hrod=35.0)

#add 2nd region to data container
mpkdata.add(r2)

# ----- initialize region 3 data ----------------------------------------------

r3 = \
    regionKineticsData(Id="3", typ="avery", x=-1.0, y=0.0, z=0.0, 
                       dependencies=["hrod"], volume=1.0)
    

#add region 3's 0cm rod height function
r3.add(#average energy released per fission
       Q=200.0,
       #one-group neutron velocity
       v=1.86E+04,        
       
       #detail how couling is generated so matrix can be constructed correctly
       coupling = [    "1",     "2",      "3"],
       
       # prompt fission coupling coefficients (sink < src)
       # coupling to:   1 < 3    2 < 3    3 < 3
       Kjk =  np.array([0.03754, 0.03756, 0.92804]),  
       Kjkd = np.array([0.03754, 0.03756, 0.92804]),
       
       # mean generation times (sink < src)
       # coupling to:   1 < 3   2 < 3   3 < 3
       Ljk =  np.array([0.6126, 0.6154, 0.3613])*1e-6,
       Ljkd = np.array([0.6126, 0.6154, 0.3613])*1e-6,
       
       #delayed neutron group data, only a single group included
       # coupling to:           1 < 3     2 < 3     3 < 3
       Bjk =      {1: np.array([0.003639, 0.003657, 0.003760])},
       lamdajk =  {1: np.array([0.40529,  0.40529,  0.40529])},
       
       #specific detector response functions
       detectors={"N": {"kp": 0.1, "lp": 10e-6, "kd": 0.1, "ld": 10e-6},
                  "S": {"kp": 1.0, "lp": 10e-6, "kd": 1.0, "ld": 10e-6},
                  "E": {"kp": 0.2, "lp": 10e-6, "kd": 0.2, "ld": 10e-6},
                  "W": {"kp": 0.2, "lp": 10e-6, "kd": 0.2, "ld": 10e-6}},
       
       #specify dependency variables
       hrod=0.0)


#add region 3's 35cm rod height function
r3.add(#average energy released per fission
       Q=200.0,
       #one-group neutron velocity
       v=1.86E+04,        
       
       #detail how couling is generated so matrix can be constructed correctly
       coupling = [    "1",     "2",      "3"],
       
       # prompt fission coupling coefficients (sink < src)
       # coupling to:   1 < 3    2 < 3    3 < 3
       Kjk =  np.array([0.03561, 0.03489, 0.92116]),  
       Kjkd = np.array([0.03561, 0.03489, 0.92116]),
       
       # mean generation times (sink < src)
       # coupling to:   1 < 3   2 < 3   3 < 3
       Ljk =  np.array([0.6145, 0.6214, 0.3595])*1e-6,
       Ljkd = np.array([0.6145, 0.6214, 0.3595])*1e-6,
       
       #delayed neutron group data, only a single group included
       # coupling to:           1 < 3     2 < 3     3 < 3
       Bjk =      {1: np.array([0.003633, 0.003649, 0.003756])},
       lamdajk =  {1: np.array([0.40529,  0.40529,  0.40529])},
       
       #specific detector response functions
       detectors={"N": {"kp": 0.1, "lp": 10e-6, "kd": 0.1, "ld": 10e-6},
                  "S": {"kp": 1.0, "lp": 10e-6, "kd": 1.0, "ld": 10e-6},
                  "E": {"kp": 0.2, "lp": 10e-6, "kd": 0.2, "ld": 10e-6},
                  "W": {"kp": 0.2, "lp": 10e-6, "kd": 0.2, "ld": 10e-6}},
       
       #specify dependency variables
       hrod=35.0)

#add 3rd region to data container
mpkdata.add(r3)

#test region interpolator
#tst = r3.evaluate(hrod=30.0)


# -----------------------------------------------------------------------------
# ----- Assemble final container ----------------------------------------------
# -----------------------------------------------------------------------------

mpkdata.validate()

#tstData = mpkdata.evaluate(hrod=30)


# -----------------------------------------------------------------------------
# ----- Solve Initial Conditions ----------------------------------------------
# -----------------------------------------------------------------------------

#detail time points transient solution should be returned
timepoints = np.linspace(0.0, 10.0, 100)

#define function that returns rod height as a function of time
def hrodfunc(t):
    return 0.0


#initialize transient solver
model = \
    transientsolver(timepoints=timepoints, kineticdata=mpkdata, hrod=hrodfunc,
                    hrod0=0.0, verbose=True)

#only solve initial conditions
model.solve(srcepi=1e-10, keffepi=1e-6, maxitr=100, transtol=1e-10,
            initialonly=False)


