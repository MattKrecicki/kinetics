# -*- coding: utf-8 -*-
"""debug_steadyinputs

Example steady state for debugging purposes

"""

import numpy as np

from ntpSystem.containers.steadyinputs import SteadyStateInputs



#------------------------------------------------------------------------------
#----- define inputs
#------------------------------------------------------------------------------


#initalize input container-----------------------------------------------------
ssInps = SteadyStateInputs()


# core geometry----------------------------------------------------------------
ssInps.AddData("core",
    noFuelElems = '127',          # total number of fuel elements in reactor core
    fuelL = 1.20,               # active core height of fuel elements
    fuelFlowA = 2.0253e-04,     # flow area of a single fuel element
    fuelHper = 3.0687e-01,      # heat perimeter of a single fuel element
    noModeratorElems = 390,     # total number of moderator elements in core
    supplyL = 1.20,             # active core height of supply channel
    supplyFlowA = 1.2566e-05,   # flow area of a single supply channel
    supplyHper = 1.2566e-02,    # heated perimenter of a single supply channel
    returnL = 1.20,             # active core height of return channel
    returnFlowA = 5.8113e-05,   # flow area of a single return channel
    returnHper = 6.1255e-02,    # heated perimenter of a single return channel
    )  # core

# reactor geometry-------------------------------------------------------------
ssInps.AddData("reactor",
    refL = 1.2,                 # length of reflector coolant channels
    noDrums = 18,               # number of control drums in reflector 
    refNoChans = 100,           # number of coolant channels in reflector 
    refR = 0.005,               # heated perimenter of a single return channel       
    nozzleL = 2.65,             # length of nozzle coolant channels
    nozzleNoChans = 232,        # number of coolant channels in nozzle
    nozzleR = 0.00257,          # nozzle coolant channel radius
    nozzleExpR = 300.,          # nozzle expansion ratio
    At = 0.005662964739541421,  # nozzle throat area, m2
    pipeL = [0.1, 0.1],         # total piping length, m
    pipeR = [0.105, 0.105],     # piping radius, m
    )  # reactor

# operational conditions-------------------------------------------------------
ssInps.AddData("operation",
    Pc = 6.89E+06,              # desired chamber pressure, Pa
    Tc = 2700,                  # desired chamber temperature, K
    Pt = 0.4e+6,               # tank pressure, Pa
    Tt = 20.0,                  # tank temperature, K
    power = None,               # absolute total power in watts
    mdot = None,                # system mass flow rate kg/s
    mdotSplitRatio = 0.5,       # fraction of total flow sent to reflector
    )  # operation

# power fractions--------------------------------------------------------------
ssInps.AddData("powerfrac",
    pfFuel = 0.986602-0.05,     # fraction of total power in fuel
    pfSupply = 0.008970,        # fraction of total power in supply 
    pfReturn = 0.0011946+0.05,  # fraction of total power in return
    pfPipe = 0.0001,          # fraction of total power in piping
    pfNozzle = 0.025,           # fraction of total power in nozzle
    pfReflector = 0.0032335,    # fraction of total power in reflector
    pfPump = 0.0,               # fraction of total power in pump                
    )  # power
 
# correlations-----------------------------------------------------------------
ssInps.AddData("correlations",
    ffFE="Churchill",
    ffMEs="McAdams",
    ffMEr="McAdams",
    ffNoz="Churchill",
    ffRef="Churchill",
    ffPipe="Blasius",
    htFE="WolfMcCarthy",
    htMEs="Taylor",
    htMEr="Taylor",
    htNoz="Taylor",
    htRef="Taylor",
    htPipe="DittusBoelter",
    dpcorrFE=1.0,
    dpcorrMEs=1.0,
    dpcorrMEr=1.0,
    dpcorrNoz=1.0,
    dpcorrRef=1.0,
    dpcorrPipe=1.0,                      
    )  # correlations


# turbopumppump data-----------------------------------------------------------
ssInps.AddData("turbopump",
    noPumps = 1,               # number of pumps in expander cycle
    pumpSs = 4.0,              # target pump specific suction
    pumpPratio = 30.0,         # guess for outlet-to-inlet pump pressure ratio
    pumpMarginPin = 0.1,       # margin inlet pump pressure for caviation
    pumpEfficiency = None,     # pump efficiency
    turbMat = None,            # turbine material
    turbDensity = 2000.0,      # turbine density [kg/m^3]
    turbStress = 1000E+06,     # turbine density [Pa]
    turbType = 'axial',        # turbine type ['axial', 'radial']
    turbD = 0.16,              # turbine rotor diameter [m]
    turbRefPoint = 'inlet',   # ref. point [inlt/outlet] for density calc.
    turbEfficiency = None,     # turbine efficiency
    turbAxTaper = 0.5,         # Blade taper
    turbAxB02B = 4.0,          # inner-to-outer radius of the disk ratio
    # strsF = 0.2,               # stress factor for turbine
    safetyF = 1.2,             # safety factor for turbine
    gearR = None,               # gear ratio turbine-to-pump speed
    cMesh = 50,                # coarse mesh division for finding optimum eff.
    fMesh = 50,                # fine mesh division for finding optimum eff.
    )  # turbopump



# validate that all attributes are defined-------------------------------------
ssInps.ValidateInputs()




