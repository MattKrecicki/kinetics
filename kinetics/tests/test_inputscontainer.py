# -*- coding: utf-8 -*-
"""test_inputscontainer

test to ensure inputs container can correctly create input object for
steady-state, operating map, and transient analysis

Created on Tue Jan 10 13:28:42 2023

@author: matt krecicki 

Script to test inputsContainer.py container
"""

import pytest

import numpy as np
from ntpSystem.containers.inputs import inputsContainer


def state_state_input():
    """test to ensure steady state input container is constructed correctly"""
    
    inps = inputsContainer(verbose=False, typ="steadystate")
    
    inps.addFuelGeom(noFuelElems = 127,        #total number of fuel elements in reactor core
                     fuelL = 1.20,             #active core height of fuel elements
                     fuelFlowA = 2.0253e-04,   #flow area of a single fuel element
                     fuelHper = 3.0687e-01)    #heat perimeter of a single fuel element
    
    inps.addModeratorGeom(noModeratorElems = 390,     #total number of moderator elements in reactor core
                          supplyL = 1.20,             #active core height of supply channel
                          supplyFlowA = 1.2566e-05,   #flow area of a single supply channel
                          supplyHper = 1.2566e-02,    #heated perimenter of a single supply channel
                          returnL = 1.20,             #active core height of return channel
                          returnFlowA = 5.8113e-05,   #flow area of a single return channel
                          returnHper = 6.1255e-02)    #heated perimenter of a single return channel
    
    inps.addReflectorGeom(refL = 1.2,          #length of reflector coolant channels
                          noDrums = 18,        #number of control drums in reflector 
                          refNoChans = 100,     #number of coolant channels in reflector 
                          refR = 0.005)     #heated perimenter of a single return channel
    
    inps.addNozzleGeom(nozzleL = 2.65,                #length of nozzle coolant channels
                       nozzleNoChans = 232,           #number of coolant channels in nozzle
                       nozzleR = 0.00257,             #nozzle coolant channel radius
                       nozzleExpR = 300.,             #nozzle expansion ratio
                       At = 0.005662964739541421)     #nozzle throat area, m2
    
    inps.addSysGeom(pipeL = 0.5,      #total piping length, m
                    pipeR = 0.0305)   #piping radius, m)
    
    inps.addTurboPumpPerformance(Pt = 0.22e+6,                         #tank pressure, Pa
                                 noPumps = 1,                          #number of pumps in expander cycle
                                 Tt = 25.0,                            #tank temperature, K
                                 turbPR = 1.7,                         #turbine pressure ratio
                                 turbEff = 0.75,                        #turbine efficiency
                                 Dt = 0.12227655780525866,             #turbine diameter
                                 Dp = 0.028284421281487063,            #pump diameter
                                 nsp = 0.22830463466887901,            #pump specific speed
                                 nst = 0.57,                           #turbine specific speed
                                 dsp = 1.75,                           #pump specific diameter
                                 dst = 2.76,                           #turbine specific diameter
                                 Hdesign = 19533.484743670328,         #pump design head
                                 omegaDesign = 6182.379794473893,      #pump design speed
                                 pumpDesignEff = 0.75,                 #pump efficiency
                                 volFlowDesign = 0.11433223641443514,  #design vol. flow rate
                                 chi = 0.5,                            #fraction of total flow sent to reflector 
                                 H = 10000.0,                          #inital pump head guess, m
                                 tbcvA = 19.0,                         #turbine bypass value coefficient A
                                 tbcvB = 7.0,                          #turbine bypass value coefficient B  
                                 pumpMaterial = "Aluminum",            #pump material
                                 turbineMaterial = "Aluminum")         #turbine material
    
    inps.addPowerFractions(pfFuel = 0.9366,
                           pfSupply = 0.00897,  
                           pfReturn = 0.0512,
                           pfPipe = 0.0,       
                           pfNozzle = 0.025,
                           pfReflector = 0.003233) 
    
    inps.addPressureDropCorrections(uFE=0.25,
                                    uMEs=0.1,
                                    uMEr=0.25,
                                    uNoz=0.1)
    
    inps.addFluidCorrelations(ffFE="Churchill",
                              ffMEs="McAdams",
                              ffMEr="McAdams",
                              ffNoz="Churchill",
                              ffPipe="Blasius",
                              htFE="WolfMcCarthy",
                              htMEs="Taylor",
                              htMEr="Taylor",
                              htNoz="Taylor",
                              htPipe="DittusBoelter")
    
    inps.addChamberConditions(Pc = 6.89e+6,   #desired chamber pressure, Pa
                              Tc = 2700.0) #desired chamber temperature, K 
    
    inps.validate()


def state_state_design_turbopump_input():
    """test to ensure steady state input container is constructed correctly"""
    
    inps = inputsContainer(verbose=False, typ="steadystate")
    
    inps.addFuelGeom(noFuelElems = 127,        #total number of fuel elements in reactor core
                     fuelL = 1.20,             #active core height of fuel elements
                     fuelFlowA = 2.0253e-04,   #flow area of a single fuel element
                     fuelHper = 3.0687e-01)    #heat perimeter of a single fuel element
    
    inps.addModeratorGeom(noModeratorElems = 390,     #total number of moderator elements in reactor core
                          supplyL = 1.20,             #active core height of supply channel
                          supplyFlowA = 1.2566e-05,   #flow area of a single supply channel
                          supplyHper = 1.2566e-02,    #heated perimenter of a single supply channel
                          returnL = 1.20,             #active core height of return channel
                          returnFlowA = 5.8113e-05,   #flow area of a single return channel
                          returnHper = 6.1255e-02)    #heated perimenter of a single return channel
    
    inps.addReflectorGeom(refL = 1.2,          #length of reflector coolant channels
                          noDrums = 18,        #number of control drums in reflector 
                          refNoChans = 100,     #number of coolant channels in reflector 
                          refR = 0.005)     #heated perimenter of a single return channel
    
    inps.addNozzleGeom(nozzleL = 2.65,                #length of nozzle coolant channels
                       nozzleNoChans = 232,           #number of coolant channels in nozzle
                       nozzleR = 0.00257,             #nozzle coolant channel radius
                       nozzleExpR = 300.,             #nozzle expansion ratio
                       At = 0.005662964739541421)     #nozzle throat area, m2
    
    inps.addSysGeom(pipeL = 0.5,      #total piping length, m
                    pipeR = 0.0305)   #piping radius, m)
    
    inps.addTurboPumpDesign(mdotMargin=0.05,
                            maxShaftSpeed=9424.0,
                            specificSuction=5.0,
                            pumpEff=0.75, 
                            turbEff=0.7,
                            turbPR=1.7)  
    
    inps.addPowerFractions(pfFuel = 0.9366,
                           pfSupply = 0.00897,  
                           pfReturn = 0.0512,
                           pfPipe = 0.0,       
                           pfNozzle = 0.025,
                           pfReflector = 0.003233) 
    
    inps.addPressureDropCorrections(uFE=0.25,
                                    uMEs=0.1,
                                    uMEr=0.25,
                                    uNoz=0.1)
    
    inps.addFluidCorrelations(ffFE="Churchill",
                              ffMEs="McAdams",
                              ffMEr="McAdams",
                              ffNoz="Churchill",
                              ffPipe="Blasius",
                              htFE="WolfMcCarthy",
                              htMEs="Taylor",
                              htMEr="Taylor",
                              htNoz="Taylor",
                              htPipe="DittusBoelter")
    
    inps.addChamberConditions(Pc = 6.89e+6,   #desired chamber pressure, Pa
                              Tc = 2700.0) #desired chamber temperature, K 
    
    inps.validate()


def operating_map_input():
    """function tests that steady state operating map input is constructed
    correctly"""
    
    inps = inputsContainer(verbose=False, typ="operatingmap")
    
    inps.addFuelGeom(noFuelElems = 127,        #total number of fuel elements in reactor core
                     fuelL = 1.20,             #active core height of fuel elements
                     fuelFlowA = 2.0253e-04,   #flow area of a single fuel element
                     fuelHper = 3.0687e-01)    #heat perimeter of a single fuel element
    
    inps.addModeratorGeom(noModeratorElems = 390,     #total number of moderator elements in reactor core
                          supplyL = 1.20,             #active core height of supply channel
                          supplyFlowA = 1.2566e-05,   #flow area of a single supply channel
                          supplyHper = 1.2566e-02,    #heated perimenter of a single supply channel
                          returnL = 1.20,             #active core height of return channel
                          returnFlowA = 5.8113e-05,   #flow area of a single return channel
                          returnHper = 6.1255e-02)    #heated perimenter of a single return channel
    
    inps.addReflectorGeom(refL = 1.2,          #length of reflector coolant channels
                          noDrums = 18,        #number of control drums in reflector 
                          refNoChans = 100,     #number of coolant channels in reflector 
                          refR = 0.005)     #heated perimenter of a single return channel
    
    inps.addNozzleGeom(nozzleL = 2.65,                #length of nozzle coolant channels
                       nozzleNoChans = 232,           #number of coolant channels in nozzle
                       nozzleR = 0.00257,             #nozzle coolant channel radius
                       nozzleExpR = 300.,             #nozzle expansion ratio
                       At = 0.005662964739541421)     #nozzle throat area, m2
    
    inps.addSysGeom(pipeL = 0.5,      #total piping length, m
                    pipeR = 0.0305)   #piping radius, m)
    
    inps.addTurboPumpPerformance(Pt = 0.22e+6,                         #tank pressure, Pa
                                 noPumps = 1,                          #number of pumps in expander cycle
                                 Tt = 25.0,                            #tank temperature, K
                                 turbPR = 1.7,                         #turbine pressure ratio
                                 turbEff = 0.75,                        #turbine efficiency
                                 Dt = 0.12227655780525866,             #turbine diameter
                                 Dp = 0.028284421281487063,            #pump diameter
                                 nsp = 0.22830463466887901,            #pump specific speed
                                 nst = 0.57,                           #turbine specific speed
                                 dsp = 1.75,                           #pump specific diameter
                                 dst = 2.76,                           #turbine specific diameter
                                 Hdesign = 19533.484743670328,         #pump design head
                                 omegaDesign = 6182.379794473893,      #pump design speed
                                 pumpDesignEff = 0.75,                 #pump efficiency
                                 volFlowDesign = 0.11433223641443514,  #design vol. flow rate
                                 chi = 0.5,                            #fraction of total flow sent to reflector 
                                 H = 10000.0,                          #inital pump head guess, m
                                 tbcvA = 19.0,                         #turbine bypass value coefficient A
                                 tbcvB = 7.0,                          #turbine bypass value coefficient B  
                                 pumpMaterial = "Aluminum",            #pump material
                                 turbineMaterial = "Aluminum")         #turbine material
    
    inps.addPowerFractions(pfFuel = 0.9366,
                           pfSupply = 0.00897,  
                           pfReturn = 0.0512,
                           pfPipe = 0.0,       
                           pfNozzle = 0.025,
                           pfReflector = 0.003233) 
    
    inps.addPressureDropCorrections(uFE=0.25,
                                    uMEs=0.1,
                                    uMEr=0.25,
                                    uNoz=0.1)
    
    inps.addFluidCorrelations(ffFE="Churchill",
                              ffMEs="McAdams",
                              ffMEr="McAdams",
                              ffNoz="Churchill",
                              ffPipe="Blasius",
                              htFE="WolfMcCarthy",
                              htMEs="Taylor",
                              htMEr="Taylor",
                              htNoz="Taylor",
                              htPipe="DittusBoelter")
    
    inps.addChamberConditions(Pc = np.linspace(2e+6, 7e+6, 10),   #desired chamber pressure, Pa
                              Tc = np.linspace(2000.0, 2700.0, 10)) #desired chamber temperature, K 
    
    inps.validate()


def transient_input():
    """test ensures that transient simulation inputs container is constructed
    correctly"""
    
    inps = inputsContainer(verbose=False, typ="transient")
    
    inps.addFuelGeom(noFuelElems = 127,        #total number of fuel elements in reactor core
                     fuelL = 1.20,             #active core height of fuel elements
                     fuelFlowA = 2.0253e-04,   #flow area of a single fuel element
                     fuelHper = 3.0687e-01)    #heat perimeter of a single fuel element
    
    inps.addModeratorGeom(noModeratorElems = 390,     #total number of moderator elements in reactor core
                          supplyL = 1.20,             #active core height of supply channel
                          supplyFlowA = 1.2566e-05,   #flow area of a single supply channel
                          supplyHper = 1.2566e-02,    #heated perimenter of a single supply channel
                          returnL = 1.20,             #active core height of return channel
                          returnFlowA = 5.8113e-05,   #flow area of a single return channel
                          returnHper = 6.1255e-02)    #heated perimenter of a single return channel
    
    inps.addReflectorGeom(refL = 1.2,          #length of reflector coolant channels
                          noDrums = 18,        #number of control drums in reflector 
                          refNoChans = 100,     #number of coolant channels in reflector 
                          refR = 0.005)     #heated perimenter of a single return channel
    
    inps.addNozzleGeom(nozzleL = 2.65,                #length of nozzle coolant channels
                       nozzleNoChans = 232,           #number of coolant channels in nozzle
                       nozzleR = 0.00257,             #nozzle coolant channel radius
                       nozzleExpR = 300.,             #nozzle expansion ratio
                       At = 0.005662964739541421)     #nozzle throat area, m2
    
    inps.addSysGeom(pipeL = 0.5,      #total piping length, m
                    pipeR = 0.0305)   #piping radius, m)
    
    inps.addTurboPumpPerformance(Pt = 0.22e+6,                         #tank pressure, Pa
                                 noPumps = 1,                          #number of pumps in expander cycle
                                 Tt = 25.0,                            #tank temperature, K
                                 turbPR = 1.7,                         #turbine pressure ratio
                                 turbEff = 0.75,                        #turbine efficiency
                                 Dt = 0.12227655780525866,             #turbine diameter
                                 Dp = 0.028284421281487063,            #pump diameter
                                 nsp = 0.22830463466887901,            #pump specific speed
                                 nst = 0.57,                           #turbine specific speed
                                 dsp = 1.75,                           #pump specific diameter
                                 dst = 2.76,                           #turbine specific diameter
                                 Hdesign = 19533.484743670328,         #pump design head
                                 omegaDesign = 6182.379794473893,      #pump design speed
                                 pumpDesignEff = 0.75,                 #pump efficiency
                                 volFlowDesign = 0.11433223641443514,  #design vol. flow rate
                                 chi = 0.5,                            #fraction of total flow sent to reflector 
                                 H = 10000.0,                          #inital pump head guess, m
                                 tbcvA = 19.0,                         #turbine bypass value coefficient A
                                 tbcvB = 7.0,                          #turbine bypass value coefficient B  
                                 pumpMaterial = "Aluminum",            #pump material
                                 turbineMaterial = "Aluminum")         #turbine material
    
    inps.addPowerFractions(pfFuel = 0.9366,
                           pfSupply = 0.00897,  
                           pfReturn = 0.0512,
                           pfPipe = 0.0,       
                           pfNozzle = 0.025,
                           pfReflector = 0.003233) 
    
    inps.addPressureDropCorrections(uFE=0.25,
                                    uMEs=0.1,
                                    uMEr=0.25,
                                    uNoz=0.1)
    
    inps.addFluidCorrelations(ffFE="Churchill",
                              ffMEs="McAdams",
                              ffMEr="McAdams",
                              ffNoz="Churchill",
                              ffPipe="Blasius",
                              htFE="WolfMcCarthy",
                              htMEs="Taylor",
                              htMEr="Taylor",
                              htNoz="Taylor",
                              htPipe="DittusBoelter")
    
    inps.addKinetics(beta=np.array([2.48080E-04, 1.31503E-03, 1.22769E-03,
                                    2.76846E-03, 1.13800E-03, 4.72114E-04]),
                     Lambda=np.array([1.33372E-02, 3.27303E-02,  1.20796E-01,
                                      3.02995E-01, 8.50277E-01,  2.85558E+00]),
                     genTime=9.35059E-05)
    
    inps.addIsotopics(U235=6.4323e+20, U238=2.57292e+21, I135=0.0, Xe135=0.0,
                      Pm149=0.0, Sm149=0.0)
    
    inps.addMaterialProperties(massFuel=300.0, cpFuel=300.0, massMod=500.0,
                               cpModerator=300.0, fuelVol=1000.0) 
    
    #define transient
    
    #inital conditions
    Tc0 = 1500.0
    Pc0 = 0.195*6.89e+6

    #end of thrust buildup 1
    Tc1 = 2000.0
    Pc1 = 0.35*6.89e+6

    dt1 = 15.0
    
    inps.defineTransient(Tci=Tc0,
                         Pci=Pc0,
                         Tcf=Tc1,
                         Pcf=Pc1,
                         dtRamp=dt1,
                         tmax=dt1+5.0,
                         noIters=200,
                         theta0=120.0)
    
    coeff1 = (Pc1 - Pc0) / dt1
     
    inps.addPressureControlRule(["linear"],           #define type of control rule
                                [[0.0, coeff1]], #coefficients that define control rule, (A*x + B)
                                np.array([dt1]))

    coeff3 = (Tc1- Tc0) / dt1

    inps.addTemperatureControlRule(["linear"],           #define type of control rule
                                   [[0.0, coeff3]], #coefficients that define control rule, (A*x + B)
                                   np.array([dt1]))
    
    inps.addDataBase(pumpHead = "pumpHead.pkl")
    
    inps.validate()


def check_negative_value():
    
    with pytest.raises(ValueError,
                       match="fuel element geometry input noFuelElems must be zero or positive and not -127"):
    
        inps = inputsContainer(verbose=False, typ="steadystate")
    
        inps.addFuelGeom(noFuelElems = -127,       
                     fuelL = 1.20,            
                     fuelFlowA = 2.0253e-04,
                     fuelHper = 3.0687e-01)   


def check_missing_value():
    
    with pytest.raises(KeyError,
                       match='fuelHper: heated perimeter of a single fuel element in units of meters, was not defined'):
    
        inps = inputsContainer(verbose=False, typ="steadystate")
    
        inps.addFuelGeom(noFuelElems = 127,     
                         fuelL = 1.20,           
                         fuelFlowA = 2.0253e-04) 

