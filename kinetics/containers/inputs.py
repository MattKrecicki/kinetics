# -*- coding: utf-8 -*-
"""input containers classes to easy use of ntpTransient code

This file contains the input container class which is used to feed information
into each of the solver classes.

@author: matt krecicki
@email: matthewkrecicki@gmail.edu

"""

H5_PATH = "ThermoPhysicalProperties.h5"

import pickle 
import numpy as np
from copy import deepcopy
from ntpSystem.errors.checkerrors import _isbool, _ispositive, _inlist,\
    _ispositiveArray, _isequallength, _iszeropositive, _isnumber, _isstr,\
        _isnonNegativeArray
from ntpSystem.functions.propertytable import PropertyTable, ChangeDataDependencies
from ntpSystem import setDataPath
from ntpSystem.functions.frictionfactor import FF_CORRELATIONS
from ntpSystem.functions.heattransfercoefficients import HT_CORRELATIONS
import ntpSystem.functions.frictionfactor as fricF
import ntpSystem.functions.heattransfercoefficients as heatT
from ntpSystem.functions.control import generalControlRule



class inputsContainer:
    """class is input container for operating map script"""
    
    SIMULATION_TYPS = ["operatingmap", "steadystate", "transient"]
    
    #STEADY-STATE KEY DICTONARIES----------------------------------------------
    
    FUELGEOM_ADD_DICT = {"noFuelElems": "total number of fuel elements in the reactor core",
                         "fuelL": "axial length of a single fuel elements, in units of meters",
                         "fuelFlowA": "flow area of a single fuel element in units of sq. meters",
                         "fuelHper": "heated perimeter of a single fuel element in units of meters"}
    
    POWFRAC_ADD_DICT = {"pfFuel": "fraction of total reactor power deposited in the fuel elements",
                        "pfSupply": "fraction of total reactor power desposited in the moderator supply channels",
                        "pfReturn": "fraction of total reactor power despoited in the moderator return channels", 
                        "pfPipe": "fraction of total reactor power deposited in the piping of the system",
                        "pfReflector": "fraction of total reactor power deposited in the radial reflector and control drums",
                        "pfNozzle": "fraction of total reactor power picked up in regeneratively cooled nozzle"}   
    
    MODGEOM_ADD_DICT = {"noModeratorElems": "total number of moderator elements in the reactor core",
                        "supplyL": "axial length of a single moderator element supply channel in units of meters",
                        "supplyFlowA": "flow area of a single moderator element supply channel in units of sq. meters",
                        "supplyHper": "heated perimeter of a single moderator element supply channel in units of meters",
                        "returnL": "axial length of a single moderator element return channel in units of meters",
                        "returnFlowA": "flow area of a single moderator element return channel in units of sq. meters",
                        "returnHper": "heated perimeter of a single moderator element return channel in units of meters"}
    
    REFGEOM_ADD_DICT = {"refL": "axial length of radial reflector in units of meters",
                        "noDrums": "total number of contorl drums in radial reflector",
                        "refNoChans": "number of coolant channels in each control drum",
                        "refR": "radius of each coolant channel in the control drums in units of meters"}
    
    NOZZLE_ADD_DICT = {"nozzleL": "lenght of the regnerative cooled jacket in units of meters",
                       "nozzleNoChans": "number of coolant channels in regneratively cooled nozzle",
                       "nozzleR": "radius of each coolant channel in the regneratively cooled nozzle in units of meters",
                       "At": "nozzle throat area in units of sq. meters",
                       "nozzleExpR": "nozzle expansion ratio"}
    
    SYS_ADD_DICT = {"pipeL": "pipe length in units of meters",
                    "pipeR": "pipe radius in units of meters"}
    
    TURBO_ADD_DICT = {"Pt": "hydrogen storage tank pressure in units of pascals",
                      "Tt": "hydorgen storage tank temperature in units of kelvin",
                      "noPumps": "number of turbopumps in engine system",
                      "turbPR": "turbine pressure ratio, unitless",
                      "turbEff": "turbine efficiency, unitless",
                      "chi": "fraction of total system mass-flow sent to the reflector circuit, unitless",
                      "H": "inital guess for pump head, in units of meters",
                      "tbcvA": "bypass value alpha coefficient, unitless",
                      "tbcvB": "bypass value beta coefficient, unitless",
                      "Dt": "diameter of turbine in units of meters",
                      "Dp": "diameter of pump in units of meters",
                      "Hdesign": "head of pump at design point, in units of meters",
                      "omegaDesign": "shaft speed of pump at design point, in units of rads/s",
                      "pumpDesignEff": "pump efficiency at design point",
                      "pumpMaterial": "material of the pump",
                      "turbineMaterial": "material of turbine",
                      "volFlowDesign": "pump volumetric flow at design point, in units of m3/s",
                      "nsp": "specific speed of pump, unitless",
                      "nst": "specific speed of turbine, unitless",
                      "dsp": "specific diameter of pump, unitless",
                      "dst": "specific diameter of turbine, unitless"}
    
    CHAMBER_DICT = {"Tc": "chamber temperature in units of kelvin",
                    "Pc": "chamber pressure in units of pascals"}
    
    DP_CORR_DICT = {"uFE": "fuel element  channel pressure drop correction fudge factor",
                    "uMEs": "moderator element supply channel pressure drop correction fudge factor",
                    "uMEr": "moderator element return channel pressure drop correction fudge factor",
                    "uNoz": "nozzle channel pressure drop correction fudge factor",}
    
    SOL_CORR_DICT = {"ffFE": "fuel element channel friction factor correlation",
                     "ffMEs": "moderator element supply channel friction factor correlation",
                     "ffMEr": "moderator element return channel friction factor correlation",
                     "ffNoz": "nozzle channel friction factor correlation",
                     "ffPipe": "pipe friction factor correlation",
                     "htFE": "fuel element channel heat transfer correlation",
                     "htMEs": "moderator element supply channel heat transfer correlation",
                     "htMEr": "moderator element return channel heat transfer correlation",
                     "htNoz": "nozzle channel heat transfer correlation",
                     "htPipe": "pipe heat transfer correlation"}
        
    #TURBOPUMP DESIGN KEY DICTONARIES------------------------------------------
    
    TURBO_DESIGN_DICT = {"mdotMargin": "fraction of total mass-flow sent through bypass line at design point, in units of kg/s",
                         "maxShaftSpeed": "maxiumum allowable shaft speed at design point, in units of rads/s",
                         "specificSuction": "specific suction at design point",
                         "pumpEff": "desired pump efficiency",
                         "turbEff": "desired turbine efficiency",
                         "turbPR": "turbine pressure ratio"}
    
    #TRANSIENT KEY DICTONARIES-------------------------------------------------
            
    KINETICS_KEYS_DICT = {"beta": "array of delayed neutron fractions",
                          "Lambda": "array of delay neutron half-lifes",
                          "genTime": "mean generation time, in units of seconds"}
    
    ISOTOPICS_KEYS_DICT = {"U235": "atom concentration of U-235 in fuel elements, in units of atoms/cm3",
                           "U238": "atom concentration of U-238 in fuel elements, in units of atoms/cm3",
                           "I135": "atom concentration of I-135 in fuel elements, in units of atoms/cm3",
                           "Xe135": "atom concentration of Xe-135 in fuel elements, in units of atoms/cm3",
                           "Pm149": "atom concentration of Pm-149 in fuel elements, in units of atoms/cm3",
                           "Sm149": "atom concentration of Sm-149 in fuel elements, in units of atoms/cm3"}
    
    MAT_PROPS_KEYS_DICT = {"cpFuel": "specific heat capacity of fuel material, in units of J/(kg K)",
                           "cpModerator": "specific heat capacity of the moderator material, in units of J/(kg K)",
                           "massFuel": "total mass of the fuel material, in units of kg",
                           "massMod": "total mass of the moderator material, in units of kg",
                           "fuelVol": "total volume of fuel elements in reactor core, in units of m3"}
    
    TRANS_DEF_DICT = {"Tci": "inital chamber temperature in units of kelvin",
                      "Pci": "inital chamber pressure in units of pascals",
                      "Tcf": "final chamber temperature in units of kelvin",
                      "Pcf": "final chamber pressure in units of pascals",
                      "dtRamp": "total time length of transient startup in units of seconds",
                      "tmax": "maximum time of transient solution, in units of seconds",
                      "noIters": "number time points results are returned on",
                      "theta0": "inital control drum position in units of degrees"}
    
    DATABASE_DICT = {"pumpHead": "python pickle object that contains pump head as a function of chamber temperature and pressure"}
    
    REQUIRED_STEADYSTATE_KEYS = []
    
    REQUIRED_STEADYSTATE_KEYS.extend(list(FUELGEOM_ADD_DICT.keys()))
    REQUIRED_STEADYSTATE_KEYS.extend(list(POWFRAC_ADD_DICT.keys()))
    REQUIRED_STEADYSTATE_KEYS.extend(list(MODGEOM_ADD_DICT.keys()))
    REQUIRED_STEADYSTATE_KEYS.extend(list(REFGEOM_ADD_DICT.keys()))
    REQUIRED_STEADYSTATE_KEYS.extend(list(SYS_ADD_DICT.keys()))
    REQUIRED_STEADYSTATE_KEYS.extend(list(TURBO_ADD_DICT.keys()))
    REQUIRED_STEADYSTATE_KEYS.extend(list(CHAMBER_DICT.keys()))
    REQUIRED_STEADYSTATE_KEYS.extend(list(NOZZLE_ADD_DICT.keys()))
    REQUIRED_STEADYSTATE_KEYS.extend(list(DP_CORR_DICT.keys()))
    
    OPTIONAL_STEADYSTATE_KEYS = deepcopy(REQUIRED_STEADYSTATE_KEYS)
    OPTIONAL_STEADYSTATE_KEYS.extend(list(TURBO_DESIGN_DICT.keys()))
    for key in list(TURBO_ADD_DICT.keys()):
        OPTIONAL_STEADYSTATE_KEYS.remove(key)
    
    REQUIRED_TRANSIENT_KEYS = deepcopy(REQUIRED_STEADYSTATE_KEYS)    
    REQUIRED_TRANSIENT_KEYS.extend(list(KINETICS_KEYS_DICT.keys()))
    REQUIRED_TRANSIENT_KEYS.extend(list(ISOTOPICS_KEYS_DICT.keys()))
    REQUIRED_TRANSIENT_KEYS.extend(list(MAT_PROPS_KEYS_DICT.keys()))
    REQUIRED_TRANSIENT_KEYS.extend(list(TRANS_DEF_DICT.keys()))
    REQUIRED_TRANSIENT_KEYS.append("temperatureControlRule")
    REQUIRED_TRANSIENT_KEYS.append("pressureControlRule")
    for key in list(CHAMBER_DICT.keys()):
        REQUIRED_TRANSIENT_KEYS.remove(key)
    
    
    def __init__(self, verbose=True, typ="steadystate"):
        """function initalizes inputs container class
        

        Parameters
        ----------
        verbose : bool, optional
            flag to set verbosity, i.e. if messages are printed during
            execution of the code. The default is True.
        typ : str, optional
            defines what kind of simulation is being run. Each simulation type
            requires different error checking and pre-processing The default
            is "steadystate".

        Returns
        -------
        None.
        
        Example
        -------
        >> inputsContainer(verbose=True, typ="steadystate")

        
        """
        
        _isbool(verbose, "verbose flag")
        _inlist(typ, "simulation type", self.SIMULATION_TYPS)
        self.verbose = verbose
        self.simulationType = typ
        self._valid = False #flag that solver classes use to make sure inputs are okay
        self.designTurboPump = False #flag used to have steady-state solver design turbo-pump
    
    
    def __printKeyDictionary(self, checkKeys, checkDict):
        """helper function that checks all keys are defined"""
        
        reqKeys = list(checkDict.keys())
        for key in reqKeys:
            if key not in checkKeys:                
                raise KeyError("{}: {}, was not defined"\
                                   .format(key, checkDict[key]))
                    
    
    def __checkFloatInputs(self, message, checkDict, checkKeys):
        """helper function to check inputs that are all floats"""
        
        for key in checkKeys:
            _inlist(key, message + " {}".format(key),
                    list(checkDict.keys()))
            _iszeropositive(getattr(self, key),
                        "fuel element geometry input {}".format(key))
        
        self.__printKeyDictionary(checkKeys, checkDict)
        
    
    def addGeom(self, **kwargs):
        """function adds system piping geometry to inputs container
        
        Example
        -------
        >> inputsContainer.addSysGeom(pipeL = 2.5, 
                                      pipeR = 0.0305)
        
        """
        
        #add kwargs to dictionary
        self.__dict__.update(kwargs)
        #check definitions
        self.__checkFloatInputs("geometry input", 
                                self.SYS_ADD_DICT,
                                list(kwargs.keys()))
    
    
    def addChamberConditions(self, **kwargs):
        """function adds chamber conditions to inputs container. Function can be
        used for both steady-state and operating map solver classes. 

        
        Examples
        --------
        
        #for steady state-solver class
        inputsContainer.addChamberConditions(Pc = 6.89e+6, Tc = 2700.0)
        
        #for operating map solver class
        inputsContainer.addChamberConditions(Pc=np.linspace(1e+6, 7e+6, 20),
                                             Tc=np.linspace(900.0, 3000.0, 20))

        """
        
        if self.simulationType == "transient":
            #nothing needs to be done here, chamber conditions are determined
            #by control laws
            pass
        
        elif self.simulationType == "steadystate":
            #add kwargs to dictionary
            self.__dict__.update(kwargs)
            #check definitions
            self.__checkFloatInputs("chamber conditions input", 
                                    self.CHAMBER_DICT,
                                    list(kwargs.keys()))
        
        elif self.simulationType == "operatingmap":
            #add kwargs to dictionary
            self.__dict__.update(kwargs)
            #preform error checking
            for key in list(kwargs.keys()):
                _inlist(key, "chamber condition input {}".format(key),
                        list(self.CHAMBER_DICT.keys()))
                _ispositiveArray(getattr(self, key),
                             "chamber condition input {}".format(key))
            #make sure are equal length
            keys = list(kwargs.keys())
            _isequallength(getattr(self, keys[0]), len(getattr(self, keys[1])),
                           "chamber conditions")
    
    
    def addFluidCorrelations(self, **kwargs):
        """function adds fluid solver friction factor and heat transfer
        correlations to inputs container
        
        example
        -------
        
        >> inputsContainer.addFluidCorrelation(ffFE="Churchill",
                                               ffMEs="McAdams",
                                               ffMEr="McAdams",
                                               ffNoz="Churchill",
                                               ffPipe="Blasius",
                                               htFE="WolfMcCarthy",
                                               htMEs="Taylor",
                                               htMEr="Taylor",
                                               htNoz="Taylor",
                                               htPipe="DittusBoelter")
        
        """
        
        #add kwargs to dictionary
        self.__dict__.update(kwargs)
        
        #preform error checking
        for key in list(kwargs.keys()):
            _inlist(key, "solution correlation input", self.SOL_CORR_DICT)
            if "ff" in key:
                _inlist(getattr(self, key), "friction factor correlation", FF_CORRELATIONS)
            elif "ht" in key:
                _inlist(getattr(self, key), "heat transfer correlation", HT_CORRELATIONS)
        
        self.__printKeyDictionary(list(kwargs.keys()), self.SOL_CORR_DICT)

    
    def addIsotopics(self, **kwargs):
        """function adds inital isotopic concentration to inputs container
        
        example
        -------
        >> inputsContainer.addIsotpics(U235=6.4323e+20, U238=2.57292e+21,
                                       I135=0.0, Xe135=0.0, Pm149=0.0,
                                       Sm149=0.0)
        
        """
        #add kwargs to dictionary
        self.__dict__.update(kwargs)
        #check definitions
        self.__checkFloatInputs("isotopic concentration inputs", 
                                self.ISOTOPICS_KEYS_DICT,
                                list(kwargs.keys()))


    def addKinetics(self, **kwargs):
        """function adds kinetics data to inputs container
        
        example
        -------
        >> inputsContainer.addKinetics(beta=np.array([2.48080E-04, 1.31503E-03,
                                                      1.22769E-03, 2.76846E-03,
                                                      1.13800E-03, 4.72114E-04]),
                                       Lambda=np.array([1.33372E-02, 3.27303E-02,
                                                        1.20796E-01, 3.02995E-01,
                                                        8.50277E-01,  2.85558E+00]),
                                       genTime=9.35059E-05)
        
        """
        #add kwargs to dictionary
        self.__dict__.update(kwargs)
        #run error checking
        self.__printKeyDictionary(list(kwargs.keys()), self.KINETICS_KEYS_DICT)
        _isnonNegativeArray(getattr(self, "beta"), "delayed neutron fractions")
        _isnonNegativeArray(getattr(self, "Lambda"),
                            "delayed neutron precusor life-lives")
        _ispositive(getattr(self, "genTime"), "neutron mean generation time")   
        
        
    def defineTransient(self, **kwargs):
        """function adds transient definition to file
        
        example
        -------
        >> inputsContainer.defineTransient(Tci=1200.0, Pci=2e+6, Tcf=2700.0,
                                           Pcf=6.89e+6, dtRamp=30.0, tmax=35.0,
                                           noIters=60, theta0=120.0)
        
        """
        #add kwargs to dictionary
        self.__dict__.update(kwargs)
        #check inputs
        self.__checkFloatInputs("transient definition inputs",
                                self.TRANS_DEF_DICT,
                                list(kwargs.keys()))
    
    
    def addTemperatureControlRule(self, funcTyps, coeffs, tends):
        """function adds temperature control rule to inputs container
        

        Parameters
        ----------
        funcTyps : list
            list of control law function types.
        coeffs : list
            list of control law function coefficients.
        tends : np.1darray
            list of time bounds for control laws.

        Returns
        -------
        None.
        
        example
        -------
        >> inputsContainer.addTemperatureControlRule(["linear"],
                                                     [[0.0, 10.0]],
                                                     np.array([15.0]))

        """
        self.temperatureControlRule =\
                generalControlRule(funcTyps, coeffs, tends)
    
    
    def addPressureControlRule(self, funcTyps, coeffs, tends):
        """function adds pressure control rule to inputs container
        

        Parameters
        ----------
        funcTyps : list
            list of control law function types.
        coeffs : list
            list of control law function coefficients.
        tends : np.1darray
            list of time bounds for control laws.

        Returns
        -------
        None.
        
        example
        -------
        >> inputsContainer.addPressureControlRule(["linear"],
                                                  [[0.0, 0.5e+6]],
                                                  np.array([5.0]))

        """
        self.pressureControlRule =\
               generalControlRule(funcTyps, coeffs, tends)

    
    def validate(self):
        """function preforms error checking to ensure all required keys have
        been defined. Function also initalizes hydrogen properties and fluid
        solution correlations"""
        
        #make sure all required keys have been initalized
        if self.simulationType == "steadystate":
            
            if self.designTurboPump is True:
                keys = self.OPTIONAL_STEADYSTATE_KEYS
            else:
                keys = self.REQUIRED_STEADYSTATE_KEYS
                
            for key in keys:
                _inlist(key, "input data {}".format(key),
                        list(self.__dict__.keys()))
        
        elif self.simulationType == "operatingmap":
            
            if self.designTurboPump is True:
                raise KeyError("turbopump design cannot be performed during \
                               operating map analysis")
            
            for key in self.REQUIRED_STEADYSTATE_KEYS:
                _inlist(key, "input data {}".format(key),
                        list(self.__dict__.keys()))
        
        elif self.simulationType == "transient":
            for key in self.REQUIRED_TRANSIENT_KEYS:
                _inlist(key, "input data {}".format(key),
                        list(self.__dict__.keys()))
        
        if self.verbose: print("inputs passed error checking...")
                
        if self.verbose: print("initalizing heat transfer and friction factor correlations...")
        self.__initalizeFluidCorrelations()
        
        if self.verbose: print("initalizing H2 properties...")
        dataPath = setDataPath(H5_PATH)
        self.__initalize_hydrogen_props(dataPath)
            
        self._valid = True
   
        
        