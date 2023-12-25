# -*- coding: utf-8 -*-
"""steadyinputs.py

Inputs containers to store data needed for the steady state and transient
simulations.

This file contains the input container class which is used to feed information
into each of the solver classes.

@author: matt krecicki
@email: matthewkrecicki@gmail.com


"""

import time
import numpy as np

from ntpSystem.errors.checkerrors import _isbool, _inlist,\
    _isnumber, _isstr, _isint, _inrange, _isarray
from ntpSystem.errors.customerrors import ssInputsError

from ntpSystem.functions.turbopumpefficiency import TURBOPUMP_MATS_LIST
from ntpSystem.functions.propertytable import PropertyTable
#    ChangeDataDependencies
from ntpSystem.functions.turbopumpefficiency import TurbopumpEfficiency,\
    TurbopumpMaterials, TurbineGeometricStressFactors, PumpCavitation
import ntpSystem.functions.frictionfactor as fricF
import ntpSystem.functions.heattransfercoefficients as heatT

TURB_MAT_TEMP = '_TEMPORARY_TURB_MAT_'

# -----------------------------------------------------------------------------
# -------- STEADY-STATE KEY DICTONARIES----------------------------------------
# -----------------------------------------------------------------------------

# key represent the attributes
# values are list, where the:
# 1st component is a description of the key
# 2nd component is the default value. None is default does not exist.
# 3rd component is the expected variable type.
# 4th component is the expected value (i.e., range or list) .
   
inf = np.inf

TemplDict = {
"inpfile":
    ["input file for ntpThermo", None, str, None],
"n": 
    ["total number of coolant channels", None, int, [1, inf]],
"chId":
    ["name of the channel Id used by ntpThermo", None, str, None],
"layer0":
   ["index of the very first layer", 0, int, [0, inf]],
"nlayers":
   ["number of layers belonging to this element type", None, int, [1, inf]],    
}

REACTOR  = {}  # empty list for all the components in the reactor
reactorItems = ['fuel', 'reflector', 'nozzle', 'mereturn', 'mesupply', 
                'pipeFuel', 'pipeMod']



# identical fields for all the components
for item in reactorItems:
    for key, val in TemplDict.items():
        REACTOR[item+"_"+key] = val
       
OPERATION = {
"Pt": 
    ["hydrogen storage tank pressure [Pa]", None, float, [0, inf]],
"Tt": 
    ["hydorgen storage tank temperature [K]", None, float, [0, inf]],
"Pc": 
    ["chamber pressure [Pa]", None, float, [0, inf]],
"Tc": 
    ["chamber temperature [K]", None, float, [0, inf]],    
"power": 
    ["total power [watts]", None, float, [0, inf]],
"mdot": 
    ["total mass flow rate [kg/s]", None, float, [0, inf]],
"At": 
    ["nozzle throat area [m^2]", None, float, [0, inf]],
"nozzleExpR": 
    ["nozzle expansion ratio", None, float, [0, inf]], 
} # OPERATION


    
TURBOPUMP = {
"noPumps": 
    ["number of turbopumps in engine system", None, int, [1, inf]],
"pumpPratio": 
    ["pump outlet-to-inlet pressure", None, float, [1, 100]],
"pumpSs": 
    ["target specific suction (20 is max.), unitless", 5.0, float, [0.1, 20]],
"pumpMarginPin": 
    ["margin inlet pump pressure needed for caviation", 0.1, float, [0, 1]],
"pumpEfficiency": 
    ["constant pump efficiency", None, float, [0, 1]],
"pumpD": 
    ["user-defined pump impeller diameter", None, float, [0, inf]],
"pumpQin":
    ["generated heat [Watts] in pump", 0.0, float, [0, inf]],    
"turbType":
    ["turbine type", "axial", str, ["axial", "radial"]],
"turbMat":
    ["material name", None, str, TURBOPUMP_MATS_LIST],
"turbDensity":
    ["material density [kg/m^3]", None, float, [0, inf]], 
"turbStress":
    ["material stress [Pa]", None, float, [0, inf]],
"turbD":
    ["turbine diameter [m]", None, float, [0, inf]],
"turbRefPoint":
    ["turbine inlet/outlet/average reference point for density", 'outlet', str, 
     ['inlet', 'outlet', 'average']],
"turbEfficiency": 
    ["constant turbine efficiency", None, float, [0, 1]],
"turbQin":
    ["generated heat [Watts] in turbine", 0.0, float, [0, inf]], 
"gearR":
    ["gear ratio", 1.0, float, [1, inf]],
"strsF":
    ["stress factor", None, float, [0, inf]],
"turbh2D":
    ["turbine h/D ratio", None, float, [0, inf]],
"safetyF":
    ["safety factor", 1.2, float, [0, inf]],
"turbAxB02B":
    ["inner-to-outer radius of the disk ratio", 4.0, float, [3, 20]],    
"turbAxTaper":
    ["Blade taper (fig.10.8 Emrich's book.)", 0.5, float, [0.4, 1.0]],  
"turbRadr2R":
    ["Blade-to-disc (r/R) radii ratio", 0.4, float, [0.0, 1.0]],  
"cMesh": 
    ["coarse mesh division for finding optimum efficiency", 30, int, [5, 10000]],
"fMesh": 
    ["fine mesh division for finding optimum efficiency", 30, int, [5, 10000]],
}  # TURBOPUMP

          
   
        
#  create a list of required data
COMPONENTS_DICT = {'reactor': REACTOR,
                   'operation': OPERATION, 'turbopump': TURBOPUMP}


COMPONENTS_LIST = list(COMPONENTS_DICT)

# optional attributes ---------------------------------------------------------
OPTIONAL_ATTRS = [
    'power', 'mdot', 'turbMat', 'turbDensity', 'turbStress', 'pumpEfficiency',
    'turbEfficiency', 'gearR', 'turbD', 'turbAxB02B', 'turbAxTaper',
    'turbRadr2R', 'turbh2D', 'strsF', 'pumpD']


class SteadyStateInputs:
    """Object container to store steady state parameters"""
        
    def __init__(self, prntMsgs=True):
        """function initalizes inputs container class
        

        Parameters
        ----------
        prntMsgs : bool, optional
            flag to indicate if messages are printed during
            execution of the code. The default is True.

        Attributes
        ----------
        key : float
            all the keys appearing in:
               FUELGEOM,  PUMPDESIGN

        Returns
        -------
        None.
        
        Example
        -------
        >> ssInps = SteadyStateInputs()

        
        """
        
        self._tic = time.perf_counter()
        _isbool(prntMsgs, "prntMsgs")
        self._prntMsgs = prntMsgs
    
        # Verify that developers define dictionaries properly
        self._ValidateDeveloperDicts()

        if self._prntMsgs:
            print("----------------------------------------------------")
            print("----        Steady State. (0 s)                -----")
            print("----------------------------------------------------")

    def AddData(self, component, **kwargs):
        """add data for a specific component
        
        Example
        -------
        >> SteadyStateInputs.AddData("reactor",
                                     noFuelElems = 127,
                                     ...)
        
        """
        
        _inlist(component, "component", COMPONENTS_LIST)
        attrsList = list(COMPONENTS_DICT[component].keys())
        
        # iterate on all the provided attributes and their values
        for attr, value in kwargs.items():
           value = self.__checkInputs(component, attr, value, attrsList) 
           setattr(self, attr, value)
        
        self._ValidateComponent(component)


    def ValidateInputs(self, H5_PATH=None):
        """check that ALL attributes are defined"""
        
        # loop over all the required components
        for cmpntKey, cmpntDict in COMPONENTS_DICT.items():
            # loop over all the attributes of the specific component
            for attr, values in cmpntDict.items():
                # check if attr is defined
                if not hasattr(self, attr):
                    # print all the allowed attributes
                    _printComponentKeys(cmpntKey)
                    raise ssInputsError("Attr=<{}> for component=<{}> is not "
                                        "defined.".format(attr, cmpntKey))
                    
        # ---------------------------------------------------------------------
        # verify that specific inputs are defined
        # ---------------------------------------------------------------------
        if (self.power is not None) and (self.mdot is not None):
            raise ssInputsError("power and mdot cannot be defined together.\n")

        if (self.turbMat is not None):
            if (self.turbDensity is not None) or (self.turbStress is not None):
                raise ssInputsError(
                    "turbMat and turbDensity / turbStress cannot be defined "
                    "together.")
        else:
            if (self.turbDensity is None) or (self.turbStress is None):
                raise ssInputsError(
                    "Both turbDensity and turbStress must be defined.")
        
        if ((self.pumpEfficiency is not None) and 
            (self.turbEfficiency is None)) or\
           ((self.pumpEfficiency is None) and 
            (self.turbEfficiency is not None)):
                raise ssInputsError(
                    "Both pumpEfficiency and turbEfficiency must be defined.")    

        if (self.gearR is not None) and (self.turbD is not None):
            raise ssInputsError("gearR and turbD cannot be defined together.\n")        

        
        if self._prntMsgs:
            print("... inputs passed error checking.({:.3f} s)"
                  .format(time.perf_counter() - self._tic))
    
        # ---------------------------------------------------------------------
        # Initalize databases
        # ---------------------------------------------------------------------
        # hydrogen properties        
        self.__initialize_hydrogen_props(H5_PATH)
        if self._prntMsgs:
            print("... success in initializing propellant properties."
                  "({:.3f} s)".format(time.perf_counter()-self._tic))
          
        # turbopump properties        
        self.__initialize_turbopump()
        if self._prntMsgs:
            print("... success in reading turbopump properties."
                  "({:.3f} s)".format(time.perf_counter()-self._tic))
            
        # heat transfer and friction factor correlations
        # self.__initialize_correlations()
        # if self._prntMsgs:
        #     print("... success loading heat transfer & friction correlations."
        #           "({:.3f} s)".format(time.perf_counter()-self._tic))


    def __initialize_hydrogen_props(self, H5_PATH):
        """function initalizes hydrogen properties"""
        
        # Read the table with hydrogen properties
        table = PropertyTable(H5_PATH)

        # table is read twice to avoid issues during execution
        H2_PT = table.read("H2")

        # table is read twice to avoid issues during execution
        H2_PH = table.read("H2_PH")


        self.H2_PT = H2_PT
        self.H2_PH = H2_PH
        

    def __initialize_turbopump(self):
        """function initalizes turbopump properties and efficiencies"""
        
        # load pump efficiency 
        self.PUMP_EFF = TurbopumpEfficiency('pump')

        # Pump cavitation curves (specific speed and thoma)
        self.PUMP_CAVIT = PumpCavitation()

        # load turbine efficiency 
        self.TURB_EFF = TurbopumpEfficiency('turbine')
        
        # define a turbopump material object
        self.TURB_MAT = TurbopumpMaterials()
        
        # turbopump stresses
        self.STRS_FCT = TurbineGeometricStressFactors(self.turbType)

        # Add new material if needed
        if self.turbMat is None:
            self.turbMat = TURB_MAT_TEMP
            self.TURB_MAT.AddMaterial(TURB_MAT_TEMP, density=self.turbDensity,
                                     sigy=self.turbStress)  # kg/m^3 & Pa

        if self.strsF is None:
            # Obtain the stress factor for the given turbine
            if self.turbType == 'axial':
                self.strsF, self.turbh2D=\
                    self.STRS_FCT.Value(self.turbAxB02B, self.turbAxTaper)
            else:  # radial turbine
                self.strsF, self.turbh2D =\
                    self.STRS_FCT.Value(r2R=self.turbRadr2R)


    def __initialize_correlations(self):
        """initalizes heat transfer and friction factor correlations"""
        
        self.HT_CORR = {
            self.htFE: getattr(heatT, self.htFE),
            self.htMEs: getattr(heatT, self.htMEs),
            self.htMEr: getattr(heatT, self.htMEr),
            self.htNoz: getattr(heatT, self.htNoz),
            self.htRef: getattr(heatT, self.htRef),
            self.htPipe: getattr(heatT, self.htPipe), }

        self.FF_CORR = {
            self.ffFE: getattr(fricF, self.ffFE),
            self.ffMEs: getattr(fricF, self.ffMEs),
            self.ffMEr: getattr(fricF, self.ffMEr),
            self.ffNoz: getattr(fricF, self.ffNoz),
            self.ffRef: getattr(fricF, self.ffRef),
            self.ffPipe: getattr(fricF, self.ffPipe), }


    def _ValidateDeveloperDicts(self):
        """check that the developers include all the dictionaries properly"""
        
        expValsLen = 4  # total of 4 components should be included
                        # [description, default val, expc var type, expc range]
        msg = 'Developer message: Dictionaries are not defined properly.'
        # loop over all the required components
        for cmpntKey, cmpntDict in COMPONENTS_DICT.items():
            # loop over all the attributes of the specific component
            for attr, values in cmpntDict.items():
                # check if attr is defined
                if len(values) != expValsLen:
                    raise ssInputsError(
                        "{}\nAttr=<{}> for component=<{}> must be of length 4 "
                        "and not {}".format(msg, attr, cmpntKey, len(values)))


    def _ValidateComponent(self, component):
        """check that ALL attributes are defined"""
        
        # loop over all the attributes of the specific component
        for attr, vals in COMPONENTS_DICT[component].items():
            # check if attr is defined
            if not hasattr(self, attr):
                
                # check if there is a default value
                if vals[1] is not None:
                    setattr(self, attr, vals[1])
                elif attr in OPTIONAL_ATTRS:
                    setattr(self, attr, vals[1])
                else:
                    # attribute is not defined and no default value
                    _printAttrNotDefined(component, attr, vals[3])
                    _printComponentKeys(component)
                    raise ssInputsError("Attr=<{}> for component=<{}> is not "
                                        "defined.".format(attr, component))    
    
                    
    def __checkInputs(self, component, attr, value, attrsList):
        """check that inputs are properly provided"""
        
        if attr not in attrsList:
            # print all the allowed attributes for the specific component
            _printComponentKeys(component)
            raise ssInputsError("Attr=<{}> for component=<{}> is not "
                                "defined.".format(attr, component))
            
        if (attr in OPTIONAL_ATTRS) and (value is None):
            return
        
        # expected values for the specific component and attribute
        exp_vals = COMPONENTS_DICT[component][attr]
        
        # check type of the variable according to the pre-defined dictionary
        msg0 = "component <{}> / attr <{}>".format(component, attr)
        try:
            if exp_vals[2] == float:
                _isnumber(value, msg0)
            elif exp_vals[2] == int:
                _isint(value, msg0)        
            elif exp_vals[2] == str:
                _isstr(value, msg0)
            elif (exp_vals[2] == list) or (exp_vals[2] == np.array):
                _isarray(value, msg0)
            else:
                raise TypeError("The variables types defined in SteadyInputs "
                                "container has an error in component {} and "
                                "attr {}".format(component, attr))
        except TypeError as detail:
            _printTypeError(component, attr, value, exp_vals[2])
            raise ssInputsError("{}\n".format(detail))

            
        # check expected value range of the variable
        try:
            if exp_vals[2] == str:
                if exp_vals[3] is not None:
                    _inlist(value, msg0, exp_vals[3])
            elif (exp_vals[2] == list) or (exp_vals[2] == np.array):
                for ival in value:
                    _inrange(ival, msg0, exp_vals[3])
                value = np.array(value)
            else:
                _inrange(value, msg0, exp_vals[3])
        except ValueError as detail:
            _printValueError(component, attr, value, exp_vals[3])
            raise ValueError("{}\n".format(detail))
        except KeyError as detail:
            _printValueError(component, attr, value, exp_vals[3])
            raise ssInputsError("{}\n".format(detail))    
            
        return value

# -----------------------------------------------------------------------------
#                    AUXILIARY FUNCITONS
# -----------------------------------------------------------------------------

def _printComponentKeys(component):
    """print keys and their description for a specific data component"""
    
    print("Component=<{}> \nmust be defined with the following attributes"
          .format(component))
    print("----------------------------------------------------------------\n")
    
    for attr, vals in COMPONENTS_DICT[component].items():
        print("... attr={}\n".format(attr))
        print("         {}".format(vals[0]))
        print("         type={} / expVal={} / defualt={}\n".format(*vals[1:]))
        

def _printValueError(component, attr, value, expvalues):
    """print description if the value is not allowed"""
    
    print("component <{}> / attr=<{}>".format(component, attr))
    print("Value error = {}".format(value)) 
    print("Expected value range = {}".format(expvalues)) 


def _printTypeError(component, attr, value, typevalue):
    """print description if the value is of wrong type"""
    
    print("component <{}> / attr=<{}>".format(component, attr))
    print("For value {} ; type error = {}".format(value, type(value))) 
    print("Expected value type = {}".format(typevalue))  


def _printAttrNotDefined(component, attr, expvalues):
    """print description if the value is of wrong type"""
    
    print("component <{}> / attr=<{}>".format(component, attr))
    print("No value was defined.") 
    print("Expected values range = {}".format(expvalues))  




            

         