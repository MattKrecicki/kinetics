# -*- coding: utf-8 -*-
"""input containers classes to easy use of ntpTransient code

This file contains the input container class which is used to feed information
into each of the solver classes.

@author: matt krecicki
@email: matthewkrecicki@gmail.edu

"""


import numpy as np
from kinetics.errors.checkerrors import _inlist, _isstr, _innotlist, \
    _isequallength, _isint, _islist, _ispositive, _iszeropositive
from kinetics.errors.customerrors import _checkdict, _pkecheck
from kinetics.functions.interpolation import LinearInterp, BiLinearInterp, \
    TriLinearInterp


class averyregion:
    
    AVERY_DICT = \
        {"Ljk": [np.ndarray, float, "region's prompt neutron mean generation time to each region",
                 "seconds/meter", True],
         
         "Ljkd": [np.ndarray, float, "region's delayed neutron mean generate time to each region",
                  "seconds/meter", True],
         
         "Kjk": [np.ndarray, float, "region's fission generation coupling coefficient",
                 "unitless", True],
         
         "Kjkd": [np.ndarray, float, "fraction of delayed neutrons emitted in region j "
                 "that produce fission in region k", "unitless",
                 True],
         
         "Bki": [np.ndarray, float, "delayed neutron fraction for region", "unitless",
                 True],
         
         "lamdaki": [np.ndarray, float, "delay neutron group decay constant",
                     "1/seconds", True],
         
         "Q": [float, None, "average recoverable energy released per fission",
               "MeV/fission", True],
         
         "v": [float, None, "one-group effective neutron velocity",
               "meters/second", True],
         
         "coupling": [list, str, "list the details coupling of each region",
                      "n/a", True],
         }
    
    
    def _checkinputs(self):
        """function runs error checking on state point of region inputs"""
        
        #run basic errror checking on inputs 
        _checkdict(self.AVERY_DICT, self)
        
        #make sure all the arrays are the correct length
        expL = len(self.coupling)
        Lkeys = ["Ljk", "Ljkd", "Kjk", "Kjkd"]
        for key in Lkeys:
            _isequallength(getattr(self, key), expL, self.AVERY_DICT[key][2])
        
        #check values that should be zero positive
        zposKeys = ["Ljk", "Ljkd", "Kjk", "Kjkd", "Bki", "lamdaki"]
        for key in zposKeys:
            for val in getattr(self, key):
                _iszeropositive(val, "{} in {}".format(val, key))
        
        #check values that must be positive
        _ispositive(self.Q, "average recoverable energy released per fission")
        _ispositive(self.v, "one-group effective neutron velocity")
        
    
    def __init__(self, **kwargs):
        """function initalizes avery's multipoint kinetics region container"""
        
        self.__dict__.update(kwargs)
        self._checkinputs()


class regionKineticsData:
    
    
    REGION_DICT = {"dependencies": [list, str, 
                                    "list of dependences for input data", "n/a",
                                    True],
                   
                   "Id": [str, None, "region name", "n/a", True],
                   
                   "typ": [str, None, "multi-point kinetics model type", "n/a",
                           True],
                   
                   "x": [float, None, "x-axis position", "meters", False],
                   
                   "y": [float, None, "y-axis position", "meters", False],
                   
                   "z": [float, None, "z-axis position", "meters", False],
                   
                   "volume": [float, None, "volume of the region", "meters^3", 
                              True]}
    
    
    def __checkinputs(self):
        """function runs basic error checking on region kinetics data"""
        
        
        #make sure all required keys are given
        
        if "typ" not in list(self.__dict__.keys()):
            raise ValueError("typ keyword must be given for kinetics data container")
        
        _inlist(self.typ, "multi-point kinetic model solution", ["avery"])
        
        

    
    def __init__(self, **kwargs):
        """function initalizes instance of region kinetics data container"""
        
        self.__dict__.update(kwargs)
        self.__checkinputs()
        
        self.states = []
        #initalize dependency lists
        for dep in self.dependencies:
            setattr(self, dep, [])
        
        self.valid = False
        
    
    def _getdict(self):
        """utility function returns input dict"""
        
        if self.typ == "avery":
            Dict = self.AVERY_DICT
            
        return Dict

    
    
    def add(self, **kwargs):
        """function allows user to add a state to the region container"""
        
        if self.typ == "avery":
            
            state = averyregion(**kwargs)
            
            for dep in self.dependencies:
                val = getattr(state, dep) 
                deplist = getattr(self, dep)
                deplist.append(val)
                setattr(self, dep, deplist)
            
            self.states.append(state)
        
        else:
            raise ValueError("only avery MPK method has been implemented")

    
    
    def __validate(self):
                
        self.valid = True
    

class multiPointKineticsInputsContainer:
    
    
    def _checkinputs(self):
        """fucntion checks inputs for class initialization"""
        
        _inlist(self.typ, "multipoint kinetics solution method", ["avery"])
        _isint(self.nregions, "number of regions in multipoint model")
        _islist(self.dependencies, "list of variables that the inputs are dependent on")
    
    
    def __init__(self, typ=False, nregions=False, dependencies=None, order=False):
        """function initalizes multipoint kinetics input container"""       
        
        self.typ = typ
        self.nregions = nregions
        self.order = order
        self.dependencies = {}
        
        #setup dependencies dict
        for di in dependencies:
            self.dependencies[di] = []        
        self._checkinputs()
        
        self.regions = []
        self.Ids = []
        self.validated = False
        
    
    
    def add(self, region):
        """function adds a region's kinetic data to container"""
                
        _innotlist(region.Id, "region Id", self.Ids)
        if self.region.typ != self.typ:
            raise ValueError("contain data typ is {}, region {}'s data typ is "
                             "{}".format(self.typ, region.Id, region.typ))
        
        self.Ids.append(region.Id)
        self.regions.append(region)
        
    
    def get(self, Id):
        """function gets a specific region's kinetic data"""
        
        return self.regions[np.where(np.array(self.Ids) == Id)[0][0]]        
    
    
    def validate(self):
        """function ensures that each region can be coupled"""
        
        #make sure all of the data is the right size
        
        #make all regions are defined
        
        
        
        expL = len(self.Ids)
        
        self.validated = True


    


class pointkineticsInputsContainer:
    """input container can be used for the following solvers:
        -kinetics.functions.pointkinetics.pke
        -kinetics.functions.pointkinetics.srcpke
        -kinetics.functions.tworegionkinetics.pke2region
    
    """
    
    # required keys for each point kinetics solver class
    ALLOWED_TYPS = ["pke2region", "srcpke2region", "pke", "spke", "invpke",
                    "srcinvpke"]
    
    #dict containing inputs descriptions and error checking info
    #format is the following:
    #   idx 0 is the first data type
    #   idx 1 is the second data type, which is required if first data type is
    #   a np.ndarray or a list
    #   idx 2 is the input key description
    #   idx 3 is the units of the input
    #   idx 4 is positive value required error checking flag
    
    DICT_PKE = \
        {"beta": [np.ndarray, float, "delayed neutron group yields", "unitless",
                  True],
                
         "lamda": [np.ndarray, float, "delayed neutron group decay constants",
                   "1/seconds", True],
                    
         "promptL": [float, None, "neutron mean generation time of the core"
                      " region", "seconds", True],
                
         "P0": [float, None, "inital reactor power", "watts", True],
         
         "volume": [float, None, "volume of the core region", "m3", True],
         
         "Q": [float, None, "average recoverable energy per fission",
               "MeV/fission", True],
                
         "nubar": [float, None, "average number of neutrons released per "
                   "fission", "n/fission", True],
                    
         "v": [float, None, "average core neutron velocity", "meters/s", True],
                    
         "rhoext": [object, None, "function that defines excess reactivity as"
                    " a function of time. Should be an instance of the general"
                    " control rules class", False],
                    
         "timepoints": [np.ndarray, float, "time points for which results "
                        "are reported by solver", "seconds", False]}
        
    
    DICT_SRC_PKE = \
        {"beta": [np.ndarray, float, "delayed neutron group yields", "unitless",
                  True],
                
         "lamda": [np.ndarray, float, "delayed neutron group decay constants",
                   "1/seconds", True],
                    
         "promptL": [float, None, "neutron mean generation time of the core"
                      " region", "seconds", True],
                
         "S0": [float, None, "source neutron generation rate", "n/s", True],
         
         "epsilon": [float, None, "Source efficiency", "fissions/source n",
                     True],
         
         "rhoi": [float, None, "initial excess reactivity", "dk/k", False],
         
         "volume": [float, None, "volume of the core region", "m3", True],
         
         "Q": [float, None, "average recoverable energy per fission",
               "MeV/fission", True],
                
         "nubar": [float, None, "average number of neutrons released per "
                   "fission", "n/fission", True],
                    
         "v": [float, None, "average core neutron velocity", "meters/s", True],
                    
         "rhoext": [object, None, "function that defines excess reactivity as"
                    " a function of time. Should be an instance of the general"
                    " control rules class", False],
                    
         "timepoints": [np.ndarray, float, "time points for which results "
                        "are reported by solver", "seconds", False]}
    
    
    DICT_TWO_PKE = \
        {"beta": [np.ndarray, float, "delayed neutron group yields", "unitless",
                  True],
                
         "lamda": [np.ndarray, float, "delayed neutron group decay constants",
                   "1/seconds", True],
                    
         "promptLc": [float, None, "neutron mean generation time of the core"
                      " region", "seconds", True],
                
         "promptLr": [float, None, "neutron mean generation time of the "
                      "reflector region", "seconds", True],
                
         "P0": [float, None, "inital reactor power", "watts", True],
                
         "volumec": [float, None, "volume of the core region", "m3", True],
                
         "volumer": [float, None, "volume of the reflector region", "m3", True],
                
         "Q": [float, None, "average recoverable energy per fission",
               "MeV/fission", True],
                
         "nubar": [float, None, "average number of neutrons released per "
                   "fission", "n/fission", True],
                    
         "vc": [float, None, "average core neutron velocity", "meters/s", True],
                    
         "vr": [float, None, "average reflector neutron velocity", "meters/s",
                True],
                    
         "rhoext": [object, None, "function that defines excess reactivity as"
                    " a function of time. Should be an instance of the general"
                    " control rules class", False],
                    
         "timepoints": [np.ndarray, float, "time points for which results "
                        "are reported by solver", "seconds", False]}
    
    
    DICT_TWO_SRC_PKE = \
        {"beta": [np.ndarray, float, "delayed neutron group yields", "unitless",
                  True],
                
         "lamda": [np.ndarray, float, "delayed neutron group decay constants",
                   "1/seconds", True],
                    
         "promptLc": [float, None, "neutron mean generation time of the core"
                      " region", "seconds", True],
                
         "promptLr": [float, None, "neutron mean generation time of the "
                      "reflector region", "seconds", True],
                
         "S0": [float, None, "neutron source strength", "n/s", True],
         
         "epsilon": [float, None, "Source efficiency", "fissions/source n",
                     True],
         
         "rhoi": [float, None, "initial excess reactivity", "dk/k", False],
                
         "volumec": [float, None, "volume of the core region", "m3", True],
                
         "volumer": [float, None, "volume of the reflector region", "m3", True],
                
         "Q": [float, None, "average recoverable energy per fission",
               "MeV/fission", True],
                
         "nubar": [float, None, "average number of neutrons released per "
                   "fission", "n/fission", True],
                    
         "vc": [float, None, "average core neutron velocity", "meters/s", True],
                    
         "vr": [float, None, "average reflector neutron velocity", "meters/s",
                True],
                    
         "rhoext": [object, None, "function that defines excess reactivity as"
                    " a function of time. Should be an instance of the general"
                    " control rules class", False],
                    
         "timepoints": [np.ndarray, float, "time points for which results "
                        "are reported by solver", "seconds", False]}
    
    INV_PKE_DICT = \
        {"beta": [np.ndarray, float, "delayed neutron fraction", "unitless", True],
         
         "lamda": [np.ndarray, float, "delay neutron group decay constant",
                   "1/seconds", True],
         
         "promptL": [float, None, "prompt neutron lifeime", "seconds", True],
              
         "volume": [float, None, "total volume of reactor control volume",
                    "meters^3", True],
         
         "nubar": [float, None, "average number of neutrons produced per fission",
                   "neutrons/fission", True],
         
         "Q": [float, None, "Average recoverable energy released per fission",
               "MeV/fission", True],
         
         "v": [float, None, "effective one-group neutron velocity", "meters/second",
               True],
         
         "timepoints": [np.ndarray, float, "time points to return solution",
                        "seconds", False],
         
         "power": [object, None, "reactor power control class", "n/a",
                   False],
         
         "typ": [str, None, "type of kinetic simulation desired", "n/a", False]}
        
        
    SRC_INV_PKE_DICT = \
        {"beta": [np.ndarray, float, "delayed neutron fraction", "unitless", True],
         
         "lamda": [np.ndarray, float, "delay neutron group decay constant",
                   "1/seconds", True],
         
         "promptL": [float, None, "prompt neutron lifeime", "seconds", True],
              
         "volume": [float, None, "total volume of reactor control volume",
                    "meters^3", True],
         
         "nubar": [float, None, "average number of neutrons produced per fission",
                   "neutrons/fission", True],
         
         "Q": [float, None, "Average recoverable energy released per fission",
               "MeV/fission", True],
         
         "v": [float, None, "effective one-group neutron velocity", "meters/second",
               True],
         
         "timepoints": [np.ndarray, float, "time points to return solution",
                        "seconds", False],
         
         "power": [object, None, "reactor power control class", "n/a",
                   False],
         
         "S0": [float, None, "initial neutron source strength", "neutrons/second",
                True],
         
         "rho0": [float, None, "initial reactivity", "dk/k", False],
         
         "epsilon": [float, None, "neutron source efficiency",
                     "fissions/neutron emitted", True],
         
         "typ": [str, None, "type of kinetic simulation desired", "n/a", False]}

    
    def __init__(self, **kwargs):
        """function initialize point kinetics inputs container"""
        
        self.__dict__.update(kwargs)
        self.__checkinputs()
        
        
    def __getdict(self):
        """utility function returns input dict"""
        
        if self.typ == "pke2region":
            Dict = self.DICT_TWO_PKE
        elif self.typ == "srcpke2region":
            Dict = self.DICT_TWO_SRC_PKE
        elif self.typ == "pke":
            Dict = self.DICT_PKE
        elif self.typ == "spke":
            Dict = self.DICT_SRC_PKE
        elif self.typ == "invpke":
            Dict = self.INV_PKE_DICT
        elif self.typ == "srcinvpke":
            Dict = self.SRC_INV_PKE_DICT
        
        return Dict
    
    
    def __checkinputs(self):
        """function runs error checking to make sure inputs correctly defined"""
        
        #make sure input container typ is defined
        if "typ" not in list(self.__dict__.keys()):
            raise ValueError("typ must be defined")
        _inlist(self.typ, "input container type", self.ALLOWED_TYPS)
        
        #get dict required for error checking
        Dict = self.__getdict()
        
        #run input dict error checking
        _checkdict(Dict, self)
        
        #run extra pke error check
        _pkecheck(self)

    
    def whatis(self, key):
        
        #make sure input key is a str
        _isstr(key, "input key {}".format(key))
        
        #check if key is allowed for input container typ
        Dict = self.__getdict()
        _inlist(key, "requested input key {}".format(key), list(Dict.keys()))
        
        #return key definition
        print("input key: {}".format(key))
        print("    Description: {}".format(Dict[key][2]))
        print("    units:       {}".format(Dict[key][3]))
  
        
        