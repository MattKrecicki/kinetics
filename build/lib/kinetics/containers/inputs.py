# -*- coding: utf-8 -*-
"""input containers classes to easy use of ntpTransient code

This file contains the input container class which is used to feed information
into each of the solver classes.

@author: matt krecicki
@email: matthewkrecicki@gmail.edu

"""


import numpy as np
from kinetics.errors.checkerrors import _inlist, _isstr, _innotlist, \
    _isequallength
from kinetics.errors.customerrors import _checkdict, _pkecheck


class regionKineticsData:
    
    
    KOBAYASHI_DICT = \
        {"Id": [str, None, "region name", "n/a", False],
         
         "lj": [float, None, "region i's inverse neutron speed",
                "seconds/meters", True],
         
         "kjk": [np.ndarray, None, "region prompt neutron coupling coefficients", "unitless", True],
         
         "kdjk": [np.ndarray, None, "region delayed neutron coupling coefficients", "unitless", True],
         
         "beta": [np.ndarray, float, "delayed neutron fraction", "unitless",
                  True],
         
         "lamda": [np.ndarray, float, "delay neutron group decay constant",
                   "1/seconds", True],
         
         "Q": [float, None, "average recoverable energy released per fission",
               "MeV/fission", True],
         
         "volume": [float, None, "volume of the region", "meters^3", True],
         
         "v": [float, None, "one-group effective neutron velocity",
               "meters/second", True]}
    
    
    AVERY_DICT = \
        {"Ljk": [np.ndarray, float, "region's inverse neutron speed to each region",
                 "seconds/meter", True],
         
         "Kjk": [np.ndarray, float, "region's fission generation coupling coefficient",
                 "unitless", True],
         
         
         "Bjk": [np.ndarray, float, "fraction of delayed neutrons emitted in region j "
                 "that produce fission in region k", "unitless",
                 True],
         
         "Bki": [np.ndarray, float, "delayed neutron fraction for region", "unitless",
                 True],
         
         "lamdaki": [np.ndarray, float, "delay neutron group decay constant",
                     "1/seconds", True],
         
         "Q": [float, None, "average recoverable energy released per fission",
               "MeV/fission", True],
         
         "volume": [float, None, "volume of the region", "meters^3", True],
         
         "v": [float, None, "one-group effective neutron velocity",
               "meters/second", True],
         
         "Id": [str, None, "region name", "n/a", False],}
    
    
    def __checkinputs(self):
        """function runs basic error checking on region kinetics data"""
        
        if "typ" not in list(self.__dict__.keys()):
            raise ValueError("typ keyword must be given for kinetics data container")
            
        _inlist(self.typ, "multi-point kinetic model solution", ["kobayashi", "avery"])
        
        _checkdict(self.__getdict(), self)

    
    def __init__(self, **kwargs):
        """function initalizes instance of region kinetics data container"""
        
        self.__dict__.update(kwargs)
        self.__checkinputs()
        
        
    def __getdict(self):
        """utility function returns input dict"""
            
        if self.typ == "kobayashi":
            Dict = self.KOBAYASHI_DICT
        elif self.typ == "avery":
            Dict = self.AVERY_DICT
            
        return Dict


class multiPointKineticsInputsContainer:
    
    
    def __init__(self, typ):
        """function initalizes multipoint kinetics input container"""
        
        _inlist(typ, "multipoint kinetics solution method", ["Kobayashi", "avery"])
        self.typ = typ
        self.validated = False
        self.regions = []
        self.Ids = []
    
    
    def add(self, **kwargs):
        """function adds a region's kinetic data to container"""
        
        region = regionKineticsData(**kwargs)
        _innotlist(region.Id, "region Id", self.Ids)
        self.Ids.append(region.Id)
        self.regions.append(region)
        
    
    def get(self, Id):
        """function gets a specific region's kinetic data"""
        
        return self.regions[np.where(np.array(self.Ids) == Id)[0][0]]        
    
    
    def validate(self):
        """function ensures that each region can be coupled"""
        
        expL = len(self.Ids)
        
        if self.typ == "avery":
            checkKeys = ["Ljk", "Kjk", "Bjk"]
        else:
            raise ValueError("Kobayashi has not been implemented yet")
        
        for Id in self.Ids:
            region = self.get(Id)
            for key in checkKeys:
                _isequallength(getattr(region, key), expL,
                               region.__getdict()[key][2])
        
        self.validated = True


class pointkineticsInputsContainer:
    """input container can be used for the following solvers:
        -kinetics.functions.pointkinetics.pke
        -kinetics.functions.pointkinetics.srcpke
        -kinetics.functions.tworegionkinetics.pke2region
    
    """
    
    # required keys for each point kinetics solver class
    ALLOWED_TYPS = ["pke2region", "srcpke2region", "pke", "spke"]
    
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
  
        
        