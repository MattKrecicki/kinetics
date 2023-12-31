# -*- coding: utf-8 -*-
"""input containers classes to easy use of ntpTransient code

This file contains the input container class which is used to feed information
into each of the solver classes.

@author: matt krecicki
@email: matthewkrecicki@gmail.edu

"""


import numpy as np
from kinetics.errors.checkerrors import _isbool, _ispositive, _inlist,\
    _ispositiveArray, _isequallength, _iszeropositive, _isnumber, _isstr,\
        _isnonNegativeArray
from kinetics.errors.customerrors import _checkdict, _pkecheck


class pointkineticsInputsContainer:
    """input container can be used for the following solvers:
        -kinetics.functions.pointkinetics.pke
        -kinetics.functions.pointkinetics.srcpke
        -kinetics.functions.tworegionkinetics.pke2region
    
    """
    
    # required keys for each point kinetics solver class
    ALLOWED_TYPS = ["pke2region"]
    
    REQ_KEYS_SRC_PKE = ['beta', 'lamda', 'promptL', 'rhoi', 'S0', 'epsilon',
                        'rhoi', 'volume', 'Q', 'rhoext', 'timepoints']
    
    REQ_KEYS_PKE = ['beta', 'lamda', 'promptL', 'P0', 'volume', 'Q', 'nubar',
                    'rhoext', 'timepoints']
    
    REQ_KEYS_TWO_PKE = ['beta', 'lamda', 'promptLc', 'promptLr', 'fcr', 'frc',
                        'P0', 'volumec', 'volumer', 'Q', 'nubar', 'vc', 'vr',
                        'rhoext', 'timepoints']
    
    #dict containing inputs descriptions and error checking info
    #format is the following:
    #   idx 0 is the first data type
    #   idx 1 is the second data type, which is required if first data type is
    #   a np.ndarray or a list
    #   idx 2 is the input key description
    #   idx 3 is the units of the input
    #   idx 4 is positive value required error checking flag
    
    DICT_TWO_PKE = {"beta": [np.ndarray, float, "delayed neutron group yields",
                             "unitless", True],
                
                    "lamda": [np.ndarray, float, "delayed neutron group decay constants",
                              "1/seconds", True],
                    
                    "promptLc": [float, None, "neutron mean generation time of the "
                                 "core region", "seconds", True],
                
                    "promptLr": [float, None, "neutron mean generation time of the "
                                 "reflector region", "seconds", True],
                
                    "P0": [float, None, "inital reactor power", "watts", True],
                
                    "volumec": [float, None, "volume of the core region", "m3",
                                True],
                
                    "volumer": [float, None, "volume of the reflector region",
                                "m3", True],
                
                    "Q": [float, None, "average recoverable energy per fission",
                          "MeV/fission", True],
                
                    "nubar": [float, None, "average number of neutrons released per "
                              "fission", "n/fission", True],
                    
                    "vc": [float, None, "average core neutron velocity", 
                           "meters/s", True],
                    
                    "vr": [float, None, "average reflector neutron velocity", 
                           "meters/s", True],
                    
                    "rhoext": [object, None, "function that defines excess reactivity as"
                               " a function of time. Should be an instance of the "
                               "general control rules class", False],
                    
                    "timepoints": [np.ndarray, float, "time points for which results "
                                   "are reported by solver", "seconds", False]}
    
    
    def __init__(self, **kwargs):
        """function initialize point kinetics inputs container"""
        
        self.__dict__.update(kwargs)
        self.__checkinputs()
        
    def __getdict(self):
        
        if self.typ == "pke2region":
            Dict = self.DICT_TWO_PKE
        
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
  
        
        