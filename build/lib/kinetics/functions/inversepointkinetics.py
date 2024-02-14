# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 22:07:53 2024

@author: matt krecicki
@email: matthewkrecicki@gmail.com

function solvers inverse point kinetics problems

"""


MEV_2_J = 1.60218e-13

import numpy as np
from scipy.integrate import quad
from kinetics.errors.checkerrors import _isnegative, _inrange, _ispositive
from kinetics.errors.customerrors import _checkdict, _pkecheck 


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


class inversepke:
    
    
    def __checkinputs(self):
        """function runs basic errpr checking on inverse point kinetics"""
        
        if self.inputs.typ != "invpke":
            raise TypeError("incorrect inputs container given: {}"\
                            .format(self.inputs.typ))
        
        _pkecheck(self.inputs)
        _checkdict(INV_PKE_DICT, self.inputs)
    
    def __init__(self, inputs):
        """function initalizes inverse point kinetics solver"""
        
        self.inputs = inputs
        self.__checkinputs()
        setattr(self.inputs, "groupsDN", len(self.inputs.beta))
        setattr(self.inputs, "betaTot", self.inputs.beta.sum())
        
        self.__checkinputs()
        
    
    def __evaluatederivative(self, factor=1e-3):
        
        dPdt = []
        
        tf = self.inputs.timepoints
          
        
        return np.array(dPdt)
        
    
    
    def __evaluateintegral(self):
        
        pass
    
    
    def __evaluatereactivity(self):
        
        pass
    
    
    def solve(self, rtol=1e-10):
        """function defined kinetics"""
        
        _ispositive(rtol, "relative tolerance of time derivative integrator")
        
        #solve derviate 
        dPdt = self.__evaluatederivative()
        
        #solver integrals
        intPt = self.__evaluateintegral()
        
        #evalute reactivity
        self.__evaluatereactivity()