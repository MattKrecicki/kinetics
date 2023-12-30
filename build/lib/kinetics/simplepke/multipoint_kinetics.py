# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 11:06:50 2023

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""

MEV_2_J = 1.60218e-13

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import root
from kinetics.errors.checkerrors import _isnegative, _inrange


class averymodel:
    
    
    def __checkinputs(self):
        pass
    
    
    def __init__(self, **kwargs):
        
        self.__dict__.update(kwargs)
        self.__checkinputs()
        
        
    def __initialconditions(self):
        pass

    
    def __generatematrix(self):
        pass
    
    
    def __dxdt(self):
        pass
    
    
    def solve(self, rtol):
        pass