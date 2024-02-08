# -*- coding: utf-8 -*-
"""
Created on Sun Feb  4 18:06:07 2024

function contains multi-point kinetics solver for avery's method

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""

MEV_2_J = 1.60218e-13

import numpy as np
from scipy.integrate import odeint
from kinetics.errors.checkerrors import _ispositive, _isnonNegativeArray, _isint



def poweriteration():
    pass





class avery:
    
    
    def __checkinputs(self):
        pass
    
    
    def __init__(self, **kwargs):
        """function sets up transient solver"""
        self.__dict__.update(**kwargs)
        self.__checkinputs()


    # ----- matrix construction classes
    
    def __constructpopulationmatrix(self):
        pass
    
    
    def __constructdelayedneutronmatrix(self):
        pass


    # ----- solve initial conditions function

    def __solveinitialconditions(self, rtol, k0=1.0, maxitr=None):
        """function solvers inital condition to kind flux and precusor conc."""
        
        _ispositive(k0, "initial guess for keff")
        _isint(maxitr, "maximum number of power iterations to obtain steady "
               "state neutron flux")        
        

    
    # ----- solve transient problem functions

    def __dxdt(self, ):
        pass
    
        
    def __postprocess(self, solution):
        pass
    

    def solve(self, rtol=1e-10, k0=1.0, maxitr=None, noproc=None):
        """function solves transient problem via avery's method"""
        
        _ispositive(rtol, "relative error tolerance for steady-state and"
                    "transient solver")
        
        #obtain initial steady-state solution
        x0 = self.__solveinitialconditions(rtol, k0, maxitr)
        
        # obtain solution over defined transient
        #solution = odeint(self.__dxdt, x0, self.inputs.timepoints, rtol=rtol)
        
        # run post-processing features
        #self.__postprocess(solution)
        
        