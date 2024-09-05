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
from kinetics.errors.checkerrors import _ispositive, _isnonNegativeArray,\
    _isint


AVERY_SOLVER_DICT = {}


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

    def __solveinitialconditions(self, srcepi=1e-10, keffepi=1e-6, maxitr=100):
        """function solvers inital condition to kind flux and precusor conc."""
        
                
        pass            

    
    # ----- solve transient problem functions

    def __dxdt(self, ):
        pass
    
        
    def __postprocess(self, solution):
        pass
    
    
    def solve(self, srcepi=1e-10, keffepi=1e-6, maxitr=100, transtol=1e-10,
              initialonly=True):
        """function solves transient problem via avery's method"""
        
        # ----- run basic error checking
        _ispositive(srcepi, "relative error tolerance for initial conditions "
                    "source convergence")
        _ispositive(keffepi, "abs error tolerance for keff convergence")      
        _isint(maxitr, "maximum number of power iterations to obtain steady "
               "state neutron flux")    
        _ispositive(transtol, "relative tolerance for tranisent solver")
        
        
        # ----- obtain initial steady-state solution
        x0 = self.__solveinitialconditions(srcepi, keffepi, maxitr)
        
        
        # ----- obtain solution over defined transient
        if initialonly is False:
            solution = odeint(self.__dxdt, x0, self.inputs.timepoints, 
                              rtol=transtol)
        
        # run post-processing features
        #self.__postprocess(solution)
        
        