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
        
        if self.verbose: print(" ... inputs error checking passed")
        

    # ----- matrix construction classes
    
    def __constructpopulationmatrix(self):
        pass
    
    
    def __constructmatrix(self, **kwargs):
        """function constructs the M matrix to predict time derivatives"""
        
        a=1


    # ----- solve initial conditions function

    def __solveinitialconditions(self, srcepi=1e-10, keffepi=1e-6, maxitr=100,
                                 printitr=10):
        """function solvers inital condition to kind flux and precusor conc."""
        
        if self.verbose: print (" ... solving initial conditions")
        
        #generate matrix for power iteration method
        if self.verbose: print(" ... constructing matrix for power iteration")
        initcond = {}
        for dep in self.kineticdata.dependencies:
            initcond[dep] = getattr(self, dep+"0")
        
        mtx0 = self.__constructmatrix(**initcond)
        
        #begin power iteration
        if self.verbose: print(" ... iterating")
        
        keff, l2, itrNo = [], [], 0
                
        while itrNo <= maxitr:
            
            #
            keffi = 1.0
            l2i = 1.0
            
            #print convergence message
            if itrNo % printitr and self.verbose:
                print(" ... iteration {} keff: {:0.5f}, l2: {:0.5e}"\
                      .format(itrNo, keffi, l2i))
            
            itrNo += 1
            
    
    # ----- solve transient problem functions

    def __dxdt(self, ):
        pass
    
        
    def __postprocess(self, solution):
        pass
    
    
    def solve(self, srcepi=1e-10, keffepi=1e-6, maxitr=100, transtol=1e-10,
              initialonly=False):
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
        
        