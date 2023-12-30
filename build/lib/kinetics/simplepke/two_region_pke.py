# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 10:28:29 2023

function solves two-region kinetics model with external source

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""

MEV_2_J = 1.60218e-13

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import newton
from kinetics.errors.checkerrors import _isnegative, _inrange

PKE_W_SRC_KEYS = {}


class twoRegionPKE:
    
    
    def _checkinputs(self):
        """function performs basic error checking"""
        pass
        
    
    def __init__(self, **kwargs):
        """function initializes PKE w/ external source solver class"""
        
        self.__dict__.update(kwargs)
        self._checkinputs()
        self.f = self.fcr * self.frc
        self.betatotal = self.beta.sum()
    
    
    def __solveinitialconditions(self):
        """function determine the inital conidtions for the time integrator"""
        
        #set the number of delayed neutron groups
        self.groupsDN = len(self.beta)
        
        #define initial conditions
        x0 = np.zeros(self.groupsDN+2)
        x0[0] = self.P0
        x0[1] = self.fcr * (self.promptLr / self.promptLc) * self.P0 #initialize reflector population
        x0[2:] = self.P0*self.beta/(self.lamda*self.promptLc) #initialize precusor conc.
        
        return x0
        
    
    def _generatematrix(self, x0, t):
        """function computes matrix"""
        
        # reset the PKE matrix
        mtxA = np.zeros((self.groupsDN+2, self.groupsDN+2))
        
        #compute total reactivity of the system
        rho = self.rho(t)
        
        # build diagonal
        np.fill_diagonal(mtxA, np.append(np.array([0.0, 0.0]), -self.lamda))        
        mtxA[0, 2:] = self.lamda
        
        # build the first column precusor conc.
        mtxA[2:self.groupsDN+2, 0] = self.beta / self.promptLc
        
        # build the first row for core neutron population
        
        nom1 = rho - self.betatotal - self.f * (1 - self.betatotal)
        denom1 = self.promptLc * (1 - self.f)
        alpha1 = nom1 / denom1
        
        nom2 = self.frc * (1 - rho)
        denom2 = self.promptLr * (1 - self.f)
        alpha2 = nom2/ denom2
        
        # build the second row for the reflector neutron population
        xi1 = (self.fcr / self.promptLc) * ((1 - rho)/(1 - self.f))
        xi2 = (1 / self.promptLr) * ((1 - rho)/(1 - self.f))
        
        mtxA[0,0] = alpha1
        mtxA[0,1] = alpha2
        
        mtxA[1,0] = xi1
        mtxA[1,1] = -xi2
        
        return mtxA
        
    
    def _dxdt(self, x0, t):
        """function produces time rate of change for each variable
        
        dX(t)/dt = AX(t)
        
        """
        
        mtx = self._generatematrix(x0, t)
        AX = np.dot(mtx, x0) 
        
        return AX
        
    
    def solve(self):
        """function solves point kinetic equations with external source
        contribution"""
        
        x0 = self.__solveinitialconditions()
        solution = odeint(self._dxdt, x0, self.timepoints, rtol=self.rtol)
        
        #post process results
        setattr(self, "Nc", solution[:,0])
        setattr(self, "Nr", solution[:,1])
        for i in range(self.groupsDN):
            setattr(self, "dg"+str(i+1), solution[:,2+i])
        
        #compute total system reactivity
        totalrho = []
        for t in self.timepoints: totalrho.append(self.rho(t))
        setattr(self, "totalrho", np.array(totalrho))
        
        
        
