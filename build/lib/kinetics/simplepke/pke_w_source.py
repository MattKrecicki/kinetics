# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 10:28:29 2023

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""

MEV_2_J = 1.60218e-13

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import newton
from kinetics.errors.checkerrors import _isnegative, _inrange


PKE_W_SRC_KEYS = {}


class PKEwSource:
    
    
    def _checkinputs(self):
        """function performs basic error checking"""
        
        _isnegative(self.rhoi, "inital negative reactivity in system")
        _inrange(self.epsilon, "external source efficiency", [0.0, 1.0],
                 upBound=True, lowBound=True)
    
    
    def __init__(self, **kwargs):
        """function initializes PKE w/ external source solver class"""
        
        self.__dict__.update(kwargs)
        self._checkinputs()
    
    
    def __numericalfindinitalpower(self, Pi):
        """function calculates inital reactor power using newton-raphson method"""

        delayed = np.sum(self.beta / (self.promptL) ) * Pi
        
        prompt = ((self.rhoi - self.beta.sum())/self.promptL) * Pi
        
        balance = prompt + delayed + self.srcPower
        
        return balance
    
    
    def __solveinitialconditions(self):
        """function determine the inital conidtions for the time integrator"""
        
        #convert MeV/fission to J/fission
        self.Qjoules = self.Q * MEV_2_J
        
        #set the number of delayed neutron groups
        self.groupsDN = len(self.beta)
        
        #compute power produced by source
        self.srcPower = (self.S0 * self.epsilon * self.Qjoules) / self.promptL
        
        #numerically obtain inital power
        self.P0, convrg = newton(self.__numericalfindinitalpower, 1.0,
                         maxiter=500, tol=1e-10, full_output=True)
        if convrg.converged is False:
            raise ValueError("numerical solution of intial power did not "
                             "converge, contact developer matt krecicki")
        
        #initialize precusor conc.
        self.C0 = self.P0*self.beta/(self.lamda*self.promptL)
        
        x0 = np.zeros(self.groupsDN+1)
        x0[0] = self.P0
        x0[1:] = self.C0
        
        return x0
        
    
    def _generatematrix(self, x0, t):
        """function computes matrix"""
        
        # reset the PKE matrix
        mtxA = np.zeros((self.groupsDN+1, self.groupsDN+1))
        
        # build diagonal
        np.fill_diagonal(mtxA, np.append(0, - self.lamda))

        # build the first row
        rho = self.rhoi + self.rho(t)
        
        mtxA[0, :] = np.append((rho-self.beta.sum())/self.promptL, self.lamda)
                
        # build the first column
        mtxA[1:self.groupsDN+1, 0] = self.beta / self.promptL
        
        return mtxA
        
    
    def _dxdt(self, x0, t):
        """function produces time rate of change for each variable in vec(x)"""
        
        mtx = self._generatematrix(x0, t)
        
        dndt = np.dot(mtx, x0)
        dndt[0] += self.srcPower
        
        return dndt
        
    
    def solve(self):
        """function solves point kinetic equations with external source
        contribution"""
        
        x0 = self.__solveinitialconditions()
        
        solution = odeint(self._dxdt, x0, self.timepoints)
        
        setattr(self, "power", solution[:,0])
        setattr(self, "srcpower", solution[:,-1])
        for i in range(self.groupsDN):
            setattr(self, "dg"+str(i+1), solution[:,1+i])
        
        #compute total system reactivity
        totalrho = []
        for t in self.timepoints: totalrho.append(self.rho(t) + self.rhoi)
        setattr(self, "totalrho", np.array(totalrho))
        
        
