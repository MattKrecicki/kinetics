# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 10:28:29 2023

function solves two-region kinetics model with external source and accounts
for photoneutron production

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""

MEV_2_J = 1.60218e-13

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import newton
from kinetics.errors.checkerrors import _isnegative, _inrange

PKE_W_SRC_KEYS = {}


class twoRegionPKEwSourcewPhotoN:
    
    
    def _checkinputs(self):
        """function performs basic error checking"""
        
        _isnegative(self.rhoi, "inital negative reactivity in system")
        _inrange(self.epsilon, "external source efficiency", [0.0, 1.0],
                 upBound=True, lowBound=True)
        
    
    def __init__(self, **kwargs):
        """function initializes PKE w/ external source solver class"""
        
        self.__dict__.update(kwargs)
        self._checkinputs()
    
    
    def __findinitalpower(self, Pi):
        """function calculates inital reactor power using newton-raphson method"""
        
        delayed = np.sum(self.beta / (self.promptLc) ) * Pi
        
        alpha1 = ((self.rhoi - self.beta.sum() - self.f*(1 - self.beta.sum())) / (1-self.f))/self.promptLc
        alpha2 = ((1 - self.rhoi)/(1 - self.f)) * (self.frc/self.promptLr)
        promptcore = alpha1 * Pi + alpha2 * self.fcr * (self.promptLr / self.promptLc) * Pi
        
        xi1 = (self.fcr / self.promptLc) * ((1 - self.rhoi)/(1 - self.f))
        xi2 = (1 / self.promptLr) * ((1 - self.rhoi)/(1 - self.f))
        promptreflector = xi1 * Pi - xi2 * self.fcr * (self.promptLr / self.promptLc) * Pi
        
        balance = promptcore + promptreflector + delayed + self.srcPower
        
        return balance
    
    
    def __solveinitialconditions(self):
        """function determine the inital conidtions for the time integrator"""
        
        #convert MeV/fission to J/fission
        self.Qjoules = self.Q * MEV_2_J
        
        #set the number of delayed neutron groups
        self.groupsDN = len(self.beta)
        
        #compute power produced by source
        self.srcPower = (self.S0 * self.epsilon * self.Qjoules) / self.promptLc
        
        #numerically obtain inital power
        P0, convrg = newton(self.__findinitalpower, 1.0, tol=1e-15,
                                 maxiter=500, rtol=1e-30, full_output=True)
        #if solution could not be found, fail run
        if convrg.converged is False:
            raise ValueError("numerical solution of intial power did not "
                             "converge")
        
        #define initial conditions
        x0 = np.zeros(self.groupsDN+2)
        x0[0] = P0
        x0[1] = self.fcr * (self.promptLr / self.promptLc) * P0 #initialize reflector population
        x0[2:] = P0*self.beta/(self.lamda*self.promptLc) #initialize precusor conc.
        
        return x0
        
    
    def _generatematrix(self, x0, t):
        """function computes matrix"""
        
        # reset the PKE matrix
        mtxA = np.zeros((self.groupsDN+2, self.groupsDN+2))
        
        #compute total reactivity of the system
        rho = self.rhoi + self.rho(t)
        
        # build diagonal
        np.fill_diagonal(mtxA, np.append(np.array([0.0, 0.0]), -self.lamda))
            
        mtxA[0, 2:] = self.lamda
        
        # build the first column precusor conc.
        mtxA[2:self.groupsDN+2, 0] = self.beta / self.promptLc
        
        # build the first row for core neutron population
        alpha1 = ((rho - self.beta.sum() - self.f*(1 - self.beta.sum())) / (1-self.f))/self.promptLc
        alpha2 = ((1 - rho)/(1 - self.f)) * (self.frc/self.promptLr)
        
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
        
        dX(t)/dt = AX(t) + B
        
        """
        
        mtx = self._generatematrix(x0, t)
        
        AX = np.dot(mtx, x0) 
        
        AX[0] += self.srcPower
        
        return AX
        
    
    def solve(self):
        """function solves point kinetic equations with external source
        contribution"""
        
        x0 = self.__solveinitialconditions()
        
        solution = odeint(self._dxdt, x0, self.timepoints, rtol=1e-30)
        
        setattr(self, "Nc", solution[:,0])
        setattr(self, "Nr", solution[:,1])
        for i in range(self.groupsDN):
            setattr(self, "dg"+str(i+1), solution[:,2+i])
        
        #compute total system reactivity
        totalrho = []
        for t in self.timepoints: totalrho.append(self.rho(t) + self.rhoi)
        setattr(self, "totalrho", np.array(totalrho))
        
        
        
