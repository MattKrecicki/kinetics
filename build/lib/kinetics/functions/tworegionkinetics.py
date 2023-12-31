# -*- coding: utf-8 -*-
"""tworegionkinetics.py

Created on Sun Dec 31 06:48:25 2023

function contains two region kinetics solvers

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""

MEV_2_J = 1.60218e-13

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import root
from kinetics.errors.checkerrors import _ispositive
from kinetics.containers.outputs import pointkineticsOutputsContainer


# -----------------------------------------------------------------------------
# ----- two region kinetics solver class 
# -----------------------------------------------------------------------------


class pke2region:
    
    
    def __init__(self, inputs):
        """function initalizes two region kinetics solver class"""
        
        self.inputs = inputs
        self.__checkinputs()
        setattr(self.inputs, "groupsDN", len(self.inputs.beta))
        setattr(self.inputs, "betaTot", self.inputs.beta.sum())
        setattr(self.inputs, "f", self.inputs.fcr * self.inputs.frc)
    
    
    def __checkinputs(self):
        """function does basic error checking"""
        
        #make sure correct input container is given to solve
        if self.inputs.typ != "pke2region":
            raise TypeError("incorrect inputs container given: {}"\
                            .format(self.inputs.typ))
    
    
    def __solveinitialconditions(self):
        """function finds initial conditions for two region kinetics"""
        
        x0 = np.zeros(self.inputs.groupsDN+2)
        
        #calculate total core neutron population
        self.inputs.n0 = (self.inputs.P0*self.inputs.nubar*self.inputs.promptLc)/(self.inputs.Q*MEV_2_J)
        
        x0[0] = self.inputs.n0
        #initialize reflector population
        x0[1] = self.inputs.fcr*(self.inputs.promptLr/self.inputs.promptLc)*self.inputs.n0
        #initialize precusor conc.
        x0[2:] = self.inputs.n0*self.inputs.beta/(self.inputs.lamda*self.inputs.promptLc)
        
        return x0
    
    def __constructmatrix(self, t):
        """function constructs matrix for estimation of time derv."""
        
        # reset the PKE matrix
        mtxA = np.zeros((self.inputs.groupsDN+2, self.inputs.groupsDN+2))
        
        #compute total reactivity of the system
        rho = self.inputs.rhoext.evaluate(t)
        
        # build diagonal
        np.fill_diagonal(mtxA, np.append(np.array([0.0, 0.0]), -self.inputs.lamda))        
        mtxA[0, 2:] = self.inputs.lamda
        
        # build the first column precusor conc.
        mtxA[2:self.inputs.groupsDN+2, 0] = self.inputs.beta / self.inputs.promptLc
        
        # build the first row for core neutron population
        alpha1 = (rho - self.inputs.betaTot - self.inputs.f * (1 - self.inputs.betaTot)) / (self.inputs.promptLc * (1 - self.inputs.f))
        alpha2 = (self.inputs.frc * (1 - rho)) / (self.inputs.promptLr * (1 - self.inputs.f))
        
        # build the second row for the reflector neutron population
        xi1 = (self.inputs.fcr / self.inputs.promptLc) * ((1 - rho)/(1 - self.inputs.f))
        xi2 = (1 / self.inputs.promptLr) * ((1 - rho)/(1 - self.inputs.f))
        
        mtxA[0,0] = alpha1
        mtxA[0,1] = alpha2
        
        mtxA[1,0] = xi1
        mtxA[1,1] = -xi2
        
        return mtxA
    
    
    def __dxdt(self, x0, t):
        """function evaluates the time derv. of the system"""
        
        return np.dot(self.__constructmatrix(t), x0)
    
    
    def __postprocess(self, solution):
        """function performs post processing routines"""
        
        #seperate out delayed neutron precusors
        dnt = solution[:,2:].T
        
        #seperate out total neutron populations
        ntc = solution[:,0]
        ntr = solution[:,1]
        
        #convert neutron populations to fluxes
        fluxc = (self.inputs.vc * ntc) / self.inputs.volumec
        fluxc /= 10000 # convert n/m2-s to n/cm2-s
        
        fluxr = (self.inputs.vc * ntr) / self.inputs.volumer
        fluxr /= 10000 # convert n/m2-s to n/cm2-s
        
        #convert core neutron population into power
        power = (ntc * self.inputs.Q * MEV_2_J) / (self.inputs.nubar * self.inputs.promptLc)
        
        #obtain time dependent external reactivity
        rho = []
        for t in self.inputs.timepoints: rho.append(self.inputs.rhoext.evaluate(t))
        rho = np.array(rho)
        
        #initalize solution container
        self.outputs = \
            pointkineticsOutputsContainer(timepoints=self.inputs.timepoints,\
                ntc=ntc, ntr=ntr, dnt=dnt, power=power, rho=rho, fluxc=fluxc,
                fluxr=fluxr, typ="pke2region")
        
    
    def solve(self, rtol=1e-10):
        """function solver two-region kinetics problem"""
        
        _ispositive(rtol, "relative tolerance of time derivative integrator")
        
        # obtain initial conditions
        x0 = self.__solveinitialconditions()
        
        # obtain solution over defined transient
        solution = odeint(self.__dxdt, x0, self.inputs.timepoints, rtol=rtol)
        
        # run post-processing features
        self.__postprocess(solution)
    
    