# -*- coding: utf-8 -*-
"""pointkinetics.py

Created on Thu Dec 28 12:32:15 2023

function contains point kinetics and source driven point kinetics solver

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""

MEV_2_J = 1.60218e-13

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import root
from kinetics.errors.checkerrors import _isnegative, _inrange, _ispositive
from kinetics.containers.outputs import pointkineticsOutputsContainer


# -----------------------------------------------------------------------------
# ----- point kinetics solver class 
# -----------------------------------------------------------------------------


class pke:
    
    
    def __init__(self, inputs):
        """function initalizes point kinetics solver class"""
        
        self.inputs = inputs
        self.__checkinputs()
        setattr(self.inputs, "groupsDN", len(self.inputs.beta))
        setattr(self.inputs, "betaTot", self.inputs.beta.sum())
    
    
    def __checkinputs(self):
        """function does basic error checking"""
        
        #make sure correct input container is given to solve
        if self.inputs.typ != "pke":
            raise TypeError("incorrect inputs container given: {}"\
                            .format(self.inputs.typ))
    
    
    def __solveinitialconditions(self):
        """function solves initial conditions"""
        
        #convert power to total neutron population
        self.inputs.n0 = \
            (self.inputs.P0 * self.inputs.nubar * self.inputs.promptL) / \
            (self.inputs.Q * MEV_2_J)
        
        #detemerine inital reactor power and precusor concentrations
        x0 = \
            self.inputs.n0*np.append(1, self.inputs.beta/self.inputs.lamda/\
                                     self.inputs.promptL)
        
        return x0
        
    
    def __constructmatrix(self, t):
        """function constructs matrix for estimation of time derv."""
        
        #construct matrix
        mtxA = np.zeros((self.inputs.groupsDN+1, self.inputs.groupsDN+1))
        
        #fill diagonal with decay constants
        np.fill_diagonal(mtxA, np.append(0, -self.inputs.lamda))
        
        mtxA[0,:] = \
            np.append((self.inputs.rhoext.evaluate(t)-self.inputs.betaTot)/\
                      self.inputs.promptL, self.inputs.lamda)
        
        # build the first column
        mtxA[1:self.inputs.groupsDN+1, 0] = self.inputs.beta/self.inputs.promptL
        
        return mtxA
    
    
    def __dxdt(self, x0, t):
        """function evaluates the time derv. of the system"""
        
        return np.dot(self.__constructmatrix(t), x0)
    
    
    def __postprocess(self, solution):
        """function performs post processing routines"""
        
        #seperate out delayed neutron precusors
        dnt = solution[:,1:].T
        
        #seperate out total neutron population
        nt = solution[:,0]
        
        #convert neutron population to flux
        flux = (self.inputs.v * nt) / self.inputs.volume
        flux /= 10000 # convert n/m2-s to n/cm2-s
        
        #convert neutron population into power
        power = (nt*self.inputs.Q*MEV_2_J)/(self.inputs.nubar*self.inputs.promptL)
        
        #obtain time dependent external reactivity
        rho = []
        for t in self.inputs.timepoints: rho.append(self.inputs.rhoext.evaluate(t))
        rho = np.array(rho)
        
        #initalize solution container
        self.outputs = \
            pointkineticsOutputsContainer(timepoints=self.inputs.timepoints,\
                nt=nt, dnt=dnt, power=power, rho=rho, flux=flux, typ="pke")
    
    
    def solve(self, rtol=1e-10):
        """function defined kinetics"""
        
        _ispositive(rtol, "relative tolerance of time derivative integrator")
        
        # obtain initial conditions
        x0 = self.__solveinitialconditions()
        
        # obtain solution over defined transient
        solution = odeint(self.__dxdt, x0, self.inputs.timepoints, rtol=rtol)
        
        # run post-processing features
        self.__postprocess(solution)


# -----------------------------------------------------------------------------
# ----- source driven point kinetics solver class 
# -----------------------------------------------------------------------------


class srcpke:
    
    
    def __init__(self, inputs):
        """function initalizes source driven point kinetics solver class"""
        
        self.inputs = inputs
        self.__checkinputs()
        setattr(self.inputs, "groupsDN", len(self.inputs.beta))
        setattr(self.inputs, "betaTot", self.inputs.beta.sum())

    
    def __checkinputs(self):
        """function does basic error checking"""
        
        #make sure correct input container is given to solve
        if self.inputs.typ != "spke":
            raise TypeError("incorrect inputs container given: {}"\
                            .format(self.inputs.typ))
    
    
    def __findinitalpopulation(self, ni):
        """function calculates inital reactor power using newton-raphson method"""
        
        delayed = np.sum(self.inputs.beta/(self.inputs.promptL))*ni
        
        prompt = ((self.inputs.rhoi-self.inputs.betaTot)/self.inputs.promptL)*ni
        
        return prompt + delayed + self.inputs.srcnt
    
    
    def __solveinitialconditions(self, rtol):
        """function solves initial conditions"""
        
        #compute neutron production rate by source
        self.inputs.srcnt = \
            self.inputs.nubar*self.inputs.S0*self.inputs.epsilon
        
        #numerically compute inital neutron population
        convrg = root(self.__findinitalpopulation, 1.0, tol=rtol,
                      method="broyden1", options={"maxiter": 1000})
        self.inputs.n0 = convrg.x
        
        #if solution could not be found, fail run
        if convrg.success is False:
            raise ValueError("inital neutron population solution did not "
                             "converge!")
        
        #initialize inital conditions
        x0 = \
            self.inputs.n0*np.append(1, \
                self.inputs.beta/self.inputs.lamda/self.inputs.promptL)
        
        return x0    
    
    
    def __constructmatrix(self, t):
        """function constructs matrix for estimation of time derv."""
        
        #construct matrix
        mtxA = np.zeros((self.inputs.groupsDN+1, self.inputs.groupsDN+1))
        
        #fill diagonal with decay constants
        np.fill_diagonal(mtxA, np.append(0, -self.inputs.lamda))
        
        mtxA[0,:] = \
            np.append((self.inputs.rhoext.evaluate(t)+\
                self.inputs.rhoi-self.inputs.betaTot)/self.inputs.promptL,
                self.inputs.lamda)
        
        # build the first column
        mtxA[1:self.inputs.groupsDN+1, 0] = self.inputs.beta/self.inputs.promptL
        
        return mtxA
    
    
    def __dxdt(self, x0, t):
        """function evaluates the time derv. of the system"""
        
        dxdt = np.dot(self.__constructmatrix(t), x0)
        dxdt[0] += self.inputs.srcnt
        
        return dxdt
    
    
    def __postprocess(self, solution):
        """function performs post processing routines"""
        
        #seperate out delayed neutron precusors
        dnt = solution[:,1:].T
        
        #seperate out total neutron population
        nt = solution[:,0]
        
        #convert neutron population to flux
        flux = (self.inputs.v * nt) / self.inputs.volume
        flux /= 10000 # convert n/m2-s to n/cm2-s
        
        #convert neutron population into power
        power = \
            (nt*self.inputs.Q*MEV_2_J)/(self.inputs.nubar*self.inputs.promptL)
        
        #obtain time dependent external reactivity
        rho = []
        for t in self.inputs.timepoints: rho.append(self.inputs.rhoext.evaluate(t))
        rho = np.array(rho) + self.inputs.rhoi
        
        #initalize solution container
        self.outputs = \
            pointkineticsOutputsContainer(timepoints=self.inputs.timepoints, \
                nt=nt, dnt=dnt, power=power, rho=rho, flux=flux, typ="spke")
    
    
    def solve(self, rtol=1e-10):
        """function defined kinetics"""
        
        _ispositive(rtol, "relative tolerance of time derivative integrator")
        
        # obtain initial conditions
        x0 = self.__solveinitialconditions(rtol)
        
        # obtain solution over defined transient
        solution = odeint(self.__dxdt, x0, self.inputs.timepoints, rtol=rtol)
        
        # run post-processing features
        self.__postprocess(solution)
    
    