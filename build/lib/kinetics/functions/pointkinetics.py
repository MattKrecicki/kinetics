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
from kinetics.containers.outputs import pointkineticscontainer


# -----------------------------------------------------------------------------
# ----- point kinetics solver class 
# -----------------------------------------------------------------------------


KEYS_PKE = {"beta": [list, float, "delayed neutron group yields, unitless"],
            
            "lamda": [list, float, "delayed neutron group decay constants, in"
                      " units of 1/seconds"],
            
            "promptL": [float, None, "neutron mean generation time, in units"
                        " of seconds"],
            
            "P0": [float, None, "inital reactor power, in units of watts"],
            
            "volume": [float, None, "volume of fissioning region, in units of"
                       "m3"],
            
            "Q": [float, None, "average recoverable energy per fission, in "
                  "units of MeV/fission"],
            
            "nubar": [float, None, "average number of neutrons released per "
                      "fission, in units of n/fission"],
            
            "rhoext": ["class", None, "function that defines excess reactivity as"
                       " a function of time. Should be an instance of the "
                       "general control rules class"],
            
            "timepoints": ["np array", float, "time points for which results "
                           "are reported by solver"],
            
            "rtol": [float, None, "relative error convergence criteria for "
                     "inital conditions and time derivative evaluation"]}


REQ_KEYS_PKE = ['beta', 'lamda', 'promptL', 'P0', 'volume', 'Q', 'nubar',
                'rhoext', 'timepoints']


class pke:
    
    
    def __init__(self, **kwargs):
        """function initalizes point kinetics solver class"""
        
        self.__dict__.update(kwargs)
        self.__checkinputs()
        self.groupsDN = len(self.beta)
        self.betaTot = self.beta.sum()
    
    
    def __checkinputs(self):
        """function does basic error checking"""
        pass
    
    
    def __solveinitialconditions(self):
        """function solves initial conditions"""
        
        #convert power to total neutron population
        self.n0 = (self.P0 * self.nubar * self.promptL) / (self.Q * MEV_2_J)
        
        #detemerine inital reactor power and precusor concentrations
        x0 = self.n0*np.append(1, self.beta/self.lamda/self.promptL)
        
        return x0
        
    
    def __constructmatrix(self, t):
        """function constructs matrix for estimation of time derv."""
        
        #construct matrix
        mtxA = np.zeros((self.groupsDN+1, self.groupsDN+1))
        
        #fill diagonal with decay constants
        np.fill_diagonal(mtxA, np.append(0, -self.lamda))
        
        mtxA[0,:] = \
            np.append((self.rhoext.evaluate(t)-self.betaTot)/self.promptL, self.lamda)
        
        # build the first column
        mtxA[1:self.groupsDN+1, 0] = self.beta / self.promptL
        
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
        flux = (self.v * nt) / self.volume
        flux /= 10000 # convert n/m2-s to n/cm2-s
        
        #convert neutron population into power
        power = (nt * self.Q * MEV_2_J) / (self.nubar * self.promptL)
        
        #obtain time dependent external reactivity
        rho = []
        for t in self.timepoints: rho.append(self.rhoext.evaluate(t))
        rho = np.array(rho)
        
        #initalize solution container
        self.solution = \
            pointkineticscontainer(timepoints=self.timepoints, nt=nt, dnt=dnt,\
                power=power, rho=rho, flux=flux, typ="PKE")
    
    
    def solve(self, rtol=1e-10):
        """function defined kinetics"""
        
        _ispositive(rtol, "relative tolerance of time derivative integrator")
        
        # obtain initial conditions
        x0 = self.__solveinitialconditions()
        
        # obtain solution over defined transient
        solution = odeint(self.__dxdt, x0, self.timepoints, rtol=rtol)
        
        # run post-processing features
        self.__postprocess(solution)


# -----------------------------------------------------------------------------
# ----- source driven point kinetics solver class 
# -----------------------------------------------------------------------------


KEYS_SRC_PKE = {"beta": [list, float, "delayed neutron group yields, unitless"],
                
                "lamda": [list, float, "delayed neutron group decay constants, in"
                          " units of 1/seconds"],
            
                "promptL": [float, None, "neutron mean generation time, in units"
                            " of seconds"],
                
                "rhoi": [float, None, "initial excess reactivity, in units "
                         "of dk/k"],
                
                "S0": [float, None, "source activity, in units of Bq"],
            
                "epsilon": [float, None, "Source efficiency, in units of "
                            "fissions/source emission"],
            
                "volume": [float, None, "volume of fissioning region, in units of m3"],
            
                "Q": [float, None, "average recoverable energy per fission, in "
                      "units of MeV/fission"],
            
                "rhoext": ["class", None, "function that defines excess reactivity as"
                           " a function of time. Should be an instance of the "
                           "general control rules class"],
            
                "timepoints": ["np array", float, "time points for which results "
                               "are reported by solver"],
            
                "rtol": [float, None, "relative error convergence criteria for "
                         "inital conditions and time derivative evaluation"]}


REQ_SRC_PKE = ['beta', 'lamda', 'promptL', 'rhoi', 'S0', 'epsilon', 'rhoi',
               'volume', 'Q', 'rhoext', 'timepoints']


class srcpke:
    
    
    def __init__(self, **kwargs):
        """function initalizes source driven point kinetics solver class"""
        
        self.__dict__.update(kwargs)
        self.__checkinputs()
        self.groupsDN = len(self.beta)
        self.betaTot = self.beta.sum()

    
    def __checkinputs(self):
        """function does basic error checking"""
        pass
    
    
    def __findinitalpopulation(self, ni):
        """function calculates inital reactor power using newton-raphson method"""
        
        delayed = np.sum(self.beta / (self.promptL) ) * ni
        
        prompt = ((self.rhoi - self.betaTot)/self.promptL) * ni
        
        return prompt + delayed + self.srcnt
    
    
    def __solveinitialconditions(self, rtol):
        """function solves initial conditions"""
        
        #compute neutron production rate by source
        self.srcnt = self.nubar * self.S0 * self.epsilon # * self.promptL
        
        #numerically compute inital neutron population
        convrg = root(self.__findinitalpopulation, 1.0, tol=rtol,
                      method="broyden1", options={"maxiter": 1000})
        self.n0 = convrg.x
        
        #if solution could not be found, fail run
        if convrg.success is False:
            raise ValueError("inital neutron population solution did not converge!")
        
        #initialize inital conditions
        x0 = self.n0*np.append(1, self.beta/self.lamda/self.promptL)
        
        return x0
        
    
    
    def __constructmatrix(self, t):
        """function constructs matrix for estimation of time derv."""
        
        #construct matrix
        mtxA = np.zeros((self.groupsDN+1, self.groupsDN+1))
        
        #fill diagonal with decay constants
        np.fill_diagonal(mtxA, np.append(0, -self.lamda))
        
        mtxA[0,:] = \
            np.append((self.rhoext.evaluate(t)+self.rhoi-self.betaTot)/self.promptL, self.lamda)
        
        # build the first column
        mtxA[1:self.groupsDN+1, 0] = self.beta / self.promptL
        
        return mtxA
    
    
    def __dxdt(self, x0, t):
        """function evaluates the time derv. of the system"""
        
        dxdt = np.dot(self.__constructmatrix(t), x0)
        dxdt[0] += self.srcnt
        
        return dxdt
    
    
    def __postprocess(self, solution):
        """function performs post processing routines"""
        
        #seperate out delayed neutron precusors
        dnt = solution[:,1:].T
        
        #seperate out total neutron population
        nt = solution[:,0]
        
        #convert neutron population to flux
        flux = (self.v * nt) / self.volume
        flux /= 10000 # convert n/m2-s to n/cm2-s
        
        #convert neutron population into power
        power = (nt * self.Q * MEV_2_J) / (self.nubar * self.promptL)
        
        #obtain time dependent external reactivity
        rho = []
        for t in self.timepoints: rho.append(self.rhoext.evaluate(t))
        rho = np.array(rho) + self.rhoi
        
        #initalize solution container
        self.solution = \
            pointkineticscontainer(timepoints=self.timepoints, nt=nt, dnt=dnt,\
                power=power, rho=rho, flux=flux, typ="SPKE")
    
    
    def solve(self, rtol=1e-10):
        """function defined kinetics"""
        
        _ispositive(rtol, "relative tolerance of time derivative integrator")
        
        # obtain initial conditions
        x0 = self.__solveinitialconditions(rtol)
        
        # obtain solution over defined transient
        solution = odeint(self.__dxdt, x0, self.timepoints, rtol=rtol)
        
        # run post-processing features
        self.__postprocess(solution)
    
    