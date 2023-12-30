# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 12:32:15 2023

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""

MEV_2_J = 1.60218e-13

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import root
from kinetics.errors.checkerrors import _isnegative, _inrange
from kinetics.containers.outputs import pointkineticsoutputs


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


class pointkinetics:
    
    
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
        
        #convert power to flux
        self.n0 = () / (self.Q * MEV_2_J)
        
        #detemerine inital reactor power and precusor concentrations
        x0 = self.n0*np.append(1, self.beta/self.lamda/self.promptL)
        
        return x0
        
    
    def __constructmatrix(self, t):
        """function constructs matrix for estimation of time derv."""
        
        #construct matrix
        mtxA = np.zeros((self.groupsDN+1, self.groupsDN+1))
        
        #fill diagonal with decay constants
        np.fill_diagonal(mtxA, np.append(0, -self.lamda))
            
        mtxA[0,:] = np.append((self.rhoext(t)-self.betaTot)/self.promptL, self.lamda)
        
        # build the first column
        mtxA[1:self.groupsDN+1, 0] = self.beta / self.promptL
        
        return mtxA
    
    
    def __dxdt(self, x0, t):
        """function evaluates the time derv. of the system"""
        
        return np.dot(self._generatematrix(x0, t), x0)
    
    
    def solve(self):
        """function defined kinetics"""
        
        # obtain initial conditions
        x0 = self.__solveinitialconditions()
        
        # obtain solution over defined transient
        self.solution = odeint(self._dxdt, x0, self.timepoints)
        
        # send array to solution object
        pointkineticsoutputs()


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


REQ_SRC_PKE = ['beta', 'lamda', 'promptL', 'rhoi', 'S0', 'epsilon', 'volume',
               'Q', 'rhoext', 'timepoints']


class sourcepointkinetics:
    
    
    def __init__(self, **kwargs):
        """function initalizes source driven point kinetics solver class"""
        
        self.__dict__.update(kwargs)
        self.__checkinputs()

    
    def __checkinputs(self):
        """function does basic error checking"""
        pass
    
    
    def __solveinitialconditions(self):
        """function solves initial conditions"""
        pass
    
    
    def __constructmatrix(self):
        """function constructs matrix for estimation of time derv."""
        pass
    
    
    def __dxdt(self, t):
        """function evaluates the time derv. of the system"""
        pass
    
    
    def solve(self):
        """function defined kinetics"""
        pass
    
    