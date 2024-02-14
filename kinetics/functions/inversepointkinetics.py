# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 22:07:53 2024

@author: matt krecicki
@email: matthewkrecicki@gmail.com

function solvers inverse point kinetics problems

"""


MEV_2_J = 1.60218e-13

import numpy as np
from scipy.integrate import quad
from kinetics.errors.checkerrors import _ispositive
from kinetics.errors.customerrors import _checkdict, _pkecheck 
from kinetics.containers.outputs import pointkineticsOutputsContainer


INV_PKE_DICT = \
    {"beta": [np.ndarray, float, "delayed neutron fraction", "unitless", True],
     
     "lamda": [np.ndarray, float, "delay neutron group decay constant",
               "1/seconds", True],
     
     "promptL": [float, None, "prompt neutron lifeime", "seconds", True],
          
     "volume": [float, None, "total volume of reactor control volume",
                "meters^3", True],
     
     "nubar": [float, None, "average number of neutrons produced per fission",
               "neutrons/fission", True],
     
     "Q": [float, None, "Average recoverable energy released per fission",
           "MeV/fission", True],
     
     "v": [float, None, "effective one-group neutron velocity", "meters/second",
           True],
     
     "timepoints": [np.ndarray, float, "time points to return solution",
                    "seconds", False],
     
     "power": [object, None, "reactor power control class", "n/a",
               False],
     
     "typ": [str, None, "type of kinetic simulation desired", "n/a", False]}


class inversepke:
    
    
    def __checkinputs(self):
        """function runs basic errpr checking on inverse point kinetics"""
        
        if self.inputs.typ != "invpke":
            raise TypeError("incorrect inputs container given: {}"\
                            .format(self.inputs.typ))
        
        _pkecheck(self.inputs)
        _checkdict(INV_PKE_DICT, self.inputs)
    
    
    def __init__(self, inputs):
        """function initalizes inverse point kinetics solver"""
        
        self.inputs = inputs
        self.__checkinputs()
        setattr(self.inputs, "groupsDN", len(self.inputs.beta))
        setattr(self.inputs, "betaTot", self.inputs.beta.sum())
        
        self.__checkinputs()
        
    
    def __evaluatederivative(self, factor=1e-3):
        """function evaluates the power derviative as a function of time"""
                
        tf = self.inputs.timepoints
        ti = tf - (tf*factor)  
        
        tmid = 0.5*(tf+ti)
        
        Pf = self.inputs.power(tf)
        Pi = self.inputs.power(ti)
        
        dPdt = (Pf - Pi) / (tf - ti)
        
        #always zero
        dPdt[0] = 0.0
        
        return dPdt, tmid
        
    
    def __integralfunctions(self, t, ti, grp):
        
        coeff = (self.inputs.lamda[grp] * self.inputs.beta[grp]) / self.inputs.betaTot
        
        return coeff * np.exp(-self.inputs.lamda[grp]*t) * self.inputs.power(ti-t)
    
    
    def __evaluateintegral(self, intrtol, intabstol):
        """function evalutes integral using scipy.integrate.quad
        
        for details on numercial integration method see:
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html
        
        """
        
        res = {}
        
        #loop over all the groups
        for grp in range(self.inputs.groupsDN):
            
            intvals = []
            
            #loop over all the time points
            for ti in self.inputs.timepoints:
                
                intgrl, _err = \
                    quad(self.__integralfunctions, 0.0, ti, args=(ti, grp),
                         epsabs=intabstol, epsrel=intrtol,
                         limit=50, points=None, weight=None, wvar=None,
                         wopts=None, maxp1=50, limlst=50, complex_func=False)
                
                intvals.append(intgrl)
            
            #append group results to container
            res["group"+str(grp+1)] = np.array(intvals)
        
        #compute total
        total = np.zeros(len(self.inputs.timepoints))
        for key in list(res.keys()): total += res[key]
        res["total"]  = total
                
        return total
    
    
    def __computedelayedcomponet(self, intDn):
        """compute delayed neutron componet coefficient"""
        
        P0 = self.inputs.power(0.0)
        res = {}
        
        #loop over all delayed neutron groups
        for grp in range(self.inputs.groupsDN):
        
            res["group"+str(grp+1)] = \
                P0 * (self.inputs.beta[grp]/self.inputs.betaTot)*\
                np.exp(-self.inputs.lamda[grp]*self.inputs.timepoints)
        
        total = np.zeros(len(self.inputs.timepoints))
        for key in list(res.keys()): total += res[key]
        
        return total
    
    
    def __evaluatereactivity(self, dPdt, intD0, intDt):
        """evaluates total reactivity as a function of time"""
        
        invPower = 1 / self.inputs.power(self.inputs.timepoints)
        
        # ----- compute coefficient
        coeff = \
            self.inputs.betaTot / (1 + self.inputs.promptL * dPdt * invPower)
        
        # ----- compute prompt term
        promptTerm = \
            (self.inputs.promptL / self.inputs.betaTot) * dPdt * invPower
        
        # ----- compute delayed term
        delayedTerm = invPower*(intD0 + intDt)
        
        # ----- compute reactivity trace
        rhot = coeff * (promptTerm + 1 - delayedTerm)
        
        return rhot
    
    
    def solve(self, intrtol=1.49e-08, intabstol=1.49e-08):
        """function solves inverse kinetics problem to obtain total reactivity
        as a function of time based on a given power trace
        

        Parameters
        ----------
        intrtol : float, optional
            relative convergence tolerance for numerical integration. The
            default is 1.49e-08.
        intabstol : float, optional
            absolute tolerance for numerical integration. The default is
            1.49e-08.

        Returns
        -------
        None.

        """
        
        _ispositive(intrtol, "relative tolerance of time derivative integrator")
        _ispositive(intabstol, "absolute tolerance of time derivative integrator")
        
        #solve power time derviative
        dPdt, tmid = self.__evaluatederivative()
        
        #solve integrals
        intDt = self.__evaluateintegral(intrtol, intabstol)
        
        #evaluate delayed terms
        intD0 = self.__computedelayedcomponet(intDt)
        
        #evalute reactivity
        rhot = self.__evaluatereactivity(dPdt, intD0, intDt)
        
        #generate output object
        self.__postprocess(rhot)
    
    
    def __postprocess(self, rhot):
        """function performs post processing routines"""
        
        #TODO: execute forward solutin with external reactivity trace
                
        #compute power
        power = self.inputs.power(self.inputs.timepoints)
        
        #compute neutron population
        nt = \
            (power*self.inputs.nubar*self.inputs.promptL)/(self.inputs.Q*MEV_2_J)
        
        #compute neutron flux
        flux = (self.inputs.v * nt) / self.inputs.volume
        flux /= 10000 # convert n/m2-s to n/cm2-s
        
        #initalize solution container
        self.outputs = \
            pointkineticsOutputsContainer(timepoints=self.inputs.timepoints,\
                nt=nt, power=power, rho=rhot, flux=flux, typ="invpke")
        
        
        