# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 22:07:53 2024

@author: matt krecicki
@email: matthewkrecicki@gmail.com

function solvers inverse point kinetics problems

formulation is based on:
    T. Ball, "The Inverse Kintics Method and Its Application to the Annular Core 
    Research Reactor", University of New Mexico Master's Thesis, (2017).

"""


MEV_2_J = 1.60218e-13

import numpy as np
import h5py
import functools
from scipy.integrate import quad
from scipy.integrate import odeint
from scipy.interpolate import UnivariateSpline
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
    
    
    def __init__(self, inputs, noSamples=1000):
        """function initalizes inverse point kinetics solver"""
        
        self.inputs = inputs
        self.__checkinputs()
        setattr(self.inputs, "groupsDN", len(self.inputs.beta))
        setattr(self.inputs, "betaTot", self.inputs.beta.sum())
        
        #evaluate 1st and 2nd derivatives
        ti = np.linspace(self.inputs.timepoints[0], self.inputs.timepoints[-1],
                         noSamples)
        
        ni = \
            (self.inputs.power(ti)*self.inputs.nubar*self.inputs.promptL)/ \
                (self.inputs.Q * MEV_2_J)
        
        self.inputs.deriv1 = \
            UnivariateSpline(ti, ni, k=1).derivative(n=1)
        self.inputs.deriv2 = \
            UnivariateSpline(ti, ni, k=2).derivative(n=2)
        
        test1 = self.inputs.deriv1(self.inputs.timepoints)
                
    
    def __solveinitialconditions(self):
        """function solves initial conditions"""

        #convert power to total neutron population
        P0 = float(self.inputs.power(self.inputs.timepoints[0]))
        
        n0 = \
            (P0*self.inputs.nubar*self.inputs.promptL)/(self.inputs.Q * MEV_2_J)
        
        self.inputs.n0 = n0
        
        #detemerine inital reactivity and precusor concentrations
        x0 = \
            n0*np.append(0.0, self.inputs.beta/self.inputs.lamda/\
                                     self.inputs.promptL)
        
        return x0
    
    
    def __delayedderivative(self, Ci, Ni):
        """function computs the neutron precusor time derivative"""
        
        dCdt = \
            (self.inputs.beta/self.inputs.promptL)*Ni - (Ci*self.inputs.lamda)
        
        return dCdt
    
    
    def __constructmatrix(self, t, state_history):
        """function constructs matrix for estimation of time derv."""
        
        Ni = \
            (float(self.inputs.power(t))*self.inputs.nubar*self.inputs.promptL)/ \
                (self.inputs.Q * MEV_2_J)
        
        #get neutron precusor density
        Ci = state_history["prev"][1:]
        
        #compute the reactivity time derivative
        dNdt = float(self.inputs.deriv1(t))
        dN2dt2 = float(self.inputs.deriv2(t)) 
        
        #compute delayed neutron precusor time derivative
        dCdt = self.__delayedderivative(Ci, Ni)
        
        ratiotest = self.inputs.x0[1:]/Ci
        
        #evaluate coefficients for reactivity time derivative
        a = self.inputs.promptL / float(self.inputs.power(t))**2
        b = -dNdt**2
        c = np.sum(self.inputs.lamda*Ci*dNdt)
        d = dNdt * dN2dt2
        e = np.sum(Ni * dCdt * self.inputs.lamda)
        
        drhodt = a * (b + c + d - e)
        
        return drhodt, dCdt
    
    
    def __dxdt(self, x0, t, state_history):
        """function evaluates the time derv. of the system"""
        
        state_history["prev"] = x0
        
        drhodt, dCdt = self.__constructmatrix(t, state_history)        
        
        dxdt = np.zeros(self.inputs.groupsDN+1)
        dxdt[0] = drhodt
        dxdt[1:] = dCdt
                
        return dxdt
        
    
    def solve(self, rtol=1e-10):
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
        factor : float, optional
            fraction of initial time step to look back to compute instaneous
            derivative of power. The default is 1e-3

        Returns
        -------
        None.

        """
        
        _ispositive(rtol, "relative tolerance of time derivative integrator")
        
        # obtain initial conditions
        x0 = self.__solveinitialconditions()
        self.inputs.x0 = x0
        
        # obtain solution over defined power trace
        dxdt = functools.partial(self.__dxdt, state_history={})
        
        solution = odeint(dxdt, x0, self.inputs.timepoints, rtol=rtol)
        
        #generate output object
        self.__postprocess(solution)
    
    
    def __postprocess(self, solution):
        """function performs post processing routines"""
        
        #seperate out delayed neutron precusors
        dnt = solution[:,1:].T
        
        #seperate out total neutron population
        rho = solution[:,0]
        
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
                nt=nt, dnt=dnt, power=power, flux=flux, rho=rho,
                typ="invpke")


class inversepke_old:
    
    
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
                
        return total, res
    
    
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
        
        return total, res
    
    
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
        
        rhop = coeff * promptTerm
        rhod = coeff * (1 - delayedTerm)
        
        return rhot, rhop, rhod
    
    
    def solve(self, intrtol=1.49e-08, intabstol=1.49e-08, factor=1e-3):
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
        factor : float, optional
            fraction of initial time step to look back to compute instaneous
            derivative of power. The default is 1e-3

        Returns
        -------
        None.

        """
        
        _ispositive(intrtol, "relative tolerance of time derivative integrator")
        _ispositive(intabstol, "absolute tolerance of time derivative integrator")
        _ispositive(factor, "factor to look back at previous time step to compute derivative")
        
        #solve power time derviative
        dPdt, tmid = self.__evaluatederivative(factor)
        
        #solve integrals
        intDt, grpDt = self.__evaluateintegral(intrtol, intabstol)
        
        #evaluate delayed terms
        intD0, grpD0 = self.__computedelayedcomponet(intDt)
        
        #evalute reactivity
        rhot, rhop, rhod = self.__evaluatereactivity(dPdt, intD0, intDt)
                
        #generate output object
        self.__postprocess(rhot, rhop, rhod)
    
    
    def __postprocess(self, rhot, rhop, rhod):
        """function performs post processing routines"""
                        
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
                nt=nt, power=power, rho=rhot,
                 flux=flux, typ="invpke")
        

SRC_INV_PKE_DICT = \
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
     
     "S0": [float, None, "initial neutron source strength", "neutrons/second",
            True],
     
     "epsilon": [float, None, "neutron source efficiency",
                 "fissions/neutron emitted", True],
     
     "typ": [str, None, "type of kinetic simulation desired", "n/a", False]}

        
class srcinversepke:
    
    
    def __checkinputs(self):
        """function runs basic errpr checking on source driven inverse point
        kinetics"""
        
        if self.inputs.typ != "srcinvpke":
            raise TypeError("incorrect inputs container given: {}"\
                            .format(self.inputs.typ))
        
        _pkecheck(self.inputs)
        _checkdict(SRC_INV_PKE_DICT, self.inputs)
        
    
    def __init__(self, inputs):
        """function initalizes inverse point kinetics solver"""
        
        self.inputs = inputs
        self.__checkinputs()
        setattr(self.inputs, "groupsDN", len(self.inputs.beta))
        setattr(self.inputs, "betaTot", self.inputs.beta.sum())
        
        
        #TODO: clean all this up after verification!
        
        #compute contribution of neutrons from external source
        self.inputs.srcnt = \
            self.inputs.nubar*self.inputs.S0*self.inputs.epsilon
        
        #convert src neutrons into power
        self.inputs.srcnt = \
            (self.inputs.srcnt*self.inputs.Q*MEV_2_J) / \
            (self.inputs.nubar * self.inputs.promptL)
        
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
                
        return total, res
    
    
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
        
        return total, res
    
    
    def __evaluatereactivity(self, dPdt, decD0, intDt):
        """evaluates total reactivity as a function of time"""
        
        invPower = 1 / self.inputs.power(self.inputs.timepoints)
        
        # ----- compute coefficient
        c1 = self.inputs.promptL * dPdt * invPower
        c2 = self.inputs.promptL*self.inputs.srcnt*invPower
        coeff = self.inputs.betaTot / (1 + c1 - c2)
        
        # ----- compute prompt term
        promptTerm = \
            (self.inputs.promptL / self.inputs.betaTot) * dPdt * invPower
        
        # ----- compute delayed term
        delayedTerm = invPower*(decD0 + intDt)
        
        # ----- compute external source term
        extsrcTerm = \
            (self.inputs.promptL) * (self.inputs.srcnt * invPower)
        
        # ----- compute reactivity trace
        rhot = coeff * (1 + promptTerm - delayedTerm - extsrcTerm)
        rhot += self.inputs.rho0
        
        return rhot
    
    
    def solve(self, intrtol=1.49e-08, intabstol=1.49e-08, factor=1e-3):
        """function solves source driven inverse kinetics problem to obtain
        total reactivity as a function of time based on a given power trace
        

        Parameters
        ----------
        intrtol : float, optional
            relative convergence tolerance for numerical integration. The
            default is 1.49e-08.
        intabstol : float, optional
            absolute tolerance for numerical integration. The default is
            1.49e-08.
        factor : float, optional
            fraction of initial time step to look back to compute instaneous
            derivative of power. The default is 1e-3.

        Returns
        -------
        None.

        """
        
        _ispositive(intrtol, "relative tolerance of time derivative integrator")
        _ispositive(intabstol, "absolute tolerance of time derivative integrator")
        _ispositive(factor, "factor to look back at previous time step to "
                    "compute derivative")
        
        #solve power time derviative
        dPdt, tmid = self.__evaluatederivative(factor)
        
        #solve integrals
        intDt, grpDt = self.__evaluateintegral(intrtol, intabstol)
        
        #evaluate delayed terms
        intD0, grpD0 = self.__computedelayedcomponet(intDt)
        
        #evalute reactivity
        rhot = self.__evaluatereactivity(dPdt, intD0, intDt)
                
        #generate output object
        self.__postprocess(rhot)
    
    
    def __postprocess(self, rhot):
        """function performs post processing routines"""
                        
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
                nt=nt, power=power, rhototal=rhot, flux=flux, typ="srcinvpke")
        
    
        
        