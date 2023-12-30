"""control.py

Define rules to control the dynamic external reactivity given to kinetics
solvers

@author: Matt Krecicki
@email: matthewkrecicki@gmail.com

"""


import typing
import numpy as np
from ntpSystem.errors.checkerrors import _isequallength, _inlist, _islist,\
    _issortedarray, _ispositiveArray
from ntpSystem.functions.state import State


# Possible functions:
    # linear:: ax + b
    # exponential:: a*exp(bx)+c
    # polynominal a*x^n + b*x^(n-1) .....
SHAPE_FUNCTIONS = ["linear", "exponential", "polynominal"]

NUM_COEFF = {"linear": 2, "exponential": 3, "polynominal": None}

# link between the functions and its intergal function
INTERGRAL_FUNCS = {"linear": "linearIntegral",
                   "exponential": "exponentialIntegral",
                   "polynominal": "polynominalIntegral"}


class shapeFunctions:
    """pre-defined shape functions"""
    
    def __init__(self):
        pass
    
    
    def linear(self, x, coeff):
        """linear function"""
        return coeff[0]*x + coeff[1]
    
    
    def linearIntegral(self, x, coeff):
        """integral of linear function"""
        return (coeff[0]/2)*x**2 + coeff[1]*x
    
    
    def exponential(self, x, coeff):
        """linear function"""
        return coeff[0]*np.exp(coeff[1]*x) + coeff[2]
    
    
    def exponentialIntegral(self, x, coeff):
        """integral of linear function"""
        return (coeff[0]/coeff[1])*np.exp(coeff[1]*x) + coeff[2]*x   
    
    
    def polynominal(self, x, coeff):
        """polynomial function"""
        return np.polyval(coeff, x)
    
    
    def polynominalIntegral(self, x, coeff):
        """integral of polynominal function"""
        return np.polyint(coeff)(x)
        
    
    
class generalControlRule(shapeFunctions):
    """Define a general control rule

    The general control rule can be applied to control
    chamber temperature and pressure only. Different time intervals can have
    different control rules. The transient always begins at 0.0 seconds.


    Parameters
    ----------
    funcTypes : list fo strings
        types of control rules, e.g., ['linear', 'linear', 'exponential']
        for each time interval.
    coeffs : list of lists
        coefficients of the corresponding function type for each time interval
    tends : list of floats
        time at which each interval ends (time bounds) without the 0.0

    Raises
    ------
    ValueError
        If ``coeffs`` or ``tends`` are not the same length as ``funcTyps``.
    TypeError
        If ``funcTyps``, ``coeffs`` or ``tends`` are not lists.


    Examples
    --------
    >>> generalControlRule(['linear', 'linear'], [[1.0, 2.0], [1.5, 2.0]],
                           [30.0, 60.0])

    """    
    
    
    def __init__(self, funcTypes, coeffs, tends):
        
        self.funcTypes = funcTypes
        self.coeffs = coeffs
        self.tends = tends  # time interval bounds
        self.__checkInputs()
        
        self.tends = np.array(self.tends)
    
    
    def __checkInputs(self):
        """check validity of variables"""
        
        # variables must be of a list type
        _islist(self.funcTypes, 'Function types')
        _islist(self.coeffs, 'Function coefficients')
        _islist(self.tends, 'Intervals end times (tends)')
        
        # lists must be all the same length
        nIntervals = len(self.funcTypes)  # number of time intervals
        _isequallength(self.coeffs, nIntervals, "Function coefficients")
        _isequallength(self.tends, nIntervals, "Intervals end times (tends)")
        
        
        for idx, fType in enumerate(self.funcTypes):
            _inlist(fType, "Function type (idx={})".format(idx),
                    SHAPE_FUNCTIONS)
            exp_ncoeff = NUM_COEFF[fType]  # expected number of coefficients
            _isequallength(self.coeffs[idx], exp_ncoeff,
                           "Function coefficients")
          
        tendsArray = np.array(self.tends)
        _issortedarray(tendsArray, "Intervals end times (tends)")
        _ispositiveArray(tendsArray, "Intervals end times (tends)")
    
    
    def evaluate(self, t):
        """evaluate of the function at a specific time point"""
        
        # Default value of the function is zero
        y = 0.0
        
        if t <= np.max(self.tends):                
            idx = np.where(t <= self.tends)[0][0]
            y = getattr(self, self.funcTypes[idx])(t, self.coeffs[idx])
                    
        return y
    
    
    def evaluateIntegral(self, Tstart, Tend):
        """evaluate the value of the integral using knowing limits"""
        
        y = 0.0
        
        if Tstart <= np.max(self.tends):
            
            idx0 = np.where(Tstart <= self.tends)[0][0]
            if Tend > self.tends[-1]:
                idx1 = len(Tend)
            else:
                idx1 = np.where(Tend <= self.tends)[0][0]
            
            
            intgType = INTERGRAL_FUNCS[self.funcTypes[idx0]]
            # evaluate the integral for Tstart
            y0 = getattr(self, intgType)(Tstart, self.coeffs[idx0])
 
            # start and end time are in the same interval
            if idx0 == idx1:
                y1 = getattr(self, intgType)(Tend, self.coeffs[idx0])
                # return the value of the integral
                return y1 - y0    
            
            # start and end time are NOT in the same interval
            # Loop over all the relevant intervals
            for idx in range(idx0, idx1+1):
                if idx > idx0:
                    intgType = INTERGRAL_FUNCS[self.funcTypes[idx]]
                    y0 = getattr(self, intgType)(self.tends[idx-1],
                                                 self.coeffs[idx])
                if idx < idx1:
                    intgType = INTERGRAL_FUNCS[self.funcTypes[idx]]
                    y1 = getattr(self, intgType)(self.tends[idx],
                                                 self.coeffs[idx])
                else:
                    intgType = INTERGRAL_FUNCS[self.funcTypes[idx]]
                    y1 = getattr(self, intgType)(Tend,
                                                 self.coeffs[idx])                    
                y += (y1 - y0)
                    
        return y 
    

