# -*- coding: utf-8 -*-
"""reactivity correlations container 


@author: matt krecicki

"""

CORR_TYPES = ["polynominal", "power", "expontial", "saturation"]

import numpy as np
from abc import ABC, abstractmethod
from ntptransient.checkerrors import _inlist, _islist, _isnumber


def initalizeCorrelation(typ, coeffs, limits):
    
    _inlist(typ, "correlation type", CORR_TYPES)
        
    if typ == "polynominal":
        corr = polynominal(coeffs, limits)
        
    elif typ == "power":
        corr = power(coeffs, limits)
        
    elif typ == "expontial":
        corr = expontial(coeffs, limits)
        
    elif typ == "saturation":
        corr = saturation(coeffs, limits)
    
    return corr


class correlationTemplate(ABC):


    def __init__(self, coeff, limits):
        """functions initalizes reactivity correlation and performs basic error
        checking
        

        Parameters
        ----------
        coeff : list
            reactivity correlation coefficients.
        limits : list
            lower and upper bounds for reactivity correlation.

        Returns
        -------
        None.

        """
        _islist(coeff, "reactivity correlation coefficient list")
        for num in coeff:
            _isnumber(num, "reactivity correlation coefficient")
            
        _islist(limits, "reactivity correlation limits")
               
        self.coeff = coeff
        self.limits = limits
    
    
    def _correctLimits(self, x):
        """function corrects the input to be within the bounds of the
        correlation
        

        Parameters
        ----------
        x : float
            reactivity correlation input value.

        Returns
        -------
        x : float
            reactivity correlation input corrected value.

        """
        if x > self.coeff[1]: x = self.coeff[1]
        elif x < self.coeff[0]: x = self.coeff[0]
        
        return x
    
    
    @abstractmethod
    def evaluate(self, x):
        """function evaluates the reactivity correlation"""


class polynominal(correlationTemplate):
        
    
    def evaluate(self, x):
        """function evaluates the reactivity correlation for a given input value
        

        Parameters
        ----------
        x : float
            reactivity correlation input value.

        Returns
        -------
        alpha : float
            output reactivity  coefficient.

        """
        x = self._correctLimits(x)
        alpha = np.polyval(self.coeff, x)
        
        return alpha
    
    
class power(correlationTemplate):

    
    def evaluate(self, x):
        """function evaluates the reactivity correlation for a given input value
        

        Parameters
        ----------
        x : float
            reactivity correlation input value.

        Returns
        -------
        alpha : float
            output reactivity  coefficient.

        """
        x = self._correctLimits(x)
        alpha = self.coeff[0] 
        
        return alpha
    

class expontial(correlationTemplate):
    
    
    def evaluate(self, x):
        pass
    
    
class saturation(correlationTemplate):
    
    
    def evaluate(self, x):
        pass
    
        
    
    
    
    