# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 10:34:06 2024

@author: matt krecicki
@email

"""

import numpy as np
from scipy.interpolate import interp1d 
from kinetics.functions.plotters import lineplot
from kinetics.errors.checkerrors import _inrange, _isbool


CONTROL_DRUM_INPS = {"integralworth": [],
                     "thetacritical": [],
                     "plot": []}


class controldrum:
    
    
    def __checkinputs(self):
        """function runs basic error checking on user defined inputs"""
        pass
    
    
    def __init__(self, **kwargs):
        """function initalizes feedback class"""
        
        self.__dict__.update(kwargs)
        self.__checkinputs()
        self.__generateworthcurve()
        
    
    def __generateworthcurve(self):
        """function generates intergal worth curve based on sigmoid function"""
        
        #generate forward function
        theta = np.linspace(180.0, 0.0, 50)
        curve = (-self.integralworth / (1 + np.exp(-1*np.linspace(-4.0, 4.0, 50)))) 
        rhofunc = interp1d(theta, curve, )
        rhovals = rhofunc(theta) - rhofunc(self.thetacritical)
        self.rhofunc = interp1d(theta, rhovals, fill_value="extrapolate")
        
        if self.plot:
            lineplot([theta], [self.rhofunc(theta)], 
                     xlabel="Drum Position, degrees", ylabel="Reactivity, dk/k",
                     grid=True, linestyles=["-"], markers=["None"], colors=["black"])
        
        #generate inverse function
        self.rhoinvfunc = interp1d(rhovals, theta, fill_value="extrapolate")
        
        #compute reactivity limits
        self.__rholimts = [self.rhoinvfunc(180.0), self.rhoinvfunc(0.0)]
        
    
    def evaluate(self, val, inverse=False):
        """function evaluates the reactivity or drum position based on the user
        input
        

        Parameters
        ----------
        val : floats
            control drum position to evaluate reactivity, in units of degrees.
        inverse : bool, optional
            flag to evaluate inverse function where the input becomes reactivity 
            in units of dk/k and the output is control drum position in units 
            of degrres. The default is False.

        Returns
        -------
        value : float
            DESCRIPTION.

        """
        _isbool(inverse, "inverse reactivity calcualtion flag")
        
        if inverse is False:
            _inrange(val, "control drum position: {}".format(val), [0.0, 180.0])
            return self.rhofunc(val)
        else:
            _inrange(val, "reactivity limits", self.__rholimits)
            return self.rhoinvfunc(val)


class feedback:
    
    def __checkinputs(self):
        """function runs basic error checking on user defined inputs"""
        pass
    
    def __init__(self, **kwargs):
        """function initalizes feedback class"""
        
        self.__dict__.update(kwargs)
        self.__checkinputs()
        
        