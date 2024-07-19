# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 16:36:11 2024

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""

import numpy as np
from scipy.interpolate import interp1d 
from kinetics.functions.plotters import lineplot



def drumworth(theta, integralworth, thetacritical=70.0, plot=False):
    """function return excess reactivity of control drum based on input angle
    and integral worth of drum. This function assumes that 0-degrees represents
    control drums fully inserted and 180-degrees represents control drums fully
    withdrawn.
    

    Parameters
    ----------
    theta : float
        control drum angle, in units of degrees.
    integralworth : float
        integral drum worth, in units of dk/k.
    thetacritical : float
        position where reactivity is equal to zero, in units of degrees.
    plot : bool
        flag that if True will plot integral drum worth curve. Typically used
        for error checking purposes. 

    Returns
    -------
    None.

    """
    
    
    #generate forward function
    theta = np.linspace(180.0, 0.0, 50)
    curve = (-integralworth / (1 + np.exp(-1*np.linspace(-4.0, 4.0, 50)))) 
    rhofunc = interp1d(theta, curve)
    rhovals = rhofunc(theta) - rhofunc(thetacritical)
    rhofunc = interp1d(theta, rhovals)
    
    if plot:
        lineplot([theta], [rhofunc(theta)], 
                 xlabel="Drum Position, degrees", ylabel="Reactivity, dk/k",
                 grid=True, linestyles=["-"], markers=["None"], colors=["black"])
    
    #generate inverse function
    rhoinvfunc = interp1d(rhovals, theta)
    
    return rhofunc, rhoinvfunc
    
    

integralworth = 10 * 0.00689 # ten dollars of reactivity 

drumworth, invdrumworth = drumworth(0.0, integralworth, thetacritical=90.0, plot=True)


rhoi = np.array([0.0, 0.0, 0.003, 0.006, -0.003])
thetai = invdrumworth(rhoi)




