"""solve PKEs

Function that solves point kinetic equations with arbitrary number
of delay groups.


"""

from odeintsolver import odeintSolver
from generatepke import GeneratePKEs, GeneratePKEwSoruce


def SolvePKEs(rho, beta, lamda, promptL, P0, timepoints, rtol=1E-10):
    """Solve the PKEs

    Parameters
    ----------
    timepoints : ndarray
        Absolute time points in seconds, e.g. [0, 0.1, 0.2, 0.3, ...]
    rho : ndarray
        constant reactivity
    beta : 1-dim ndarray
        delayed neutron fractions
    lamda : 1-dim ndarray
        delayed neutron decay constants
    prompt : float
        prompt neutron generation time in seconds
    P0 : float
        initial flux or power

    Returns
    -------
    xt : 2-dim ndarray
        time-dependent vector

    """
    
    # generate the initial vector and reset the PKE matrix
    x0, mtx = GeneratePKEs(rho, beta, lamda, promptL, P0)
    
    # reset odeintSolver
    pke = odeintSolver(rtol)
    
    # obtain the solution for x
    xt = pke.solve(x0, mtx, timepoints)
    
    return xt


def SolvePKEwSource(rhoi, beta, lamda, promptL, S0, Q, epsilon,
                    timepoints, rtol=1E-10):
    """Solve the PKEs

    Parameters
    ----------
    
    rhoi : float
        initial negative reactivity
    beta : 1-dim ndarray
        delayed neutron fractions
    lamda : 1-dim ndarray
        delayed neutron decay constants
    prompt : float
        prompt neutron generation time in seconds
    S0 : float
        initial source strength in units of neutrons/second.
    timepoints : ndarray
        Absolute time points in seconds, e.g. [0, 0.1, 0.2, 0.3, ...]
        
        
    Returns
    -------
    xt : 2-dim ndarray
        time-dependent vector

    """
    
    # generate the initial vector and reset the PKE matrix
    x0, mtx = GeneratePKEwSoruce(rhoi, beta, lamda, promptL, S0, Q, epsilon)
    
    # reset odeintSolver
    pke = odeintSolver(rtol)
    
    # obtain the solution for x
    xt = pke.solve(x0, mtx, timepoints)
    
    return xt
    
    

