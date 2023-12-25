"""odeintsolver

Integrate a system of ordinary differential equations with RK45 adaptive time
mesh scheme


"""

import numpy as np
from scipy.integrate import odeint



class odeintSolver:
    """Solve using scipy odeint RK45 adaptive time mesh scheme"""
    
    def __init__(self, rtol=1E-10):
        """function initalizes odeint solver
        

        Parameters
        ----------
        rtol : float, optional
            Relative difference convergence criteria, default is 1E-10.

        Returns
        -------
        The matrix-vector product of the PKE matrix and initial vector x.

        """

        self.rtol = rtol
    
    def dxdt(self, x0, t, mtx):
        """function produces time rate of change for each variable in vec(x)"""
        return np.dot(mtx, x0)
    
    def solve(self, x0, mtx, t):
        """solve change in concentration"""
        return odeint(self.dxdt, x0, t, args=(mtx,), rtol=self.rtol)
    
    
class odeintSolverwSoruce:
    """Solve using scipy odeint RK45 adaptive time mesh scheme"""
    
    
    def __init__(self, rtol=1E-10):
        """function initalizes odeint solver
        

        Parameters
        ----------
        rtol : float, optional
            Relative difference convergence criteria, default is 1E-10.

        Returns
        -------
        The matrix-vector product of the PKE matrix and initial vector x.

        """

        self.rtol = rtol
    
    
    def dxdt(self, x0, t, mtx):
        """function produces time rate of change for each variable in vec(x)"""
        
        return np.dot(mtx, x0)
    
    
    def solve(self, x0, mtx, t):
        """solve change in concentration"""
        return odeint(self.dxdt, x0, t, args=(mtx,), rtol=self.rtol)
    