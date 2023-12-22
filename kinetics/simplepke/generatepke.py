"""generate PKE matrix

Function that creates the point kinetic equations with arbitrary number
of delay groups.

"""

MEV_2_J = 1.60218e-13

import numpy as np
from kinetics.errors.checkerrors import _isnegative, _inrange


def GeneratePKEs(rho, beta, lamda, promptL, P0):
    """Generate the PKEs in a matrix form

    Parameters
    ----------
    timePoints : ndarray
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
    x0 : 1-dim ndarray
        An initial vector (at steady state)
    mtxA : 2-dim ndarray
        PKE matrix


    """

    # Number of delay neutron groups
    groupsDN = len(beta)

    # Calculate beta total:
    betaTot = beta.sum()

    # Initial conditions
    x0 = P0*np.append(1, beta/lamda/promptL)

    # reset the PKE matrix
    mtxA = np.zeros((groupsDN+1, groupsDN+1))

    # build diagonal
    np.fill_diagonal(mtxA, np.append(0, -lamda))

    # build the first row
    mtxA[0, :] = np.append((rho-betaTot)/promptL, lamda)

    # build the first column
    mtxA[1:groupsDN+1, 0] = beta / promptL

    return x0, mtxA


def GeneratePKEwSoruce(rhoi, beta, lamda, promptL, S0, Q, epsilon):
    """Generate the PKEs with external source in a matrix form

    Parameters
    ----------
    timePoints : ndarray
        Absolute time points in seconds, e.g. [0, 0.1, 0.2, 0.3, ...]
    rhoi : ndarray
        inital negativity reactivity of the system
    beta : 1-dim ndarray
        delayed neutron fractions
    lamda : 1-dim ndarray
        delayed neutron decay constants
    prompt : float
        prompt neutron generation, in units of seconds
    S0 : float
        source strength
    Q : float
        energy released per fission, in units of MeV.
    epsilon : float
        source efficiency, i.e. fraction of neutrons emitted by source that
        cause a fission event.

    Returns
    -------
    x0 : 1-dim ndarray
        An initial vector (at steady state)
    mtxA : 2-dim ndarray
        PKE matrix


    """
    
    #formulation only works for initally subcritical configuration
    _isnegative(rhoi, "inital negative reactivity in system")
    _inrange(epsilon, "external source efficiency", [0.0, 1.0], upBound=False,
             lowBound=True)
    
    # Number of delay neutron groups
    groupsDN = len(beta)
    
    # Calculate beta total:
    betaTot = beta.sum()
    
    #compute inital power based on source efficiency
    Q *= MEV_2_J
    
    Sp = (epsilon * S0 * Q) / promptL
    
    P0 = promptL * Sp / abs(rhoi)
    
    # Initial conditions
    x0 = P0*np.append(1, beta/lamda/promptL)
    
    # reset the PKE matrix
    mtxA = np.zeros((groupsDN+1, groupsDN+1))

    # build diagonal
    np.fill_diagonal(mtxA, np.append(0, -lamda))

    # build the first row
    mtxA[0, :] = np.append((rhoi-betaTot)/promptL, lamda)
    
    #include source power in matrix
    
    # build the first column
    mtxA[1:groupsDN+1, 0] = beta / promptL

    return x0, mtxA
