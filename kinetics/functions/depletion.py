# -*- coding: utf-8 -*-
"""Solve bateman equations for Xe, Sm, and U chains

Created on Wed Dec  1 18:12:37 2021

@author: matt krecicki

"""

import numpy as np
from scipy.linalg import expm
from scipy.integrate import odeint


def dndt(Nt: list, t: np.ndarray, power: float, volume: float, Tc: float, Pc: float):
    """function calculates time rate of change of each isotope using the 
    pade approximation. 
    

    Parameters
    ----------
    Nt : list
        isotopic concentration in atoms/cm3.
    t : float
        points in time, in units of seconds.
    power : float
        reactor power, in units of Watts.
    volume : float
        volume of fuel region, in units of cm3.
    Tc : float
        chamber temperature in units of kelvin.
    Pc : float
        chamber pressure in units of pascal.

    Notes
    -------
    the isotope arrays will always be indexed in the following way:
        [I-135, Xe-135, Pm-149, Sm-149, U-235, U-238]
    

    Returns
    -------
    dndt : np.ndarray
        time rate of change of each isotope in units of atoms/cm3-s.

    """
    #      531350,         541350,         611490,         621490,         922350       92238
    dcMtx = np.array([
        [-2.93061e-05, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000],  # 531350
        [2.93061e-05, -2.10657e-05, 0.000000000, 0.000000000, 0.000000000, 0.000000000],  # 541350
        [0.000000000, 0.000000000, -3.626825E-6, 0.000000000, 0.000000000, 0.000000000],  # 611490
        [0.000000000, 0.000000000, 3.626825E-6, 0.000000000, 0.000000000, 0.000000000],  # 621490
        [0.000000000, 0.000000000, 0.000000000, 0.000000000, -3.12085e-17, 0.000000000],  # 922350
        [0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, -4.91597e-18]  # 922380
    ])

    sigc, sigf, FY, Q = crossSections(Tc, Pc)
    flx = power2flux(power, volume, sigf, Nt, Q)
    trnsMtx = transmutationMatrix(sigc, sigf, FY)
    deplMtx = (dcMtx + trnsMtx * flx)

    dndt = np.dot(deplMtx, Nt)

    return dndt


def transmutationMatrix(siga: np.ndarray, sigf: np.ndarray, FY: np.ndarray):
    """function builds transmutation matrix
    

    Parameters
    ----------
    sigc : np.ndarray
        n,gamma microscopic cross sections, in units of cm2.
    sigf : np.ndarray
        fission microscopic cross sections, in units of cm2.
    FY : np.ndarray
        fission yields

    Returns
    -------
    trnsMtx : np.ndarray
        transmutation matrix.

    """
    trnsMtx = np.zeros((6, 6))  # create empty matrix
    np.fill_diagonal(trnsMtx, -1 * siga)  # add absorption cross section

    fissMtx = np.tile(sigf, (6, 1))
    fissMtx *= FY

    trnsMtx += fissMtx

    return trnsMtx


def power2flux(power: float, volume: float, sigf: np.ndarray, Nt: np.ndarray,
               Q: np.ndarray):
    """function converts reactor power to neutron flux
    

    Parameters
    ----------
    power : float
        total reactor power, in units of Watts.
    volume : float
        burnable region volume, in units of cm3.
    sigf : np.ndarray
        microscopic fission cross section.
    Nt : np.ndarray
        isotopic concentration, in units of atoms/cm3.
    Q : np.ndarray
        energy released per fission, in units of J/fission.

    Returns
    -------
    flux : float
        neutron flux, in units of n/cm2-s.

    """
    flux = power / (sigf * Nt * Q * volume).sum()

    return flux


def crossSections(Tc: float, Pc: float):
    """function obtains cross sections as a function of chamber temperature and
    pressure. 
    

    Parameters
    ----------
    Tc : float
        Chamber temperature, in units of kelvin.
    Pc : float
        Chamber pressure in units of pascals.

    Returns
    -------
    sigc : np.ndarray
        n,gamma microscopic cross sections, in units of cm2.
    sigf : np.ndarray
        fission microscopic cross sections, in units of cm2.
    mtxFY : np.ndarray
        fission yield matrix.
    Q : np.ndarray
        energy released per fission.

    """

    sigc = np.array([6.8, 250537.62, 132.47, 6968.75, 5.0, 8.0]) * 1e-24
    sigf = np.array([0.0, 0.0, 0.0, 0.0, 97., 3.8]) * 1e-24

    mtxFY = np.array([
        [0.0000, 0.0000, 0.0000, 0.0000, 0.06306, 0.06306],  # 531350
        [0.0000, 0.0000, 0.0000, 0.0000, 0.00248, 0.00248],  # 541350
        [0.0000, 0.0000, 0.0000, 0.0000, 0.01100, 0.01100],  # 611490
        [0.0000, 0.0000, 0.0000, 0.0000, 0.00000, 0],  # 621490
        [0.0000, 0.0000, 0.0000, 0.0000, 0.00000, 0],  # 922350
        [0.0000, 0.0000, 0.0000, 0.0000, 0.00000, 0]])  # 922380

    Q = 1.60218e-13 * np.array([0.0, 0.0, 0.0, 0.0, 202.44, 202.44])  # J/fission

    return sigc, sigf, mtxFY, Q


def exampleSolver(n0, power, volume, Tc, Pc, timepoints, rtol=1E-9):
    """function is an example solution for testing"""
    solution = odeint(dndt, tuple(n0), timepoints,
                      args=(power, volume, Tc, Pc), rtol=rtol)

    return solution
