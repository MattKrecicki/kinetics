# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 23:41:36 2024

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""

from kinetics.errors.checkerrors import _inrange


def LinearInterp(x0, x1, vals0, vals1, x):
    """Linear interpolation

    Given the dependnecies x0 (left) and x1 (right) a linear interpolation is
    going to be performed on vectors/floats vals0 and vals1 for a specific x.

    Parameters
    ----------
    x0 : float
        a value of the dependnecy (e.g. temperature)
    x1 : float
        a value of the dependnecy
    vals0 : float or array
        values corresponding to a dependnecy `x0`
    vals1 : float or array
        values corresponding to a dependnecy `x1`
    x : float
        the value of the dependency for the interpolation
    extrpFlag : boolean, default True
        extrapolation flag to indicate whether extrapolation is allowed.

    Raises
    ------
    ValueError
        If x is outside the range [x0, x1] and extrapolation flag is False

    Examples
    --------
    >>> LinearInterp(500., 600., 0.01, 0.02, 550, False)
    ... 0.015
    >>> LinearInterp(500., 600., 0.01, 0.02, 490, True)
    ... 0.009

    """

    _inrange(x, "value of x", [x0, x1])
    xd = (x - x0) / (x1 - x0)
    vals = (1 - xd)*vals0 + xd*vals1

    return vals


def BiLinearInterp(x0, x1, y0, y1, vals00, vals10, vals01, vals11, x, y):
    """Bi-Linear interpolation on a x and y grid

    The grid and values are provided in the schematics below:

        (x0, y0)--------(x1, y0)           vals00 -------- vals10
            |   (x,y)   |                    |     vals       |
            |           |         ==>        |                |
        (x0, y1)-------(x1, y1)            vals01 -------- vals11

     (x, y) are interpolated from knowing x0, x1, y0, and y1, and the weights
     are used to obtain the interpolated `vals` from known values (e.g. vals00)
     that correpond to the grid structure.

    Parameters
    ----------
    x0 : float
        a value of the x-dependnecy (e.g. temperature)
    x1 : float
        a value of the x-dependnecy (e.g. temperature)
    y0 : float
        a value of the y-dependnecy (e.g. pressure)
    y1 : float
        a value of the y-dependnecy (e.g. pressure)
    vals00 : float or array
        values corresponding to dependnecy pair (x0, y0)
    vals10 : float or array
        values corresponding to dependnecy pair (x0, y0)
    vals01 : float or array
        values corresponding to dependnecy pair (x0, y0)
    vals11 : float or array
        values corresponding to dependnecy pair (x0, y0)
    x : float
        the value of the x-dependency for the interpolation
    y : float
        the value of the y-dependency for the interpolation
    extrpFlag : boolean, default True
        extrapolation flag to indicate whether extrapolation is allowed.

    Raises
    ------
    ValueError
        If x is outside the range [x0, x1] or y is outside the range [y0, y1]
        and extrapolation flag is False (i.e., extrapolation is not allowed).

    Examples
    --------
    >>> BiLinearInterp(500., 600., 4, 5, 0.01, 0.02, 0.03, 0.04, 550, 4.5)
    ... 0.025
    >>> BiLinearInterp(500., 600., 4, 5, 0.01, 0.02, 0.03, 0.04, 500, 4)
    ... 0.01
    >>> BiLinearInterp(500., 600., 4, 5, 0.01, 0.02, 0.03, 0.04, 600, 5)
    ... 0.04
    >>> BiLinearInterp(500., 600., 4, 5, 0.01, 0.02, 0.03, 0.04, 450, 3.9)
    ... 0.0029999999999999975

    """

    _inrange(x, "value of x", [x0, x1])
    _inrange(y, "value of y", [y0, y1])

    xd = (x - x0) / (x1 - x0)
    yd = (y - y0) / (y1 - y0)

    # Interpolation over the first x-dependency for different y-values
    vals0 = (1 - xd)*vals00 + xd*vals10
    vals1 = (1 - xd)*vals01 + xd*vals11

    # Interpolation for the second y-dep. based on the x-interpolated values
    vals = (1 - yd)*vals0 + yd*vals1

    return vals


def TriLinearInterp(x0, x1, y0, y1, z0, z1, vals000, vals100, vals010, vals110,
                    vals001, vals101, vals011, vals111, x, y, z):
    """Tri-Linear interpolation on (x, y, z) grid

    The grid and values are provided in the schematics below:


  vals011 *------------* vals111
         /|           /|
        / |          / |
vals001* -----------*  |vals101
       |  |         |  |
       |  |         |  |
   010 |  *---------|--* vals110
       | /          | /
       |/           |/
       * -----------*
     vals000        vals100

     (x, y, z) are interpolated from knowing [x0, x1], [y0, y1], [z1, z2]
     and the weights are used to obtain the interpolated `vals`
     from known values (e.g. vals001, vals111, ...)
     that correpond to the parallelepiped/cubic-based grid structure.

    Parameters
    ----------
    x0 : float
        a value of the x-dependnecy (e.g. temperature)
    x1 : float
        a value of the x-dependnecy (e.g. temperature)
    y0 : float
        a value of the y-dependnecy (e.g. pressure)
    y1 : float
        a value of the y-dependnecy (e.g. pressure)
    z0 : float
        a value of the z-dependnecy (e.g. burnup)
    z1 : float
        a value of the z-dependnecy (e.g. burnup)
    vals000 : float or array
        values corresponding to f(x0, y0, z0)
    vals100 : float or array
        values corresponding to f(x1, y0, z0)
    vals010 : float or array
        values corresponding to f(x0, y1, z0)
    vals110 : float or array
        values corresponding to f(x1, y1, z0)
    vals001 : float or array
        values corresponding to f(x0, y0, z1)
    vals101 : float or array
        values corresponding to f(x1, y0, z1)
    vals011 : float or array
        values corresponding to f(x0, y1, z1)
    vals111 : float or array
        values corresponding to f(x1, y1, z1)
    x : float
        the value of the x-dependency for the interpolation
    y : float
        the value of the y-dependency for the interpolation
    z : float
        the value of the z-dependency for the interpolation
    extrpFlag : boolean, default True
        extrapolation flag to indicate whether extrapolation is allowed.

    Raises
    ------
    ValueError
        If x is outside the range [x0, x1] or y is outside the range [y0, y1]
        and extrapolation flag is False (i.e., extrapolation is not allowed).

    Examples
    --------
    >>> TriLinearInterp(500., 600., 4, 5, 0.01, 0.02, 0.03, 0.04, 550, 4.5)
    ... 0.025
    >>> TriLinearInterp(500, 600, 4, 5, 0, 10, 0.01, 0.02, 0.03, 0.04,
                    0.05, 0.06, 0.07, 0.08, 550, 4.5, 5,
                    extrpFlag=True)
    ... 0.045

    """

    _inrange(x, "value of x", [x0, x1])
    _inrange(y, "value of y", [y0, y1])
    _inrange(z, "value of z", [z0, z1])

    xd = (x - x0) / (x1 - x0)
    yd = (y - y0) / (y1 - y0)
    zd = (z - z0) / (z1 - z0)

    # Interpolation over the first x-dependency for diff y and z values
    vals00 = (1 - xd)*vals000 + xd*vals100
    vals01 = (1 - xd)*vals001 + xd*vals101
    vals10 = (1 - xd)*vals010 + xd*vals110
    vals11 = (1 - xd)*vals011 + xd*vals111

    # Interpolation over the y-dep. for diff. z-values
    vals0 = (1 - yd)*vals00 + yd*vals10
    vals1 = (1 - yd)*vals01 + yd*vals11

    # Interpolation for the z-dep.
    vals = (1 - zd)*vals0 + zd*vals1

    return vals