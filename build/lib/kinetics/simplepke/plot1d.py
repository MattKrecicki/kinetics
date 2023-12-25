# -*- coding: utf-8 -*-
"""plot1d

Basic plotting capability.

Created on Tue July 19 11:10:00 2022 @author: Dan Kotlyar
Last updated om Tue July 19 11:30:00 2022 @author: Dan Kotlyar
email: dan.kotlyar@me.gatech.edu
"""

import matplotlib.pyplot as plt

# Default values
numSteps = 100  # number of time-steps
FONT_SIZE = 16  # font size for plotting purposes




def plot1D(xvals, yvals, xlabel=None, ylabel=None, fontsize=FONT_SIZE,
           marker="--*", markerfill=False, markersize=6):
    """Plot the 1D slab neutron flux distribution.

    The function is meant to be executed after the one-group slab diffusion
    solver is ran.

    Parameters
    ----------
    power : bool, optional
        if ``power`` is `True` or not included the
        power distribution is plotted, otherwise not
    xvals : ndarray
        x-axis values
    yvals : ndarray
        y-axis values
    xlabel : str
        x-axis label with a default ``Length, meters``
    ylabel : str
        y-axis label with a default ``Normalized Flux``
    fontsize : float
        font size value
    markers : str or list of strings
        markers type
    markerfill : bool
        True if the marking filling to be included and False otherwise
    markersize : int or float
        size of the marker with a default of 8.

    """

    if xlabel is None:
        xlabel = "Slab length, meters"
    if ylabel is None:
        ylabel = "Normalized flux"

    if markerfill:
        mfc = "white"  # marker fill color
    else:
        mfc = None

    plt.ylabel(ylabel)
    plt.xlabel(xlabel)

    plt.plot(xvals, yvals, marker, mfc=mfc, ms=markersize)
    plt.grid()
    plt.rc('font', size=fontsize)      # text sizes
    plt.rc('axes', labelsize=fontsize)  # labels
    plt.rc('xtick', labelsize=fontsize)  # tick labels
    plt.rc('ytick', labelsize=fontsize)  # tick labels
