# -*- coding: utf-8 -*-
"""
@author: Matt Krecicki


function for plotting

"""

import numpy as np
import matplotlib.pyplot as plt
from kinetics.errors.checkerrors import _isstr, _isbool, _inlist, _isint, \
    _isequallength, _isarray, _islist


def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
        

def _checktripleaxis(x, y1, y2, y3, label1, label2, label3, xlabel, ylabel1, 
                     ylabel2, ylabel3, xscale, yscale1, yscale2, yscale3, loc,
                     c1, c2, c3, marker1, fill1, alpha1, linestyle1, marker2,
                     fill2, alpha2, linestyle2, marker3, fill3, alpha3, 
                     linestyle3, markeredgewidth, save, dpi, ylim1, ylim2,
                     ylim3):
    
    _isarray(x, "x-axis values")
    _isarray(y1, "left-hand side axis values")
    for i in range(len(y1)): _isarray(y1[i], "left-hand side axis values")
    for i in range(len(y1)): _isequallength(y1[i], len(x), "left-hand side axis values")
    _isarray(y2, "right-hand side axis values")
    for i in range(len(y2)): _isarray(y2[i], "right-hand side axis values")
    for i in range(len(y2)): _isequallength(y2[i], len(x), "right-hand side axis values")
    _isarray(y3, "second right hand side axis values")
    for i in range(len(y3)): _isequallength(y3[i], len(x), "second right-hand side axis values")


def _fillintripleaxis(y1, y2, y3, c1, c2, c3, marker1, fill1, alpha1, 
                      linestyle1, marker2, fill2, alpha2, linestyle2, marker3,
                      fill3, alpha3, linestyle3):
    
    if not c1:
        c1 = []
        for i in range(len(y1)): c1.append("C"+str(i))
            
    if not c2:
        c2 = []
        for i in range(len(y2)): c2.append("C"+str(i+len(c1)))
    
    if not c3:
        c3 = []
        for i in range(len(y3)): c3.append("C"+str(i+len(c1)+len(c2)))
    
    if not linestyle1:
        linestyle1 = ["-"] * len(y1)
    if not linestyle2:
        linestyle2 = ["-"] * len(y2)
    if not linestyle3:
        linestyle3 = ["-"] * len(y3)
    
    if not marker1:
        marker1 = ["o"] * len(y1)
    if not marker2:
        marker2 = ["x"] * len(y2)
    if not marker3:
        marker3 = ["v"] * len(y3)
    
    if not fill1:
        fill1 = ["none"] * len(y1)
    if not fill2:
        fill2 = ["none"] * len(y2)
    if not fill3:
        fill3 = ["none"] * len(y3)
        
    if not alpha1:
        alpha1 = [1.0] * len(y1)
    if not alpha2:
        alpha2 = [1.0] * len(y2)
    if not alpha3:
        alpha3 = [1.0] * len(y3)
    
    return c1, c2, c3, linestyle1, linestyle2, linestyle3, marker1, marker2,\
        marker3, fill1, fill2, fill3, alpha1, alpha2, alpha3
        
    
def tripleaxis(x, y1, y2, y3, label1, label2, label3, xlabel, ylabel1, ylabel2,
               ylabel3, xscale='linear', yscale1='linear', yscale2='linear',
               yscale3='linear', loc='center right', c1=False, c2=False,
               c3=False, marker1=False, fill1=False, alpha1=False,
               linestyle1=False, marker2=False, fill2=False, alpha2=False,
               linestyle2=False, marker3=False, fill3=False, alpha3=False,
               linestyle3=False, markeredgewidth=0.75, save=False, dpi=1800,
               ylim1=False, ylim2=False, ylim3=False):

    _checktripleaxis(x, y1, y2, y3, label1, label2, label3, xlabel, ylabel1,\
        ylabel2, ylabel3, xscale, yscale1, yscale2, yscale3, loc, c1, c2, c3,\
            marker1, fill1, alpha1, linestyle1, marker2, fill2, alpha2,\
                linestyle2, marker3, fill3, alpha3, linestyle3,\
                    markeredgewidth, save, dpi, ylim1, ylim2, ylim3)
    
    c1, c2, c3, linestyle1, linestyle2, linestyle3, marker1, marker2,\
        marker3, fill1, fill2, fill3, alpha1, alpha2, alpha3 = \
            _fillintripleaxis(y1, y2, y3, c1, c2, c3, marker1, fill1, alpha1,\
                              linestyle1, marker2, fill2, alpha2, linestyle2,\
                                  marker3, fill3, alpha3, linestyle3)
    
    fig, host = plt.subplots()
    fig.subplots_adjust(right=0.75)

    par1 = host.twinx()
    par2 = host.twinx()

    par2.spines["right"].set_position(("axes", 1.2))
    make_patch_spines_invisible(par2)
    par2.spines["right"].set_visible(True)
    
    lns = []
    for i in range(len(y1)):
        p = host.plot(x, y1[i], c1[i], label=label1[i], marker=marker1[i], 
                     fillstyle=fill1[i], alpha=alpha1[i],
                     linestyle=linestyle1[i], markeredgewidth=markeredgewidth,
                     c=c1[i])
        lns.append(p)
        
    for i in range(len(y2)):
        p = par1.plot(x, y2[i], c2[i], label=label2[i], marker=marker2[i], 
                     fillstyle=fill2[i], alpha=alpha2[i],
                     linestyle=linestyle2[i], markeredgewidth=markeredgewidth,
                     c=c2[i])
        lns.append(p)

    for i in range(len(y3)):
        p = par2.plot(x, y3[i], c3[i], label=label3[i], marker=marker3[i], 
                     fillstyle=fill3[i], alpha=alpha3[i],
                     linestyle=linestyle3[i], markeredgewidth=markeredgewidth,
                     c=c3[i])
        lns.append(p)

    host.set_xlabel(xlabel)
    host.set_ylabel(ylabel1)
    par1.set_ylabel(ylabel2)
    par2.set_ylabel(ylabel3)
    
    if ylim1:
        host.set_ylim(ylim1[0], ylim1[1])
    if ylim2:
        par1.set_ylim(ylim2[0], ylim2[1])
    if ylim3:
        par2.set_ylim(ylim3[0], ylim3[1])

    ls = lns[0]
    for i in range(len(lns)-1):
        ls = ls + lns[i+1]
    labs = [l.get_label() for l in ls]
    host.legend(ls, labs, loc=loc)
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=dpi)
    plt.show()
    plt.close()
    

def _correctdualaxis(L1, L2, c1, c2, marker1, fill1, alpha1, linestyle1,
                     marker2, fill2, alpha2, linestyle2):
    
    if type(alpha1) == bool:
        alpha1 = [1.0] * L1
    if type(alpha2) == bool:
        alpha2 = [1.0] * L2
    if type(fill1) == bool:
        fill1 = ['none'] * L1
    if type(fill2) == bool:
        fill2 = ['none'] * L2
    if type(linestyle1) == bool:
        linestyle1 = ['-'] * L1
    if type(linestyle2) == bool:
        linestyle2 = ['-'] * L2
    if type(marker1) == bool:
        marker1 = ['o'] * L1
    if type(marker2) == bool:
        marker2 = ['o'] * L2
    if type(c1) == bool:
        c1 = []
        for i in range(L1):
            c1.append("C"+str(i))
    if type(c2) == bool:
        c2 = []
        for i in range(L2):
            c2.append("C"+str((i+L1)))
    
    return c1, c2, marker1, fill1, alpha1, linestyle1, marker2, fill2, alpha2,\
        linestyle2


def dualaxis(x, y1, y2, label1, label2, xlabel, ylabel1, ylabel2, xscale='linear',
             yscale1='linear', yscale2='linear', loc='center right', c1=False,
             c2=False, marker1=False, fill1=False, alpha1=False, linestyle1=False,
             marker2=False, fill2=False, alpha2=False, linestyle2=False,
             markeredgewidth=1.0, linewidth=0.25, save=False, dpi=1800,
             ylim1=None, ylim2=None):
    
    c1, c2, marker1, fill1, alpha1, linestyle1, marker2, fill2, alpha2,\
        linestyle2 = _correctdualaxis(len(y1), len(y2), c1, c2, marker1, fill1,\
            alpha1, linestyle1, marker2, fill2, alpha2, linestyle2)
            
    fig, ax = plt.subplots()
    lns = []
    
    for i in range(len(y1)):
        li = ax.plot(x, y1[i], label=label1[i], marker=marker1[i], 
                     fillstyle=fill1[i], alpha=alpha1[i], linestyle=linestyle1[i],
                     markeredgewidth=markeredgewidth, c=c1[i], lw=linewidth)
        lns.append(li)
    
    ax2=ax.twinx()
    for i in range(len(y2)):    
        li = ax2.plot(x, y2[i], label=label2[i], marker=marker2[i], 
                     fillstyle=fill2[i], alpha=alpha2[i], linestyle=linestyle2[i],
                     markeredgewidth=markeredgewidth, c=c2[i], lw=linewidth)
        lns.append(li)

    ax.set_yscale(yscale1)
    ax2.set_yscale(yscale2)
    
    ax.set_xscale(xscale)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel1)
    
    ax2.set_ylabel(ylabel2)
    
    if ylim1 is not None:
        ax.set_ylim(ylim1)
    if ylim2 is not None:
        ax2.set_ylim(ylim2)
    
    ls = lns[0]
    for i in range(len(lns)-1):
        ls = ls + lns[i+1]
    labs = [l.get_label() for l in ls]
    ax.legend(ls, labs, loc=loc)
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=dpi)
    plt.show()
    plt.close()
        

def _checklineplotinputs(x, y, label, xlabel, ylabel, title, grid, save, dpi,
                         markers, fills, linestyles, colors, ylim):
    #check arrays
    _isarray(x, 'x-axis values')
    _isarray(y, 'y-axis values')
    
    #compare array lengths
    if len(x) != len(y):
        raise ValueError('x-axis and y-axis lengths not equal')
    for i in range(len(x)):
        if len(x[i]) != len(y[i]):
            raise ValueError('{} index entries are not equal'.format(i))
    
    #check strings
    if xlabel:
        _isstr(xlabel, 'x axis label')
    if ylabel:
        _isstr(ylabel, 'y axis label')
    if title:
        _isstr(title, 'title')

    #check bool
    _isbool(grid, 'grid flag')

    if save:
        _isstr(save, 'plot save name')
        if save[-4:] != '.png':
            raise ValueError('figure file name must end in .png not: {}'.format(save[-4:]))
    
    if markers:
        _isarray(markers, "line marker types")
        _isequallength(markers, len(x), "line marker types")
        for i in markers: _isstr(i, "marker type")
    
    if fills:
        _isarray(fills, "marker fill types")
        _isequallength(fills, len(x), "marker fill types")
        for i in fills: _isstr(i, "fill type")
    
    if linestyles:
        _isarray(linestyles, "line styles")
        _isequallength(linestyles, len(x), "line marker types")
        for i in linestyles: _isstr(i, "line style")
    
    if colors:
        _isarray(colors, "line colors")
        _isequallength(colors, len(x), "colors")
        for i in colors: _isstr(i, "color")        
    
    if ylim:
        _islist(ylim, "y-axis limits")
        _isequallength(ylim, 2, "y-axis limits")
        if ylim[1] <= ylim[0]:
            raise ValueError("lower ylimit is greater than upper limit")
    
        
def lineplot(x, y, label=False, xlabel=False, ylabel=False, title=False, grid=False,
             save=False, dpi=1600, markersize=1.5, markers=False, fills=False,
             show=True, linestyles=False, yscale=False, xscale=False, colors=False,
             ySize=10, xSize=10, figuresize=None, linewidth=1.0, ylim=False,
             loc="best", xticks=False):
    '''function plots simple line plot.
    

    Parameters
    ----------
    x : array-like 
        array or list containing x-axis value arrays. 
    y : array-like
        array or list containing y-axis value arrays. 
    label : list, optional
        list og labels for each curve. The default is False.
    xlabel : str, optional
        x-axis label. The default is False.
    ylabel : str, optional
        y-axis label. The default is False.
    title : str, optional
        figure title. The default is False.
    grid : bool, optional
        flag to plot grid lines on figure. The default is False.
    save : str, optional
        figure save name, if given figure will be saved. The default is False.
    dpi : int, optional
        save image quality, lowest recommended setting is 1600. The default is
        1600.
    markersize : float, optional
        size of marker for each data point. The default is 1.5.
    marker : str, optional
        data point marker style. The default is True.
    show : bool, optional
        flag to show plot after creation. The default is True.
    colors : list
        list of colors for each line plot
    ySize : int
        text size of y-axis label
    xSize : int
        text size of x-axis label
    ylim : list
        list of lower and upper y-axis limits

    Returns
    -------
    None.

    '''
    _checklineplotinputs(x, y, label, xlabel, ylabel, title, grid, save, dpi,
                         markers, fills, linestyles, colors, ylim)
    if not colors:
        colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9',
                  'C10', 'C11', 'C12', 'C13']
    
    if not markers:
        markers =['o'] * len(x)
    if not fills:
        fills =['none'] * len(x)
    if not linestyles:
        linestyles = ["-"] * len(x)
    
    for i in range(len(x)):
        if label:
            plt.plot(x[i], y[i], marker=markers[i], fillstyle=fills[i], 
                     markeredgewidth=markersize, label=label[i], c=colors[i],
                     linestyle=linestyles[i], linewidth=linewidth)
        else:
            plt.plot(x[i], y[i], marker=markers[i], fillstyle=fills[i], 
                     markeredgewidth=markersize, c=colors[i], 
                     linestyle=linestyles[i], linewidth=linewidth)
    if label:
        plt.legend(loc=loc)
    if title:
        plt.title(title)
    if xlabel:
        plt.xlabel(xlabel, fontsize=xSize)
    if ylabel:
        plt.ylabel(ylabel, fontsize=ySize)
    if grid:
        plt.grid()
    if xscale:
        plt.xscale(xscale)
    if yscale:
        plt.yscale(yscale)
    if ylim:
        plt.ylim(ylim[0], ylim[1])
    
    if xticks:
        plt.xticks(xticks)
    
    if figuresize is not None:
        plt.figure(figsize=figuresize)
    
    plt.tight_layout()
    if save:
        _isstr(save, 'figure save name')
        plt.savefig(save, dpi=dpi)
    if show:
        plt.show()
    plt.close()    
    

def _check2Dinputs(x, y, z, xlabel, ylabel, zlabel, title, grid, save, cmap,
                   dpi):
    '''function checks inputs for 2D heat map plotter
    
    
    Parameters
    ----------
    x : np.array, 2D
        x axis values.
    y : np.array, 2D
        y axis values.
    z : np.array, 2D
        heat map data values.
    xlabel : str
        x axis label.
    ylabel : str
        y axis label.
    zlabel : str
        color bar label.
    title : str
        plot title.
    grid : bool
        figure grids.
    save : str
        figure save name.
    camp : str
        color map for plotting.
    dpi : int
        fieldity of saved figure.

    Returns
    -------
    None.

    '''
    
    if title:
        _isstr(title, 'plot title')
    if xlabel:
        _isstr(xlabel, 'x-axis label')
    if ylabel:
        _isstr(ylabel, 'y-axis label')
    if zlabel:
        _isstr(zlabel, 'z-axis label')
    _isbool(grid, 'grids')
    cmaps = ['viridis', 'RdBu', 'jet', 'plasma', 'magma', 'cool']
    _inlist(cmap, 'acceptable color schemes', cmaps)
    if save:
        _isstr(save, 'plot save name')
        if save[-4:] != '.png':
            raise ValueError('figure file name must end in .png not: {}'.format(save[-4:]))
        
    _isint(dpi, 'plot resolution')
    

def plot2Dfield(x, y, z, xlabel=False, ylabel=False, zlabel=False, 
                title=False, grid=True, save=False, cmap='viridis', dpi=1600,
                zlimits=False, Format='%.0f', xIntrpGrid=False,
                yIntrpGrid=False):
    '''function plots 2D field as a heatmap
    

    Parameters
    ----------
    x : np.ndarray, 2D
        x axis position value.
    y : np.ndarray, 2D
        y axis position value.
    z : np.ndarray, 2D
        value for plotting.
    xlabel : str, optional
        x-axis label. The default is False.
    ylabel : str, optional
        y-label. The default is False.
    zlabel : str, optional
        z-label. The default is False.
    title : str, optional
        figure title. The default is False.
    grid : bool, optional
        flag for plotting grid lines. The default is True.
    save : bool, optional
        name for figure save file. The default is False.
    cmap : str, optional
        color scheme. The default is 'viridis'.
    dpi : int, optional
        save figure resolution. The default is 1600.

    Returns
    -------
    None.

    '''
    _check2Dinputs(x, y, z, xlabel, ylabel, zlabel, title, grid, save, cmap, 
                   dpi)
    
    fig, ax = plt.subplots()
    if zlimits:
        vmax = zlimits[1]
        vmin = zlimits[0]
    else:
        vmax = np.max(z)
        vmin = np.min(z)
    
    c = ax.pcolormesh(x, y, z, cmap=cmap, vmin=vmin, vmax=vmax)
    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
        
    # set the limits of the plot to the limits of the data
    ax.axis([np.min(x), np.max(x), np.min(y), np.max(y)])
    
    clb = fig.colorbar(c, ax=ax, format=Format)
    if zlabel:
        clb.set_label(zlabel)
    fig.tight_layout()
    if save:
        plt.savefig(save, dpi=dpi)
    plt.show()
    plt.close()


def scatterplot(x, y, label=False, xlabel=False, ylabel=False, title=False, grid=False,
             save=False, dpi=1600, markersize=1.5, markers=False, fills=False,
             show=True, yscale=False, xscale=False, colors=False,
             ySize=10, xSize=10, figuresize=None, linewidth=1.0, ylim=False,
             loc="best", xticks=False):
    '''function plots simple line plot.
    

    Parameters
    ----------
    x : array-like 
        array or list containing x-axis value arrays. 
    y : array-like
        array or list containing y-axis value arrays. 
    label : list, optional
        list og labels for each curve. The default is False.
    xlabel : str, optional
        x-axis label. The default is False.
    ylabel : str, optional
        y-axis label. The default is False.
    title : str, optional
        figure title. The default is False.
    grid : bool, optional
        flag to plot grid lines on figure. The default is False.
    save : str, optional
        figure save name, if given figure will be saved. The default is False.
    dpi : int, optional
        save image quality, lowest recommended setting is 1600. The default is
        1600.
    markersize : float, optional
        size of marker for each data point. The default is 1.5.
    marker : str, optional
        data point marker style. The default is True.
    show : bool, optional
        flag to show plot after creation. The default is True.
    colors : list
        list of colors for each line plot
    ySize : int
        text size of y-axis label
    xSize : int
        text size of x-axis label
    ylim : list
        list of lower and upper y-axis limits

    Returns
    -------
    None.

    '''
    _checklineplotinputs(x, y, label, xlabel, ylabel, title, grid, save, dpi,
                         markers, fills, False, colors, ylim)
    if not colors:
        colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9',
                  'C10', 'C11', 'C12', 'C13']
    
    if not markers:
        markers =['o'] * len(x)
    if not fills:
        fills =['none'] * len(x)
    
    for i in range(len(x)):
        if label:
            plt.scatter(x[i], y[i], marker=markers[i], 
                     label=label[i], c=colors[i],
                     )
        else:
            plt.scatter(x[i], y[i], marker=markers[i], 
                      c=colors[i], 
                     )
    if label:
        plt.legend(loc=loc)
    if title:
        plt.title(title)
    if xlabel:
        plt.xlabel(xlabel, fontsize=xSize)
    if ylabel:
        plt.ylabel(ylabel, fontsize=ySize)
    if grid:
        plt.grid()
    if xscale:
        plt.xscale(xscale)
    if yscale:
        plt.yscale(yscale)
    if ylim:
        plt.ylim(ylim[0], ylim[1])
    
    if xticks:
        plt.xticks(xticks)
    
    if figuresize is not None:
        plt.figure(figsize=figuresize)
    
    plt.tight_layout()
    if save:
        _isstr(save, 'figure save name')
        plt.savefig(save, dpi=dpi)
    if show:
        plt.show()
    plt.close()   
    

def _checkStackedSubLineplot(noPlots, x, y, xlabel, ylabel, markersize,
                             linewidth, hspace, linestyles, markers, fills, colors,
                             grids, show, save, dpi, ylim, loc, xscale, yscale):
    """simple error checking"""
    
    _isint(noPlots, "number of stacked plots")
    

def _correctStackedSubLineplot(noPlots, linestyles, markers, fills, colors,
                               grids, show, save, ylim, loc, xscale, yscale):
    
    if fills is False:
        fills = ['none'] * noPlots
    
    if colors is False:
        colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9',
                  'C10', 'C11', 'C12', 'C13']
    
    
    return linestyles, markers, fills, colors, grids, show, save, ylim, loc,\
        xscale, yscale
    
    
def stackedSubLineplot(noPlots, x, y, xlabel, ylabel, markersize=1.5, margin=False,
                       linewidth=1.0, hspace=0.5, linestyles=False, markers=False,
                       fills=False, colors=False, labels=False, grids=False, 
                       show=True, save=False, dpi=1600, ylims=False, locs=False,
                       xscale=False, yscale=False, fontsizes=False, alphas=False,
                       figsize=(5, 5), xTextSize=10, yTextSize=False,
                       xticksize=False, yticksize=False):
    
    #check for inputs errors---------------------------------------------------
    _checkStackedSubLineplot(noPlots, x, y, xlabel, ylabel, markersize,\
        linewidth, hspace, linestyles, markers, fills, colors, grids, show, save, dpi,
        ylims, locs, xscale, yscale)
    
    #correct not specified inputs----------------------------------------------
    linestyles, markers, fills, colors, grids, show, save, ylim, loc, xscale,\
        yscale = _correctStackedSubLineplot(noPlots, linestyles, markers,\
            fills, colors, grids, show, save, ylims, locs, xscale, yscale)
            
    #create plot---------------------------------------------------------------
    
    fig, axs = plt.subplots(noPlots, 1, sharex=True)
    
    fig.set_figwidth(figsize[0])
    fig.set_figheight(figsize[1])
    
    
    fig.subplots_adjust(hspace=hspace)
    
    for i in range(noPlots):
        
        axi = axs[i]
        
        for j in range(len(y[i])):
            axi.plot(x[i][j], y[i][j], color=colors[i][j], marker=markers[i][j],
                     linestyle=linestyles[i][j], fillstyle=fills[i][j], 
                     label=labels[i][j], markeredgewidth=markersize,
                     linewidth=linewidth, alpha=alphas[i][j])
                
        #apply grid 
        if grids[i]: axi.grid(True)
        
        #apply ylabel
        axi.set_ylabel(ylabel[i], fontsize=fontsizes[i])
        
        #apply legend to each graph
        axi.legend(loc=locs[i])
    
        #size y-axis 
        if ylim[i] is False:
            mini = np.min(y[i]) * (1 - margin[i])
            maxi = np.max(y[i]) * (1 + margin[i])
            axi.set_ylim((mini, maxi))
        else:
            axi.set_ylim(ylim[i])
        
        #set x/y-scale 
        if xscale[i]: axi.set_xscale(xscale[i])
        if yscale[i]: axi.set_yscale(yscale[i])
        
        #set y-label tick size
        if yticksize:
            axi.tick_params(axis='y', labelsize=yticksize[i])
    
    axi.set_xlabel(xlabel, fontsize=xTextSize)
    
    if xticksize:
        axi.tick_params(axis='x', labelsize=xticksize)
    
    fig.tight_layout()
    
    #save figure to file
    if save:
        _isstr(save, 'figure save name')
        plt.savefig(save, dpi=dpi)
    
    #show plot
    if show: plt.show()
    
    #clear memory of plot
    plt.close()   
    







