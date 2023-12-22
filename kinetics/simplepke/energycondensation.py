"""energycondensation.py

Energy condensation on multi-group maco- and micro- parameters.
The condensation procedure can be applied on vectors, such as absorption xs,
as well as on scattering matrices.


"""

import numpy as np
from numpy.matlib import repmat

from xsInterface.errors.checkerrors import _isint, _isarray, _is1darray,\
    _is2darray, _inrange, _ispositive, _isequallength, _isnonNegativeArray,\
    _issortedarray


def EnergyCondensation(ndim, ng, boundsE, attr, flux, cutoffE):
    """Energy condensation method from finer-to-coarser energy grid

    Parameters
    ----------
    ndim : int
        number of dimensions. 1 for a vector and 2 for a matrix.
    ng : int
        number of energy groups.
    boundsE : array
        group boundaries in descending order. All boundaries must be included.
    attr : 1-dim or 2-dim array
        macro or micro 1-dim (e.g., absorption) or 2-dim (e.g., scattering) to
        be condensed.
    flux : 1-dim array
        flux/weighting multi-group function
    cutoffE : 1-dim array
        energy cutoffs

    Raises
    ------
    TypeError
        If ``ndim``, ``ng``,  are not integers.
        If ``boundsE``, ``attr``, ``flux``, ``cutoffE`` are not arrays.
    ValueError
        If ``ndim`` is not 1 or 2, ``ng`` is not positive.
        If the length of ``boundsE`` is not ng+1.
        If the length of ``attr`` or ``flux`` is not ng.
        If ``boundsE`` or ``cutoffE`` contain negative values.
        If ``boundsE`` is not sorted in a decreasing order.

    Examples
    --------
    >>> EnergyCondensation(ndim, ng, boundsE, attr, flux, cutoffE)

    """

    # -------------------------------------------------------------------------
    # Variables are properly defined
    # -------------------------------------------------------------------------
    boundsE, attr, flux, cutoffE =\
        _checkErrors(ndim, ng, boundsE, attr, flux, cutoffE)

    # -------------------------------------------------------------------------
    # Find the indices that correpond to the cutoff energies
    # -------------------------------------------------------------------------
    idxE = np.zeros(len(cutoffE), dtype=int)  # max. number of indices
    for idx, valE in enumerate(cutoffE):
        if (valE < boundsE).all():
            raise ValueError(
                "Energy value <{}> in cutoffE is below the lower <{}> bound"
                .format(valE, min(boundsE)))
        if (valE > boundsE).all():
            raise ValueError(
                "Energy value <{}> in cutoffE is above the upper <{}> bound"
                .format(valE, max(boundsE)))
        # energy value is within the bounds
        condIdx = np.argwhere(boundsE <= valE)
        idxE[idx] = condIdx[0, 0]
    idxE = np.append(0, idxE)  # add the zeroth index
    idxE = np.append(idxE, len(boundsE))  # add the last index
    # save only the unique indices
    idxE = np.unique(idxE)

    # -------------------------------------------------------------------------
    # Condense a vector or a scattering matrix
    # -------------------------------------------------------------------------
    condNG = int(len(idxE) - 1)  # number of condensed energy groups
    if ndim == 1:
        condAttr = np.zeros(condNG)  # reset attribute with the condensed grid
        for ig in range(condNG):  # loop over all the condensed groups
            i0 = idxE[ig]  # left bounding index
            i1 = idxE[ig+1]  # right bounding index
            condAttr[ig] = (attr[i0:i1]*flux[i0:i1]).sum()/(flux[i0:i1]).sum()

    idxE[-1] -= 1
    if ndim == 2:  # scattering matrices
        condAttr = np.zeros((condNG, condNG))
        for ig in range(condNG):
            gFromL = idxE[ig]  # left boundary
            gFromR = idxE[ig+1]  # right boundary
            for jg in range(condNG):
                gToL = idxE[jg]
                gToR = idxE[jg+1]
                fluxMtx = repmat(flux[gFromL:gFromR], gToR-gToL, 1)
                fluxMtx = fluxMtx.transpose()
                sctRR = (attr[gFromL:gFromR, gToL:gToR]*fluxMtx).sum()
                condAttr[ig, jg] = sctRR / (flux[gFromL:gFromR]).sum()

    condE = boundsE[idxE]  # energy bounds of the condensed energy structure
    return condAttr, condE


def _checkErrors(ndim, ng, boundsE, attr, flux, cutoffE):
    """check that parameters are properly provided"""
    # ---------------------------------------------------------------------
    # Variables types
    # ---------------------------------------------------------------------
    _isint(ndim, "Dimensions")
    _isint(ng, "Number of energy groups")
    _inrange(ndim, "Dimensions", [1, 2])
    _ispositive(ng, "Number of energy groups")

    _isarray(boundsE, "Energy boundaries")
    _isarray(attr, "Attribute")
    _isarray(flux, "Flux")
    _isarray(cutoffE, "Energy cutoffs")
    # convert to arrays if lists
    boundsE = np.array(boundsE)
    attr = np.array(attr)
    flux = np.array(flux)
    cutoffE = np.array(cutoffE)

    _is1darray(boundsE, "Energy boundaries")
    _is1darray(flux, "Flux")
    _is1darray(cutoffE, "Energy cutoffs")

    # check dimensions of attribute
    if ndim == 1:
        _is1darray(attr, "Attribute")
    else:
        _is2darray(attr, "Attribute")

    # ---------------------------------------------------------------------
    # Variables values
    # ---------------------------------------------------------------------
    _isequallength(boundsE, ng+1, "Energy boundaries")
    _isequallength(attr, ng, "Attribute")
    _isequallength(flux, ng, "Flux")
    _isnonNegativeArray(boundsE, "Energy boundaries")
    _isnonNegativeArray(cutoffE, "Energy cutoffs")
    _issortedarray(boundsE[::-1], "Energy boundaries")

    return boundsE, attr, flux, cutoffE
