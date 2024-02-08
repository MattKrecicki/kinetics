# -*- coding: utf-8 -*-
"""checkerrors

Set of muted functions to examine various errors.

The error checking provides the ability to monitor whether variables
are of a certain type, within a certain values range, or part of a list.

@author: matt krecicki
@email: matthewkrecicki@gmail.com

"""

import numbers
import numpy as np


def _isnumber(var, description):
    """checks if the variable is a scalar"""
    if not isinstance(var, numbers.Real):
        raise TypeError("{} must be a scalar and not {}".
                        format(description, var))


def _isint(var, description):
    """checks if the variable is an integer"""
    if not isinstance(var, numbers.Integral):
        raise TypeError("{} must be an integer and not {}".
                        format(description, var))


def _isstr(var, description):
    """checks if the variable of string type"""
    if not isinstance(var, str):
        raise TypeError("{} must be string and not {}"
                        .format(description, var))


def _isbool(var, description):
    """checks if the variable of boolean type"""
    if not isinstance(var, bool):
        raise TypeError("{} must be bool and not {}"
                        .format(description, var))


def _islist(var, description):
    """checks if the variable of list type"""
    if not isinstance(var, list):
        raise TypeError("{} must be list and not {}"
                        .format(description, var))


def _is1dlist(var, description):
    """checks if the variable of list type"""
    if not isinstance(var, list):
        raise TypeError("{} must be list and not {}"
                        .format(description, var))
    for ivar in var:
        if isinstance(ivar, list):
            raise TypeError("{} must be a 1d list and not {}"
                            .format(description, var))


def _is2dlist(var, description):
    """checks if the variable of list type"""
    if not isinstance(var, list):
        raise TypeError("{} must be list and not {}"
                        .format(description, var))
    for ivar in var:
        if not isinstance(ivar, list):
            raise TypeError("{} must be a 2d list and not {}"
                            .format(description, var))
        for jvar in ivar:
            if isinstance(jvar, list):
                raise TypeError("{} must be a 2d list and not {}"
                                .format(description, var))


def _isobject(var, objtype, description):
    """checks if the variable of a specific object type"""
    if not isinstance(var, objtype):
        raise TypeError("{} must be {} and not {}"
                        .format(description, objtype, var))


def _isarray(var, description):
    """checks if the variable of an array type"""
    if not isinstance(var, (np.ndarray, list)):
        raise TypeError("{} must be an array and not {}"
                        .format(description, var))


def _isndarray(var, description):
    """checks if the variable of an ndarray type"""
    if not isinstance(var, np.ndarray):
        raise TypeError("{} must be a ndarray and not {}"
                        .format(description, var))


def _ispositive(var, description):
    """checks if the variable is positive"""
    _isnumber(var, description)
    if not var > 0:
        raise ValueError("{} must be positive and not {}"
                         .format(description, var))


def _isnegative(var, description):
    """checks if the variable is negative"""
    _isnumber(var, description)
    if not var < 0:
        raise ValueError("{} must be negative and not {}"
                         .format(description, var))


def _iszeropositive(var, description):
    """checks if the variable is zero or positive"""
    _isnumber(var, description)
    if not var >= 0:
        raise ValueError("{} must be zero or positive and not {}"
                         .format(description, var))


def _isnonnegative(var, description):
    """checks if the variable is positive"""
    _isnumber(var, description)
    if not var >= 0:
        raise ValueError("{} must be non-negative and not {}"
                         .format(description, var))


def _ispositiveArray(var, description):
    """checks if the variable is positive"""
    _isndarray(var, description)
    if not (var > 0).all():
        raise ValueError("{} must be positive and not {}"
                         .format(description, var))


def _isnonNegativeArray(var, description):
    """checks if the variable is positive"""
    _isndarray(var, description)
    if not (var >= 0).all():
        raise ValueError("{} must be positive and not {}"
                         .format(description, var))


def _isBoundArray(var, bounds, description):
    """checks if the variable is positive"""
    _isndarray(var, description)
    if (var > bounds[1]).any() or (var < bounds[0]).any():
        raise ValueError("{} must be bounded {} by [{}]"
                         .format(description, var, bounds))


def _is1darray(var, description):
    """checks if the array is 1D"""
    if np.array(var).ndim != 1 or len(var) == 0:
        raise TypeError("{} must be 1D array and not {}"
                        .format(description, var))


def _is2darray(var, description):
    """checks if the array is 2D"""
    _isndarray(var, description)
    if np.array(var).ndim != 2:
        raise TypeError("{} must be 2D array and not {}"
                        .format(description, var))


def _isndimarray(var, description, ndim):
    """checks if the array is n-Dim"""
    _isndarray(var, description)
    if np.array(var).ndim != ndim:
        raise TypeError("{} must be {}-dim array and not {}"
                        .format(description, ndim, var))


def _exp2dshape(var, expshape, description):
    """checks if the 2d ndarray has a certain shape"""
    _is2darray(var, description)
    if np.array(var).shape != expshape:
        raise ValueError("{} must have a shape of {} and not {}"
                         .format(description, expshape, var))


def _allzero(var, description):
    """checks if the vector is all zero"""
    _isndarray(var, description)
    if (var == 0.0).all():
        raise ValueError("{} cannot have only zero values {}"
                         .format(description, var))


def _anyzero(var, description):
    """checks if the vector contains any zeroes"""
    _isndarray(var, description)
    if (var == 0.0).any():
        raise ValueError("{} cannot have any zero values {}"
                         .format(description, var))


def _anynegative(var, description):
    """checks if the vector contains any zeroes"""
    _isndarray(var, description)
    if (var < 0.0).any():
        raise ValueError("{} cannot have any negative values {}"
                         .format(description, var))


def _isequallength(var, expLen, description):
    """checks if the length of the 1D array equal to the an expected value"""
    if len(var) != expLen:
        raise ValueError("{} must have {} components and not {}"
                         .format(description, expLen, var))


def _inrange(var, description, limits, upBound=True, lowBound=True):
    """checks if the variable is in a certain values range"""
    _isnumber(var, description)
    if upBound and lowBound:
        if not limits[0] <= var <= limits[1]:
            raise ValueError("{} must be in the [{},{}] range and not {}"
                             .format(description, limits[0], limits[1], var))
    elif upBound:
        if not limits[0] < var <= limits[1]:
            raise ValueError("{} must be in the ({},{}] range and not {}"
                             .format(description, limits[0], limits[1], var))
    elif lowBound:
        if not limits[0] <= var < limits[1]:
            raise ValueError("{} must be in the [{},{}) range and not {}"
                             .format(description, limits[0], limits[1], var))
    else:
        if not limits[0] < var < limits[1]:
            raise ValueError("{} must be in the ({},{}) range and not {}"
                             .format(description, limits[0], limits[1], var))


def _inlist(var, description, keyslist):
    """checks if the variable is in a certain list"""
    if var not in keyslist:
        raise KeyError("{}=<{}> does not exist!!! in the option list {}"
                       .format(description, var, keyslist))
        
        
def _innotlist(var, description, keyslist):
    """checks if the variable is not in a certain list"""
    if var in keyslist:
        raise KeyError("{}=<{}> exists!!! in the option list {}"
                       .format(description, var, keyslist))


def _compare2lists(list1, list2, description1, description2):
    """sort and compare two lists"""
    _islist(list1, description1)
    _islist(list2, description2)
    list1.sort()
    list2.sort()
    if list1 != list2:
        raise ValueError("{} {}\n is not equal to \n {} {}"
                         .format(description1, list1, description2, list2))


def _isuniquelist(list0, description):
    """check that components in a list do not appear more than once"""
    # insert the list to the set
    list_set = set(list0)
    # convert the set to the list
    unique_list = (list(list_set))
    if len(unique_list) != len(list0):
        raise ValueError("{} {} contains duplicate components".format(
            description, list0))


def _issortedarray(var, description):
    """check that an array is sorted"""
    if not (np.diff(var) >= 0).all():
        raise ValueError("{} must be sorted {}"
                         .format(description, var))


def _arriscloseInList(myarr, description, list_arrays):
    """test for approximate equality (for floating point types)"""
    idx0 = -1
    for idx, elem in enumerate(list_arrays):
        if (elem.size == myarr.size) and np.allclose(elem, myarr):
            idx0 = idx
            break
    if idx0 > -1:
        return idx0
    else:
        raise ValueError("{} {} does not exist in {}"
                         .format(description, myarr, list_arrays))


def _arrIdenticalInList(myarr, description, list_arrays):
    """test for identical equality"""
    idx0 = -1
    for idx, elem in enumerate(list_arrays):
        if (elem.size == myarr.size) and np.array_equal(elem, myarr):
            idx0 = idx
            break
    if idx0 > -1:
        return idx0
    else:
        raise ValueError("{} {} does not exist in {}"
                         .format(description, myarr, list_arrays))


def _isSortedArray(var, description):
    """checks if the array is sorted and unique"""
    _isndarray(var, description)
    if len(var) > 1:
        newvar = var[1:] - var[0:-1]
    else:
        newvar = var
    if not (newvar > 0).all():
        raise ValueError("{} must be sorted and unique {}"
                         .format(description, var))