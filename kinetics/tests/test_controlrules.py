# -*- coding: utf-8 -*-
"""test_propertytable

Tests to compare that the package produces identical results to
results obtain with the matlab version

Created on Mon Mar  9 19:27:04 2020 @author: Dan Kotlyar
Last updated on Sun April 05 07:30:00 2020 @author: Dan Kotlyar

Script to test propertytable.py function
"""

import pytest

import numpy as np
from ntpSystem.functions.control import generalControlRule


def test_reset():
    """check that reseting works properly"""

    with pytest.raises(TypeError,
                       match="Function*"):
        generalControlRule(funcTypes = 'BAD_VAR_TYPE', 
                           coeffs = [[2.0, 1.0], [3.0, 3.0, 1.0]],
                           tends = [30.0, 60.0, 70.0])

    with pytest.raises(KeyError,
                       match="Function*"):
        generalControlRule(funcTypes = ['linear', 'BAD_FUNCTION'], 
                           coeffs = [[2.0, 1.0], [3.0, 3.0, 1.0]],
                           tends = [30.0, 60.0])

    with pytest.raises(ValueError,
                       match="Intervals*"):
        generalControlRule(funcTypes = ['linear', 'exponential'], 
                           coeffs = [[2.0, 1.0], [3.0, 3.0, 1.0]],
                           tends = [30.0, 60.0, 70.0])

    with pytest.raises(ValueError,
                       match="Function*"):
        generalControlRule(funcTypes = ['linear', 'exponential'], 
                           coeffs = [[2.0], [3.0, 3.0, 1.0]],
                           tends = [30.0, 60.0])

    with pytest.raises(ValueError,
                       match="Intervals*"):
        generalControlRule(funcTypes = ['linear', 'exponential'], 
                           coeffs = [[2.0, 1.0], [3.0, 3.0, 1.0]],
                           tends = [60.0, 30.0])


def test_evaluate():
    """check the evaluate method"""

    cntlRules = generalControlRule(funcTypes = ['linear', 'exponential'], 
                                   coeffs = [[2.0, 1.0], [3.0, 3.0, 1.0]],
                                   tends = [30.0, 60.0])
    
    t1 = 20.0
    ypred1 = cntlRules.evaluate(t1)
    # analytic solution
    yexp1 = 2*t1 + 1
    assert ypred1 == pytest.approx(yexp1)

    t1 = 50.0
    ypred1 = cntlRules.evaluate(t1)
    # analytic solution
    yexp1 = 3*np.exp(3*t1)+1
    assert ypred1 == pytest.approx(yexp1)


def test_evaluateIntegral():
    """check the evaluate integral method"""

    cntlRules = generalControlRule(funcTypes = ['linear', 'exponential'], 
                                   coeffs = [[2.0, 1.0], [0.3, 0.3, 1.0]],
                                   tends = [30.0, 60.0])
    
    t0 = 1.0
    t1 = 29.0
    ypred0 = cntlRules.evaluateIntegral(Tstart=t0, Tend=t1)
    # analytic solution
    yexp0 = (t1**2 + t1) - (t0**2 + t0)
    
    assert ypred0 == pytest.approx(yexp0)

    t0 = 5.0
    t1 = 55.0
    ypred1 = cntlRules.evaluateIntegral(Tstart=t0, Tend=t1)
    # analytic solution
    yexp1 = (30.0**2 + 30.0) - (5.0**2 + 5.0) +\
        np.exp(0.3*55.0)+ 55.0 - (np.exp(0.3*30.0)+ 30.0) 
    assert ypred1 == pytest.approx(yexp1)
