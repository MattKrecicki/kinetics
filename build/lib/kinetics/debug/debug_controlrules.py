# -*- coding: utf-8 -*-
"""debug_controlrules.py

Define rules to control the dynamic behavior during startups.
The control rules defined here are only meant to control chamber
pressure and temperature.


"""

import numpy as np
from kinetics.functions.control import generalControlRule

# Test two arbitrary ramp intervals
# -----------------------------------------------------------------------------

# linear function
# 2x + 1
# integral of the function is:
# x^2 + x    

# exponential function
# 3*exp(3x)+1
# integral of the function is:
# exp(3x)+ x 

cntlRules = generalControlRule(funcTypes = ['linear', 'exponential'], 
                            coeffs = [[2.0, 1.0], [0.3, 0.3, 1.0]],
                            tends = [30.0, 60.0])

# -----------------------------------------------------------------------------
#                TEST THE VALUES FOR THE FUNCTION
# -----------------------------------------------------------------------------

# numerically solve
nvals = 5
t0 = np.linspace(5.0, 25.0, nvals)
ypred0 = np.zeros(nvals)
for idx, tval in enumerate(t0):
    ypred0[idx] = cntlRules.evaluate(tval)
# analytic solution
yexp0 = 2*t0 + 1

mean_erry0 = np.mean(100*(1 - ypred0/yexp0))
print("Mean error 1st interval = {} %".format(mean_erry0))

t1 = np.linspace(35.0, 55.0, nvals)
ypred1 = np.zeros(nvals)
for idx, tval in enumerate(t1):
    ypred1[idx] = cntlRules.evaluate(tval)
# analytic solution
yexp1 = 0.3*np.exp(0.3*t1)+1

mean_erry1 = np.mean(100*(1 - ypred1/yexp1))
print("Mean error 2nd interval = {} %".format(mean_erry1))


# -----------------------------------------------------------------------------
#                TEST THE VALUES FOR THE INTEGRAL OF THE FUNCTIONS
# -----------------------------------------------------------------------------

# numerically solve for the 1st interval
t0 = 1.0
t1 = 29.0
ypred0 = cntlRules.evaluateIntegral(Tstart=t0, Tend=t1)
# analytic solution
yexp0 = (t1**2 + t1) - (t0**2 + t0)

mean_erry0 = np.mean(100*(1 - ypred0/yexp0))
print("Mean error 1st interval = {} %".format(mean_erry0))

# numerically solve for the 1st and 2nd intervals
t0 = 5.0
t1 = 55.0
ypred1 = cntlRules.evaluateIntegral(Tstart=t0, Tend=t1)
# analytic solution
yexp1 = (30.0**2 + 30.0) - (5.0**2 + 5.0) +\
    np.exp(0.3*55.0)+ 55.0 - (np.exp(0.3*30.0)+ 30.0) 

mean_erry0 = np.mean(100*(1 - ypred1/yexp1))
print("Mean error 2nd interval = {} %".format(mean_erry0))

