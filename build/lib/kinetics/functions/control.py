"""control.py

Define rules to control the dynamic behavior during startups.
The control rules defined here are only meant to control chamber
pressure and temperature.


"""


import typing
import numpy as np
from ntpSystem.errors.checkerrors import _isequallength, _inlist, _islist,\
    _issortedarray, _ispositiveArray
from ntpSystem.functions.state import State


# Possible functions:
    # linear:: ax + b
    # exponential:: a*exp(bx)+c
SHAPE_FUNCTIONS = ["linear", "exponential"]
NUM_COEFF = {"linear": 2, "exponential": 3}
# link between the functions and its intergal function
INTERGRAL_FUNCS = {"linear": "linearIntegral",
                   "exponential": "exponentialIntegral"}

class shapeFunctions:
    """pre-defined shape functions"""
    
    def __init__(self):
        pass
    
    
    def linear(self, x, coeff):
        """linear function"""
        return coeff[0]*x + coeff[1]
    
    
    def linearIntegral(self, x, coeff):
        """integral of linear function"""
        return (coeff[0]/2)*x**2 + coeff[1]*x
        
    def exponential(self, x, coeff):
        """linear function"""
        return coeff[0]*np.exp(coeff[1]*x) + coeff[2]
    
    
    def exponentialIntegral(self, x, coeff):
        """integral of linear function"""
        return (coeff[0]/coeff[1])*np.exp(coeff[1]*x) + coeff[2]*x     
    
    
class generalControlRule(shapeFunctions):
    """Define a general control rule

    The general control rule can be applied to control
    chamber temperature and pressure only. Different time intervals can have
    different control rules. The transient always begins at 0.0 seconds.


    Parameters
    ----------
    funcTypes : list fo strings
        types of control rules, e.g., ['linear', 'linear', 'exponential']
        for each time interval.
    coeffs : list of lists
        coefficients of the corresponding function type for each time interval
    tends : list of floats
        time at which each interval ends (time bounds) without the 0.0

    Raises
    ------
    ValueError
        If ``coeffs`` or ``tends`` are not the same length as ``funcTyps``.
    TypeError
        If ``funcTyps``, ``coeffs`` or ``tends`` are not lists.


    Examples
    --------
    >>> generalControlRule(['linear', 'linear'], [[1.0, 2.0], [1.5, 2.0]],
                           [30.0, 60.0])

    """    
    
    
    def __init__(self, funcTypes, coeffs, tends):
        
        self.funcTypes = funcTypes
        self.coeffs = coeffs
        self.tends = tends  # time interval bounds
        self._checkInputs()
        
        self.tends = np.array(self.tends)
    
    def _checkInputs(self):
        """check validity of variables"""
        
        # variables must be of a list type
        _islist(self.funcTypes, 'Function types')
        _islist(self.coeffs, 'Function coefficients')
        _islist(self.tends, 'Intervals end times (tends)')
        
        # lists must be all the same length
        nIntervals = len(self.funcTypes)  # number of time intervals
        _isequallength(self.coeffs, nIntervals, "Function coefficients")
        _isequallength(self.tends, nIntervals, "Intervals end times (tends)")
        
        
        for idx, fType in enumerate(self.funcTypes):
            _inlist(fType, "Function type (idx={})".format(idx),
                    SHAPE_FUNCTIONS)
            exp_ncoeff = NUM_COEFF[fType]  # expected number of coefficients
            _isequallength(self.coeffs[idx], exp_ncoeff,
                           "Function coefficients")
          
        tendsArray = np.array(self.tends)
        _issortedarray(tendsArray, "Intervals end times (tends)")
        _ispositiveArray(tendsArray, "Intervals end times (tends)")
        
    def evaluate(self, t):
        """evaluate of the function at a specific time point"""
        
        # Default value of the function is zero
        y = 0.0
        
        if t <= np.max(self.tends):                
            idx = np.where(t <= self.tends)[0][0]
            y = getattr(self, self.funcTypes[idx])(t, self.coeffs[idx])
                    
        return y
    
    
    def evaluateIntegral(self, Tstart, Tend):
        """evaluate the value of the integral using knowing limits"""
        
        y = 0.0
        
        if Tstart <= np.max(self.tends):
            
            idx0 = np.where(Tstart <= self.tends)[0][0]
            if Tend > self.tends[-1]:
                idx1 = len(Tend)
            else:
                idx1 = np.where(Tend <= self.tends)[0][0]
            
            
            intgType = INTERGRAL_FUNCS[self.funcTypes[idx0]]
            # evaluate the integral for Tstart
            y0 = getattr(self, intgType)(Tstart, self.coeffs[idx0])
 
            # start and end time are in the same interval
            if idx0 == idx1:
                y1 = getattr(self, intgType)(Tend, self.coeffs[idx0])
                # return the value of the integral
                return y1 - y0    
                
            # start and end time are NOT in the same interval
            # Loop over all the relevant intervals
            for idx in range(idx0, idx1+1):
                if idx > idx0:
                    intgType = INTERGRAL_FUNCS[self.funcTypes[idx]]
                    y0 = getattr(self, intgType)(self.tends[idx-1],
                                                 self.coeffs[idx])
                if idx < idx1:
                    intgType = INTERGRAL_FUNCS[self.funcTypes[idx]]
                    y1 = getattr(self, intgType)(self.tends[idx],
                                                 self.coeffs[idx])
                else:
                    intgType = INTERGRAL_FUNCS[self.funcTypes[idx]]
                    y1 = getattr(self, intgType)(Tend,
                                                 self.coeffs[idx])                    
                y += (y1 - y0)
                    
        return y 



# -----------------------------------------------------------------------------
#                 OLD CONTROL RULES (NOT USED ANYMORE)
#                   Written by Vignesh Manickam
# -----------------------------------------------------------------------------


class ControlRule:
    """
    Initializes ControlRule class to handle rotation of the
    control drums and mass flow transients.
    """

    def __init__(self, default: float = 0.0):
        self.default = default

    def __add__(self, other):
        if isinstance(other, int) and other == 0:
            return self
        if not isinstance(other, ControlRule):
            raise ValueError('Addition not defined for ControlRule and type: {}'.format(type(other)))
        rules = self.rules if isinstance(self, CompositeControlRule) else [self]
        if isinstance(other, CompositeControlRule):
            rules.extend(other.rules)
        else:
            rules.append(other)
        return CompositeControlRule(rules=rules)

    def __radd__(self, other):
        return self + other

    def rule_applies(self, t: float, state: State):
        raise NotImplementedError

    def drum_speed(self, t: float, state: State):
        raise NotImplementedError

    def demand_pressure(self, t: float, state: State):
        raise NotImplementedError

    def demand_temperature(self, t: float, state: State):
        raise NotImplementedError


class LinearControlRule(ControlRule):
    """
    LinearControlRule class allows for ramping the speed of
    the control drum position based on y=mx+b. Here the user
    can define the start time, ``t_min`` and finishing time,
    ``t_max`` of the control drum rotation. The rule_applies
    function verifies that these time bounds are within the time
    of the overall simulation.
    """

    def __init__(self, coeff: float, const: float, t_min: float = None,
                 t_max: float = None, default: float = 0.0):
        super().__init__(default=default)
        self.coeff = coeff
        self.const = const
        self.t_min = t_min
        self.t_max = t_max

    def __repr__(self):
        t_boundary_string = ('None' if self.t_min is None else '{:.1f}'.format(self.t_min) + ', ') + \
                            ('None' if self.t_max is None else '{:.1f}'.format(self.t_max))
        return 'LinearControlRule({:.1f}, {:.1f}, {})'.format(self.coeff, self.const, t_boundary_string)

    def rule_applies(self, t: float, state: State):
        return ((self.t_min is not None and t >= self.t_min) or self.t_min is None) and \
               ((self.t_max is not None and t <= self.t_max) or self.t_max is None)

    def drum_speed(self, t: float, state: State):
        if self.rule_applies(t, state):
            return self.coeff * t + self.const
        return self.default

    def demand_pressure(self, t: float, state: State):
        if self.rule_applies(t, state):
            return self.coeff * t + self.const
        return self.default

    def demand_temperature(self, t: float, state: State):
        if self.rule_applies(t, state):
            return self.coeff * t + self.const
        return self.default


class CompositeControlRule(ControlRule):
    """
    CompositeControlRule allows for the linear combination
    of multiple LinearControlRules. That is the user can define
    multiple control drum rotations in one single simulation. In order
    to this the user, just needs to ``+`` the LinearControlRules together
    when defining drum speed.
    """

    def __init__(self, rules: typing.Tuple[ControlRule]):
        self.rules = rules

    def rule_applies(self, t: float, state: State):
        return any(rule.rule_applies(t, state) for rule in self.rules)

    def drum_speed(self, t: float, state: State):
        return sum(rule.drum_speed(t, state) for rule in self.rules)

    def demand_pressure(self, t: float, state: State):
        return sum(rule.demand_pressure(t, state) for rule in self.rules)

    def demand_temperature(self, t: float, state: State):
        return sum(rule.demand_pressure(t, state) for rule in self.rules)
