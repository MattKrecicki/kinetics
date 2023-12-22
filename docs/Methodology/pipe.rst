.. _meth-pipe:


Pipe
====

General Pipe
------------

The models for pipining are based on a control volume approach.
The methodology was obtained from :ref:`J. K. Witter, 1993<pipe-references>` (Pages 42-44, Eq. 2.4). 

Several governing assumptions are adopted to model the time-dependent behaviour in a pipe:

- We rely on conversation of mass formalism;
- The velocity of the flow in the pipe is much smaller than Mach speeds, and thus the flow is assumed to be constant at any time step.
- This quasi-static state, incompressible assumption causes the time derivative to be zero.
- Finally, the simplifying assumption is that the time rate of change of the control volume enthalpy is approximately the same as that of the exit junction enthalpy.

The final form of the energy equation used by the transient model is represented as:

.. math::

	 \frac{dh}{dt}=\frac{1}{M}\left[\dot{m} (h_{in}-h_{out})+ Q - V\frac{dP}{dt} \right]




where, 
:math:`M` mass of the propellant (e.g., hydrogen) in the volume [:math:`kg`], 
:math:`Q` heat generation [:math:`watts`],
:math:`V` control volume [:math:`m^3`],
:math:`\dot{m}` mass flow rate [:math:`kg/sec`],

:math:`R` is the universal gas constant, and
:math:`T_c` is the chamber temperature  [Kelvin].
units are all SI.



Tee components
--------------

Tee junctions are used to connect the turbine bypass control valve lines to the main piping.
The methodology was obtained from :ref:`J. K. Witter, 1993<pipe-references>` (Pages 46-48, Eq. 2.15). 
Several assumptions are made to simplify the analysis:

- Tee junction properties are the same as the control volume to which it is attached.
- Pressure drop is calculated from the straight line path. A turning loss coefficient, that can be specified by the user, is applied to the tee junction.


.. math::

	 \frac{dh}{dt}=\frac{1}{M}\left[\dot{m_a}h_{in,a}+\dot{m_b}h_{in,b}-\dot{m}h_{out}+ Q - V\frac{dP}{dt} \right]


where, 
:math:`h_{in,a}` Enthalpy of hydrogen from join split inlet 'a', in units of Joules, 
:math:`h_{in,a}` Enthalpy of hydrogen from join split inlet 'b', in units of Joules, 
:math:`\dot{m_a}` mass Flow Rate into tee from junction A, in units of kg/s, 
:math:`\dot{m_b}` mass Flow Rate into tee from junction B, in units of kg/s, 
:math:`\dot{m}=\dot{m_a}+\dot{m_b}`   
 
.. _pipe-references:

References
----------

J. K. Witter, 1993. Modeling for the simulation and control of nuclear reactor rocket systems, MIT Ph.D. Disseration, (1993), Pages 42-44, Eq. 2.4	