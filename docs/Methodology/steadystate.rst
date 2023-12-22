.. _meth-ss:


Steady State
------------

The steady state solution is obtained through numerical iteration to balance the engine’s turbomachinery to satisfy the cryogenic tank and chamber boundary conditions.
The rocket’s geometry, including number of fuel and moderating elements, along with the relevant dimensions (e.g., engine component lengths and diameters) must be defined. 
The desired :math:`T_c` and P_c are used to calculate the overall initial system mass flow rate using Eq. (1), which presents the isentropic flow rate equation
(:ref:`Sutton and Biblarz, 2001<meth-references>`).





.. math::

   Q_{sys} = \dot{m}_{sys}\left(h_c-h_{tank} \right) \quad \textrm{Eq.(1)}
   
   
   
 d
 d
 d
 d
 d
 d
 d
 d
 d
 d
 d
 d
 
 d
 d
 d
 d
 

   
 
.. _meth-references:

References
----------

Sutton, G.P., Biblarz, O., 2001. Isentropic Flow Through Nozzles, in: Rocket Propulsion Elements. John Wiley & Sons Inc., Canada, pp. 45–102.

	