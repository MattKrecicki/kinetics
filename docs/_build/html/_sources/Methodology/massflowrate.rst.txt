.. _meth-mdot:


Nozzle Flow and Throat Condition
--------------------------------

Methodology was obtained from :ref:`Sutton and Biblarz, 2001<mdot-references>` (See pages 47 - 59, eq. 3-24). 
The method relies on the assumption that the flow at the nozzle is choked, this is 
not always true for very low chamber temperatures and pressure and in the
future should be expanded to account for subsonic flows in the nozzle.


If sonic velocity is reached at any point within a steady flow system, it
is impossible for a pressure disturbance to travel past the location of sonic or
supersonic flow. It is not possible to increase the throat velocity
or the flow rate in the nozzle by further lowering the exit pressure or even
evacuating the exhaust section. This important condition is often described as
choking the flow. It is always established at the throat and not the nozzle exit
plane. It is equal to the mass flow at any
section within the nozzle.





.. math::

	 \dot{m} = A_t P_1 \gamma \frac{\sqrt{[2/(\gamma+1)]^{\frac{\gamma+1}{\gamma-1}}}}{\sqrt{\gamma R T_1}}



where, 
:math:`A_t` is the throat area [:math:`m^2`], 
:math:`P_c` is the stagnation chamber pressure [Pascals],
:math:`\gamma` is the ratio of the specific heats, 
:math:`R` is the universal gas constant, and
:math:`T_c` is the chamber temperature  [Kelvin].

The terms for :math:`R` and :math:`\gamma` are determined by



.. math::

	R = \frac{R^{'}}{MW}



where, :math:`MW` represents the molecular weight of the propellant (:math:`H_2`)


.. math::

	\gamma = \frac{c_p}{c_v}
	

where, :math:`c_p` and :math:`c_v` represent the specific heat capacity at constant 
pressure and volume, respectively. The value for :math:`c_p` is found from the
hydrogen property tables, while :math:`c_v` is simply the difference
between :math:`c_p` and :math:`R`.

Meanwhile, the mach number at the throat, :math:`M_t`, must be determined to evaluate the nozzle type 
as sonic, subsonic, and supersonic nozzles have differing throat velocities, exit
velocities, and pressure ratios. The purpose of this is to ensure that the mach number is 1.0 
and therefore sonic, meaning moving at the velocity of sound.


.. math::

	M_t = \frac{v_t}{a_t}


:math:`v_t` is the fluid velocity at the throat and :math:`a_t` is the
velocity of sound in the chamber.



.. math::

	v_t = \sqrt{\frac{2 \gamma}{ \gamma + 1} R T_c}
	

.. math::
	
	
	a_t = {\sqrt{\gamma R T_t}}



:math:`T_t` represents the fluid temperature at the throat [Kelvin]



.. math::

	T_t = \frac{2 T_c}{\gamma + 1} 



 


   
 
.. _mdot-references:

References
----------

Sutton, G.P., Biblarz, O., 2001. Isentropic Flow Through Nozzles, in: Rocket Propulsion Elements. John Wiley & Sons Inc., Canada, 7th Edition. See pages 47 - 59, eq. 3-24.

	