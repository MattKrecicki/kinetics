.. _meth-turbine:


**Turbine Efficiency**

~~~~~~~~~~~~~~~~~~

The approach numerically obtains the specific diameter of the turbine by 
guessing a specific diameter, :math:`D_{s_{guess}}`, and comparing the error in the calculated 
turbine efficiency and the desired turbine efficiency. The methodology was
obtained from :ref:`Balje, Part A. 1962<turbopump-references>`. First, the turbine inlet efficiency is found
using hydrogen property tables with the input values of the inlet pressure and enthalpy. An additional term,
:math:`\gamma_{\rho}`, the ratio of outlet and inlet densities, is used as an argument
to Newton solve for the specific diameter. To begin solving for the parameters for the turbine efficiency a
few parameters must first be calculated, including :math:`\frac{h}{D}`:


.. math::

	\frac{h}{D} = \frac{1}{D_s} ({\frac{1}{D_s [460 + (\frac{n_s D_{s}^2}{10})^2]}})^\frac{1}{3}


Next, the nozzle angle, :math:`\alpha_2` [radians] is needed:


.. math::

	\alpha_2 = \arcsin({\frac{\frac{\gamma_3}{\gamma_2}}{\pi \frac{h}{D} (1 - \frac{h}{D}) D_{s}^2 \psi_N \sqrt{2g}}})



where, 
:math:`\frac{\gamma_2}{\gamma_3}` is representative of the ratio of densities present across the rotor and 
:math:`\psi_N` is the velocity coefficient of the nozzle with a given value of 0.96, whereas the velocity coefficient
of the rotor, :math:`\psi_R`, is calculated.



.. math::

	\psi_R = [1 - \frac{\cot{\alpha_2}}{30} + \frac{n_s D_s \sqrt{1 - 2 \frac{h}{D} + 2(\frac{h}{D})^2}}{\sin{\alpha_2} \psi_N 4620}][1 - \frac{0.24}{\frac{h}{D}}]
	
	

Finally, the efficiency of the turbine, :math:`\eta_t` can be found:



.. math::

	\eta_t = \frac{n_s D_s \sqrt{1 - 2(\frac{h}{D}) + 2(\frac{h}{D})^2}}{77}[1 + \psi_R][\psi_N \cos{\alpha_2}
	- \frac{n_s D_s \sqrt{1 - 2(\frac{h}{D}) + 2(\frac{h}{D})^2}}{154}]

	- \frac{n_{s}^3 D_{s}^5 16 \beta^{\star} (1 - 2 \frac{h}{D})^5}{154^3}



where :math:`\beta^{\star}` denotes the friction coefficient of a wheel-disk which has a value of :math:`3 \times 10^{-3}` 
for a large Reynolds number.



.. _turbine-references:

References
----------

O.E. Balje, "A study on design criteria and matching of turbomachines Part B: Compressor and Pump Performance and Matching of Turbocomponets", Journal of Engineering for Power, (1962).

