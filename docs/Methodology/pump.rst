.. _meth-pump:


Pump
---------
The below methodology defines equations for the pump efficiency generation, primarily focusing on the
loss factors that prohibits an efficiency, :math:`\eta`, of 1. The equations described were obtained from
:ref:`Balje, Nichols, McPherson<pump-references>`. The main terms needed to describe a turbopump are known as the
similarity parameters, which are represented by :math:`N_s`, specific speed, :math:`D_s`, specific diameter, :math:`S`, suction
specific speed, and :math:`R_{e}^*`, the pump Reynolds number. It should be noted that the reference listed utilizes
imperial units for calculations.


.. math::

	N_s = \frac{N \sqrt{Q}}{H^{\frac{3}{4}}}

The parameters for the specific speed are :math:`N`, the desired rotative speed, :math:`Q`, the volumetric flow rate,
and :math:`H`, the head rise. Additionally, the general form of the specific diameter is represented with the
impeller diameter, :math:`D`, as well.

.. math::

	D_s = \frac{D H^{\frac{1}{4}}}{\sqrt{Q}}


Furthermore, the pump reynolds number utilizes the fluid density, :math:`\rho`, absolute viscosity, :math:`\mu`, and :math:`U`
which is a term described as proportional to the product of the impeller diameter and the desired rotative speed.

.. math::

	R_{e}^* = \frac{D U \rho}{\mu}


Finally, the pump specific speed is described by :ref:`Balje, Nichols, McPherson<pump-references>` 
using :math:`H_{SV}` which is utilized for cavitation, but in this scenario can be ignored and treated as :math:`H`.

.. math::

	S = \frac{N \sqrt{Q}}{H^{\frac{1}{4}}}
	
	
Additional terms such as the blade number, :math:`Z` are also utilized in further calculations using :math:`d` which is
the eye diameter as well as :math:`\beta_2` and :math:`\beta_1` which are the blade outlet and inlet angles, respectively.

.. math::

	Z = k(\frac{D_2 + d}{D_2 - d})sin(\frac{\beta_1 + \beta_2}{2})
	
	
	
.. math::

	k = 10sin(\frac{\beta_1 + \beta_2}{2})


The overall efficiency is represented by the ratio of outlet and inlet head coefficients:
:math:`q_{th}`, theoretical head coefficient,
:math:`q_{il}`, impellor loss,
:math:`q_{ad}`, annular diffuser,
:math:`q_{sc}`, scroll,
:math:`q_{sd}`, straight diffuser,
:math:`q_m`, mixing,
:math:`q_{df}`, disc friction,
:math:`q_{rf}`, wear-ring friction, and
:math:`q_{rl}`, wear-ring leakage.


.. math::

	\eta = \frac{q_{out}}{q_{in}} = \frac{q_{th} - q_{il} - q_{ad} - q_{sc} - q_{sd} - q_m}{q_{th} + q_{df} + q_{rf} + q_{rl}}



.. math::

	q_{th} = \frac{.95(1 - K_1 \phi cot\beta_2)}{1 + \frac{1 + .6sin\beta_2}{.5Z(1-.2833\delta)}}
	
	
	
Here, the denominator is denoted by :math:`m` which is known as the slip factor. Additionally, :math:`C_{m1}`
and :math:`C_{m2}` are the meridional velocities at the impeller inlet and outlet, respectively.


.. math::

	K_1 = \frac{C_{m2}}{C_{m1}}
	
	
.. math::

	\phi = \frac{C_{m1}}{U_2}
	
	
.. math::

	\delta = \frac{d}{D_2}
	
	
Furthermore, the head coefficients are further defined using :math:`\tau` which is the ratio of the hub diameter,
:math:`d_h`, to the eye diameter, as well as :math:`\frac{t}{D}` which is the ratio of the hub diameter to the eye diameter.

.. math::

	q_{il} = \frac{.0264 \phi^{1.75} [1 - \delta (\frac{1 + \tau^2}{2})^{\frac{1}{2}}] \
	(\frac{K_{1}^{2} + K_1 +1}{3})^{1.5}}{sin^{4} \beta_2 R_{e}^{* \frac{1}{4}} (1 - \tau^{2})^{1.25}} \
	(Z[\frac{ 2(\pi J - Z \frac{t}{D})sin \beta }{Z \delta^2}
	
	+ \frac{ \pi(1-\tau^2) }{ 2(\frac{K_{1}^2 + K_1 + 1}{3})^{\frac{1}{2}} (\pi J - Z \frac{t}{D}) }])^{1.25} \
	+ \frac{N_{id}}{2} [\phi^2 (1 - \frac{K_{1}^2}{sin^2 \beta_2}) + \delta^2]
	



.. math::

	J = \frac{1 - \delta (\frac{1 + \tau^2}{2})^\frac{1}{2}}{K_1 - 1}[(\frac{K_{1}^2 + K_1 + 1}{3})^\frac{1}{2} - 1] \
	+ \delta(\frac{1 + \tau^2}{2})^\frac{1}{2}
	

The equations for the head coefficients continue with that occuring due to the annular diffuser utilizing :math:`\lambda`, the ratio of
the outer diameter to the inner diameter of the annular diffuser.


.. math::

	q_{ad} = \frac{.0854 \phi^{1.75} (\lambda - 1) }{R_{e}^{* \frac{1}{4}} \delta^{2.5}} \
	[1 + \frac{1}{m^2} (\frac{1}{\phi K_1} -

	\frac{1}{tan\beta_2})^2]^{1.375}(\frac{1}{\lambda^2} + \frac{1}{\lambda} + 1)^{.375} \
	K_{1}^{3}[1 - \frac{q_{rl}}{q_{th}}]^{1.75}


The remaining numerator head loss coefficients are :math:`q_{sc}`, :math:`q_{sd}`, and :math:`q_m` which are the head loss coefficients
for the scroll, straight diffuser, and mixing, respectively.


.. math::

	q_{sc} = \frac{.155 \phi^{1.75} (1 - \frac{q_{rl}}{q_{th}})^{1.75} K_{1}^{2.375}}{\delta^{1.25} \
	{R_{e}^{* \frac{1}{4}}(1 - \tau^2) \lambda^{1.4375}}}[1 + \frac{1}{m^2}(\frac{1}{\phi K_1} - cot\beta_2)^2]^{1.188} 
	
	
	
.. math::
	
	q_{sd} = \frac{\phi^2 (1 - \eta_{sd})}{2}(1 - \frac{q_{rl}}{q_{th}})^2(\frac{K_{1}^2 [1 + \frac{1}{m^2}(\frac{1}{\phi K_1} \
	- cot\beta_2)^2]}{(\lambda + 2[\frac{m (1 - \tau^2) \delta^2}{4 K_1 (\frac{1}{\phi K_1} - cot\beta_2)} + \
	(\frac{m \lambda (1 - \tau^2) \delta^2}{4 K_1 (\frac{1}{\phi K_1} - cot\beta_2)})^{\frac{1}{2}}])^2} - 1)
	


.. math::

	q_m = .125(1 - \frac{1}{m})^2
	
	
	
The remaining pump head loss coefficients are representative of the shaft power.


.. math::

	q_{rl} = \frac{2 K_1 q_{th}}{\phi (1 - \tau^2)(\frac{k_2}{k_1} f + 1.5)^\frac{1}{2}}(2[q_{th} - q_{il} - \frac{1}{2 m^2} (1 + \
	\frac{K_{1}^2 \phi^2}{sin^2 \beta_2} - \frac{2 K_1}{tan\beta_2} ) - \frac{1 - \delta^2}{8} ])^\frac{1}{2}
	
	
	
.. math::

	f = \frac{.316}{R_{e}^{* \frac{1}{4}}}
	
	k_2 = \frac{b}{d} = \frac{b}{2R_1}
	
	k_1 = \frac{2 \Lambda}{d}
	
	
Here, :math:`f` represents the friction factor and :math:`k` is the clearance width ratio with parameters :math:`b`, width [ft],
:math:`R`, pump radius, :math:`\Lambda`, and wear-ring clearance.  



.. math::

	q_{df} = \frac{2.2(1 - \delta^2)}{\phi \delta^2 (1 - \tau^2)}
	
	
	q_{rf} = \frac{8k_2}{ln(\frac{1}{1 + k_1})R_{e}^{*} \phi \delta (1 - \tau^2)}
	


In order to fully utilize these parameters for system development there must be optimization of the combined terms.
By combining the below terms and substituting into the above equations develops an expression for efficiency that can
be used for the development of the pump system.

.. math::

	\Psi = \delta[\phi (1 - \tau^2)]^\frac{1}{2}
	
	q_{out} = \frac{11.75}{N_{s}^2 D_{s}^2}
	
	\Psi = \frac{4.85}{N_{s}^{\frac{1}{2}} D_{s}^{\frac{3}{2}}}
	
	N_s = 228.5 \frac{\Psi}{q_{out}^{\frac{3}{4}}}
	
	D_s = .473 \frac{q_{out}^{\frac{1}{4}}}{\Psi}
	
	\delta = \frac{\psi}{[\phi(1 - \tau^2)]^{\frac{1}{2}}}
	

.. _pump-references:

References
----------

O. E. Balje, K. E. Nichols, and D. G. McPherson, "Study of Turbine and Turbopump Design Parameters," final report, vol. 4,
"Low Specific Speed Turbopump Study, S / T D No. 17.35," Department
of the Navy, Office of Naval Research, Contract No. NONR 2292
