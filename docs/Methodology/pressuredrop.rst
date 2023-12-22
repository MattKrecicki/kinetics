.. _meth-pressuredrop:

Pressure loss calculations
==========================

Description
-----------

The package allows to predict the differential pressure drop
componenets. These functions allow to calculate the pressure lossses due
to static head, acceleration, friction, and form losses. These pressure
losses are calculated according to the inlet and outlet or average
properties.

Computational methodology
-------------------------

The following relation presents the **pressure loss between points
`1` and `2`** at steady-state:

:math:`\begin{equation}
\begin{split}
\Delta P =\Delta P_{acc}+\Delta P_{gravity}+\Delta P_{fric}+\Delta P_{form}
\end{split}
\end{equation}`

**Acceleration pressure losses** 
:math:`\begin{equation}
\begin{split}
&\Delta P_{acc}=G^2 \left(\frac{1}{\rho_2}-\frac{1}{\rho_1} \right) 
\end{split}
\end{equation}`

**Static head/gravity pressure losses** 
:math:`\begin{equation}
\begin{split}
&\Delta P_{gravity}=\bar{\rho}gL
\end{split}
\end{equation}`

**Friction pressure losses** 
:math:`\begin{equation}
\begin{split}
\Delta P = \phi_{heated}\times 2f \frac{L}{D_e} \frac{G^2}{\bar{\rho}}
\end{split}
\end{equation}`

where, - :math:`f` is the friction factor - :math:`G` is the mass
velocity in kg/m\ :math:`^2`/s - :math:`\phi_{heated}` is the heated
channel multiplier - :math:`L` is the length of the axial layer -
:math:`D_e` is the hydraulic diameter of the channel - :math:`\rho` is
the density of the fluid

