"""Creates reactivity coefficients based on inputs from Serpent

Updated on Tues Aug 31 20:10:00 2021 @author: Vigneshwar Manickam
Updated on Mon Dec 06 20:10:00 2021 @author: Vigneshwar Manickam
"""

from numpy.polynomial import polynomial as p
from ntptransient.utilities import load_static


def ddens_dtemp(gas_temperature: float, gas_pressure: float):
    """Returns the partial derivative of gas density

    Parameters
    ----------
    gas_temperature: float
        Temperature of gas in K
    gas_pressure: float
        Pressure of gas in MPa

    Returns
    -------
    ddens_dtemp: float
        Returns the partial derivative of gas density
    """

    return -2.37738852E2 * (gas_pressure / (gas_temperature * gas_temperature))


def alpha_fuel(fuel_temperature: float):
    """Cached version of reactivity coefficient due to fuel temperature
    """
    if load_static.REACTIVITY_POLYNOMIAL_1 is None:
        load_static.load_reactivity_polynomial()
    d = load_static.REACTIVITY_POLYNOMIAL_1
    coeff = p.polyval(fuel_temperature, d)

    return coeff


def alpha_fuel_gas_density(fuel_gas_density: float, gas_temperature: float, gas_pressure: float):
    """Cached version of reactivity coefficient due to fuel gas density
    """
    if load_static.REACTIVITY_POLYNOMIAL_2 is None:
        load_static.load_reactivity_polynomial()
    d = load_static.REACTIVITY_POLYNOMIAL_2
    coeff = p.polyval(fuel_gas_density, d)

    return coeff * ddens_dtemp(gas_temperature=gas_temperature, gas_pressure=gas_pressure)


def alpha_mes_gas_density(mes_gas_density: float, gas_temperature: float, gas_pressure: float):
    """Cached version of reactivity coefficient due to moderator supply gas density
    """
    if load_static.REACTIVITY_POLYNOMIAL_3 is None:
        load_static.load_reactivity_polynomial()
    d = load_static.REACTIVITY_POLYNOMIAL_3
    coeff = p.polyval(mes_gas_density, d)

    return coeff * ddens_dtemp(gas_temperature=gas_temperature, gas_pressure=gas_pressure)


def alpha_mer_gas_density(mer_gas_density: float, gas_temperature: float, gas_pressure: float):
    """Cached version of reactivity coefficient due to moderator return gas density
    """
    if load_static.REACTIVITY_POLYNOMIAL_4 is None:
        load_static.load_reactivity_polynomial()
    d = load_static.REACTIVITY_POLYNOMIAL_4
    coeff = p.polyval(mer_gas_density, d)

    return coeff * ddens_dtemp(gas_temperature=gas_temperature, gas_pressure=gas_pressure)


def alpha_drum(drum_angle: float):
    """Cached version of reactivity coefficient due to drum angle position
    """
    if load_static.REACTIVITY_POLYNOMIAL_5 is None:
        load_static.load_reactivity_polynomial()
    d = load_static.REACTIVITY_POLYNOMIAL_5
    coeff = p.polyval(drum_angle, d)

    return coeff


def alpha_xenon(Xe_135: float):
    """Cached version of reactivity coefficient due to Xe-135 concentration
    """

    if load_static.REACTIVITY_POLYNOMIAL_6 is None:
        load_static.load_reactivity_polynomial()
    d = load_static.REACTIVITY_POLYNOMIAL_6
    coeff = p.polyval(Xe_135, d)

    return coeff


def alpha_samarium(Sm_149: float):
    """Cached version of reactivity coefficient due to Sm-149 concentration
    """

    if load_static.REACTIVITY_POLYNOMIAL_7 is None:
        load_static.load_reactivity_polynomial()
    d = load_static.REACTIVITY_POLYNOMIAL_7
    coeff = p.polyval(Sm_149, d)

    return coeff


def reactivity(alpha: float, x: float, x_initial: float) -> float:
    return alpha * (x - x_initial)
