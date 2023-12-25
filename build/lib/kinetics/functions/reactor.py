"""Returns the transient equation for reactor kinetics and reactivity
due to drum rotation

Updated on Tues Jul 28 09:10:00 2020 @author: Vigneshwar Manickam
Updated on Mon Dec 06 20:10:00 2021 @author: Vigneshwar Manickam
"""

import numpy

from ntptransient import reactivity


#################################################
#                POINT  KINETICS               #
#################################################

def total_neutron_deriv(beta: float, gen_time: float, power: float, precursor_constants: numpy.ndarray,
                        precursor_density: numpy.ndarray, rho_fuel_temp: float, rho_fuel_gas_density: float,
                        rho_mes_gas_density: float, rho_mer_gas_density: float, rho_Xe_135: float, rho_Sm_149: float,
                        rho_con_drum: float) -> float:
    """Calculates time derivative of total neutron population (i.e. reactor power)

    Parameters
    ----------
    beta: float
        Delayed Neutron Fraction
    gen_time: float
        Effective generation time in seconds. Must be positive
    power: float
        Total Reactor power in Watts. Must be positive
    precursor_constants: numpy.ndarray
        Array of precursor constants
    precursor_density: numpy.ndarray
        Precursor Concentration
    rho_fuel_temp: float
        Reactivity due to fuel temperature
    rho_fuel_gas_density: float
        Reactivity due to hydrogen gas density in the fuel element
    rho_mes_gas_density: float
        Reactivity due to hydrogen gas density in the moderator supply channel
    rho_mer_gas_density: float
        Reactivity due to hydrogen gas density in the moderator return channel
    rho_Xe_135: float
        Reactivity due to Xe-135 concentration
    rho_Sm_149: float
        Reactivity due to Sm-149 concentration
    rho_con_drum: float
        Reactivity due to control drum position

    Returns
    -------
        float, time derivative of power
    """

    total_rho = rho_fuel_temp + rho_fuel_gas_density + rho_mes_gas_density + rho_mer_gas_density + rho_Xe_135 + \
                rho_Sm_149 + rho_con_drum

    return (((total_rho - beta) / gen_time) * power) + \
           numpy.inner(precursor_constants, precursor_density)


def delay_neutron_deriv(beta_vector: numpy.ndarray, gen_time: float, power: float, precursor_constants: numpy.ndarray,
                        precursor_density: numpy.ndarray) -> numpy.ndarray:
    """Compute time derivative of delayed neutron population

    Parameters
    ----------
    beta_vector: ndarray
        Delayed Neutron Fraction for each delayed neutron group
    gen_time: float
        Effective generation time in seconds
    power: float
        Reactor power in Watts
    precursor_constants: ndarray
        Decay Constant for each delayed neutron group in sec^-1
    precursor_density: ndarray
        Population of each delayed neutron group

    Returns
    -------
    dcdt: ndarray
        Time derivative of delayed neutron population

    Note
    ----
    * ``beta``, ``gen_time``, ``power``, ``precursor_constants``, ``precursor_densities``
    must be positive values
    """

    return beta_vector * power / gen_time - precursor_constants * precursor_density


def fuel_temp_deriv(power: float, mass_fuel: float, heat_cap_fuel: float, heat_coeff: float, temp_fuel: float,
                    temp_h2_fuel: float) -> float:
    """Compute time derivative of fuel temperature,

    Parameters
    ----------
    power: float
        Reactor power in Watts
    mass_fuel: float
        Mass of fuel in kg
    heat_cap_fuel: float
        Specific Heat capacity of fuel in J/kg/K
    heat_coeff: float
        Heat transfer coefficient of fuel and hydrogen in J/K/sec
    temp_fuel: float
        Fuel Temperature in K
    temp_h2_fuel: float
        Hydrogen Gas Temperature in FE in K

    Returns
    -------
    dT_fueldt: float
        Time derivative of fuel temperature in K/sec

    Note
    ----
    * All parameters in this function must be positive values
    """

    return (power / (mass_fuel * heat_cap_fuel)) - \
           ((heat_coeff / (mass_fuel * heat_cap_fuel)) * (temp_fuel - temp_h2_fuel))


def temp_fuel_reactivity_deriv(power: float, mass_fuel: float, heat_cap_fuel: float, heat_coeff: float,
                               temp_fuel: float, temp_h2_fuel: float) -> float:
    """Compute time derivative reactivity due to fuel temperature,

    Parameters
    ----------
    power: float
        Reactor power in Watts
    mass_fuel: float
        Mass of fuel in kg
    heat_cap_fuel: float
        Specific Heat capacity of fuel in J/kg/K
    heat_coeff: float
        Heat transfer coefficient of fuel and hydrogen in J/K/sec
    temp_fuel: float
        Fuel Temperature in K
    temp_h2_fuel: float
        Hydrogen Gas Temperature in FE in K
    Returns
    -------
    drho_fuel_temp_dt: float
        Time derivative of fuel reactivity temperature in dk/sec
    Note
    ----
    * All parameters in this function must be positive values
    """

    return reactivity.alpha_fuel(temp_fuel) * \
           fuel_temp_deriv(power=power, mass_fuel=mass_fuel, heat_cap_fuel=heat_cap_fuel, heat_coeff=heat_coeff,
                           temp_fuel=temp_fuel, temp_h2_fuel=temp_h2_fuel)


def manual_con_drum_reactivity_deriv(drum_speed: float, drum_angle: float) -> float:
    """Compute time derivative of reactivity due to drum position,

    Parameters
    ----------
    drum_angle: float
        Position of drums measured in degrees
    drum_speed: float
        Speed of drum rotation as defined in ControlRule Class
    Note
    ----
    * All ``beta`` and ``drum_angle`` must be positive values.
    * ``drum_angle`` must be between 0 to 180 degrees.
    """

    """Verifying intial conditions are defined properly"""
    if 0 > drum_angle > 180:
        raise ValueError("Drum Angle must be within 0 to 180 degrees")

    return reactivity.alpha_drum(drum_angle) * drum_speed


#################################################
#                CONTROL   KINETICS             #
#################################################

def desired_inverse_period(desired_power: float, power: float, time_interval: float) -> float:
    """Calculates the desired reactor period

    Parameters
    ----------
    desired_power: float
        User-defined reactor power in Watts. Must be positive
    power: float
        Reactor power in Watts. Must be positive
    time_interval: float
        Duration of the sample interval. Must be positive.

    Returns
    -------
        float, desired inverse reactor period
    """

    return numpy.log(desired_power / power) / time_interval


def measured_inverse_period(power: float, beta: float, gen_time: float, precursor_constants: numpy.ndarray,
                            precursor_density: numpy.ndarray, rho_fuel_temp: float, rho_fuel_gas_density: float,
                            rho_mes_gas_density: float, rho_mer_gas_density: float, rho_Xe_135: float,
                            rho_Sm_149: float, rho_con_drum: float) -> float:
    """Calculates the actual reactor period

    Parameters
    ----------
    power: float
        Total Reactor power in Watts. Must be positive
    beta: float
        Delayed Neutron Fraction. Must be positive
    gen_time: float
        Effective generation time in seconds. Must be positive
    precursor_constants: numpy.ndarray
        Array of precursor constants
    precursor_density: numpy.ndarray
        Precursor Concentration
    rho_fuel_temp: float
        Reactivity due to fuel temperature
    rho_fuel_gas_density: float
        Reactivity due to hydrogen gas density in the fuel element
    rho_mes_gas_density: float
        Reactivity due to hydrogen gas density in the moderator supply channel
    rho_mer_gas_density: float
        Reactivity due to hydrogen gas density in the moderator return channel
    rho_Xe_135: float
        Reactivity due to Xe-135 concentration
    rho_Sm_149: float
        Reactivity due to Sm-149 concentration
    rho_con_drum: float
        Reactivity due to control drum position

    Returns
    -------
        float, actual inverse reactor period
    """

    total_rho = rho_fuel_temp + rho_fuel_gas_density + rho_mes_gas_density + rho_mer_gas_density + rho_Xe_135 + \
                rho_Sm_149 + rho_con_drum

    return ((total_rho - beta) / gen_time) + (numpy.inner(precursor_constants, precursor_density) / power)


def effective_precursor_constant(precursor_constants: numpy.ndarray, precursor_density: numpy.ndarray) -> float:
    """Calculates effective precursor constant for Alternate MIT Period-Generated Control Law

    Parameters
    ----------
    precursor_constants: numpy.ndarray
        Array of precursor constants
    precursor_density: numpy.ndarray
        Precursor Concentration

    Returns
    -------
        float, effective precursor constant
    """

    return numpy.inner(precursor_constants ** 2, precursor_density) / \
        numpy.inner(precursor_constants, precursor_density)


def auto_con_drum_reactivity_deriv(drum_angle, beta, beta_vector, gen_time,
                                   power, desired_power, precursor_constants,
                                   precursor_density, rho_fuel_temp,
                                   rho_fuel_gas_density, rho_mes_gas_density,
                                   rho_mer_gas_density, rho_Xe_135,
                                   rho_Sm_149, rho_con_drum, time_ramp, t,
                                   rho_feedback_deriv):
    """Compute time derivative of reactivity due to automatic drum rotation,

     Parameters
    ----------
    drum_angle: float
        Position of drums measured in degrees
    drum_speed: float
        Speed of drum rotation
    beta: float
        Delayed Neutron Fraction
    beta_vector: ndarray
        Delayed Neutron Fraction for each delayed neutron group
    gen_time: float
        Effective generation time in seconds. Must be positive
    power: float
        Total Reactor power in Watts. Must be positive
    desired_power: float
        User-defined reactor power in Watts. Must be positive
    power_initial: float
        Initial Reactor power in Watts. Must be positive
    precursor_constants: numpy.ndarray
        Array of precursor constants
    precursor_density: numpy.ndarray
        Precursor Concentration
    rho_fuel_temp: float
        Reactivity due to fuel temperature
    rho_fuel_gas_density: float
        Reactivity due to hydrogen gas density in the fuel element
    rho_mes_gas_density: float
        Reactivity due to hydrogen gas density in the moderator supply channel
    rho_mer_gas_density: float
        Reactivity due to hydrogen gas density in the moderator return channel
    rho_Xe_135: float
        Reactivity due to Xe-135 concentration
    rho_Sm_149: float
        Reactivity due to Sm-149 concentration
    rho_con_drum: float
        Reactivity due to control drum position
    time_ramp: float
        Duration of power ramp in seconds. Must be positive.
    t: float
        Time of ODEINT iterator in seconds. Must be positive
    mass_fuel: float
        Mass of fuel in kg
    heat_cap_fuel: float
        Specific Heat capacity of fuel in J/kg/K
    heat_coeff: float
        Heat transfer coefficient of fuel and hydrogen in J/K/sec
    mass_h2: float
        Mass of hydrogen in kg
    heat_cap_h2: float
        Specific Heat capacity of hydrogen in J/kg/K
    mass_flow: float
        Total hydrogen mass flow rate in in FE in kg/sec
    temp_fuel: float
        Fuel Temperature in K
    temp_h2_fuel: float
        Hydrogen Gas Temperature in FE in K
    temp_in: float
        Inlet Coolant Temperature in K

    Note
    ----
    * All ``beta`` and ``drum_angle`` must be positive values.
    * ``drum_angle`` must be between 0 to 180 degrees.
    """

    """Verifying initial conditions are defined properly"""
    if 0 > drum_angle > 180:
        raise ValueError("Drum Angle must be within 0 to 180 degrees")

    global power_initial

    if t == 0:
        power_initial = power
    else:
        pass

    total_rho = rho_fuel_temp + rho_fuel_gas_density + rho_mes_gas_density + \
        rho_mer_gas_density + rho_Xe_135 + rho_Sm_149 + rho_con_drum

    # Effective precursor constant
    lamda_eff = effective_precursor_constant(precursor_constants=precursor_constants,
                                             precursor_density=precursor_density)

    # Measured inverse period
    if t == 0:
        wt = 0
    else:
        wt = \
            measured_inverse_period(power=power, beta=beta, gen_time=gen_time,\
                precursor_constants=precursor_constants, 
                precursor_density=precursor_density, rho_fuel_temp=rho_fuel_temp,
                rho_fuel_gas_density=rho_fuel_gas_density, \
                rho_mes_gas_density=rho_mes_gas_density, 
                rho_mer_gas_density=rho_mer_gas_density, rho_Xe_135=rho_Xe_135,
                rho_Sm_149=rho_Sm_149, rho_con_drum=rho_con_drum)

    if desired_power > power_initial:
        if power >= desired_power:
            ws = 0
            w_dot = (ws - wt) / time_ramp
        else:
            ws = \
                desired_inverse_period(desired_power=desired_power,
                                       power=power_initial,
                                       time_interval=time_ramp)
            w_dot = (ws - wt) / time_ramp

        control_signal = ((beta - total_rho) * ws) - (lamda_eff * total_rho) - \
                         numpy.inner(beta_vector, (precursor_constants - lamda_eff)) - rho_feedback_deriv + \
                         (gen_time * (w_dot + ws * (ws + lamda_eff)))

        drum_speed = control_signal / reactivity.alpha_drum(drum_angle)

        drho_con_drum_dt = reactivity.alpha_drum(drum_angle) * drum_speed
        
    else:
        if power <= desired_power:
            ws = 0
            w_dot = (ws - wt) / time_ramp
        else:
            ws = desired_inverse_period(desired_power=desired_power, power=power_initial, time_interval=time_ramp)
            w_dot = (ws - wt) / time_ramp

        control_signal = ((beta - total_rho) * ws) - (lamda_eff * total_rho) - \
                         numpy.inner(beta_vector, (precursor_constants - lamda_eff)) - rho_feedback_deriv + \
                         (gen_time * (w_dot + ws * (ws + lamda_eff)))

        drum_speed = control_signal / reactivity.alpha_drum(drum_angle)

        drho_con_drum_dt = reactivity.alpha_drum(drum_angle) * drum_speed

    return control_signal, drho_con_drum_dt, drum_speed
