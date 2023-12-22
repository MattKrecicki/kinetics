\"""state.py

Creates State Array to pass into ODEINT for the NTP Transient Equations

"""

import enum
import numpy


class StateComponent(enum.IntEnum):
    """StateComponent enumeration defines the array-index components
    of the state vector Initializes an array that holds each reactor state
    for the simulation. For example Hydrogen Gas in Fuel Element is index 7
    in the array.
    """
    Power = 0
    PrecursorDensity1 = 1
    PrecursorDensity2 = 2
    PrecursorDensity3 = 3
    PrecursorDensity4 = 4
    PrecursorDensity5 = 5
    PrecursorDensity6 = 6
    TFuel = 7
    I135 = 8
    Xe135 = 9
    Pm149 = 10
    Sm149 = 11
    U235 = 12
    U238 = 13
    RhoFuelTemp = 14
    RhoFuelGasDensity = 15
    RhoMESGasDensity = 16
    RhoMERGasDensity = 17
    RhoXe135 = 18
    RhoSm149 = 19
    DrumAngle = 20
    RhoConDrum = 21
    J4AEnthalpy = 22
    J5AEnthalpy = 23
    J4BEnthalpy = 24
    J5BEnthalpy = 25
    J6Enthalpy = 26
    J7Enthalpy = 27
    J8Enthalpy = 28
    MixEnthalpy = 29
    J10Enthalpy = 30
    J11Enthalpy = 31
    J12Enthalpy = 32
    DemandPressure = 33
    DemandTemperature = 34


class NonIntComponent(enum.IntEnum):
    """NonIntComponent enumeration defines the array-index components
    of the non-integrated quantities vector:

    t, mass_flow, turbine_mass_flow, bypass_mass_flow, j3a_pressure, j3a_temperature, j4a_pressure, j4a_temperature,
    j5a_pressure, j5a_temperature, j3b_pressure, j3b_temperature, j4b_pressure, j4b_temperature, j5b_pressure,
    j5b_temperature, j6_pressure, j6_temperature, j7_pressure, j7_temperature, j8_pressure, j8_temperature,
    j9_pressure, j9_temperature,  j10_pressure, j10_temperature, j11_pressure, j11_temperature, j12_pressure,
    j12_temperature, pump_speed, tbcv_position, pump_p, turbine_p
    """
    Time = 0
    MassFlow = 1
    TurbineMassFlow = 2
    BypassMassFlow = 3
    J3APressure = 4
    J3ATemperature = 5
    J4APressure = 6
    J4ATemperature = 7
    J5APressure = 8
    J5ATemperature = 9
    J3BPressure = 10
    J3BTemperature = 11
    J4BPressure = 12
    J4BTemperature = 13
    J5BPressure = 14
    J5BTemperature = 15
    J6Pressure = 16
    J6Temperature = 17
    J7Pressure = 18
    J7Temperature = 19
    J8Pressure = 20
    J8Temperature = 21
    J9Pressure = 22
    J9Temperature = 23
    J10Pressure = 24
    J10Temperature = 25
    J11Pressure = 26
    J11Temperature = 27
    J12Pressure = 28
    J12Temperature = 29
    PumpSpeed = 30
    TBCVPosition = 31
    PumpPower = 32
    TurbinePower = 33


class State:
    __slots__ = ('power', 'precursor_densities', 't_fuel', 'I_135', 'Xe_135', 'Pm_149', 'Sm_149', 'U_235', 'U_238',
                 'rho_fuel_temp', 'rho_fuel_gas_density', 'rho_mes_gas_density', 'rho_mer_gas_density', 'rho_Xe_135',
                 'rho_Sm_149', 'drum_angle', 'rho_con_drum', 'j4a_enthalpy', 'j5a_enthalpy', 'j4b_enthalpy',
                 'j5b_enthalpy', 'j6_enthalpy', 'j7_enthalpy', 'j8_enthalpy', 'mix_enthalpy', 'j10_enthalpy',
                 'j11_enthalpy', 'j12_enthalpy', 'demand_pressure', 'demand_temperature')

    def __init__(self, power: float, precursor_densities: numpy.ndarray, t_fuel: float, I_135: float, Xe_135: float,
                 Pm_149: float, Sm_149: float, U_235: float, U_238: float, rho_fuel_temp: float,
                 rho_fuel_gas_density: float, rho_mes_gas_density: float, rho_mer_gas_density: float, rho_Xe_135: float,
                 rho_Sm_149: float, drum_angle: float, rho_con_drum: float, j4a_enthalpy: float, j5a_enthalpy: float,
                 j4b_enthalpy: float, j5b_enthalpy: float, j6_enthalpy: float, j7_enthalpy: float, j8_enthalpy: float,
                 mix_enthalpy: float, j10_enthalpy: float, j11_enthalpy: float, j12_enthalpy: float,
                 demand_pressure: float, demand_temperature: float):
        """
        Args:
            power
            precursor_densities
            t_fuel
            I_135
            Xe_135
            Pm_149
            Sm_149
            U_235
            U_238
            rho_fuel_temp
            rho_fuel_gas_density
            rho_mes_gas_density
            rho_mer_gas_density
            rho_Xe_135
            rho_Sm_149
            drum_angle
            rho_con_drum
            j4a_enthalpy
            j5a_enthalpy
            j4b_enthalpy
            j5b_enthalpy
            j6_enthalpy
            j7_enthalpy
            j8_enthalpy
            mix_enthalpy
            j10_enthalpy
            j11_enthalpy
            j12_enthalpy
            demand_pressure
            demand_temperature
        """
        self.power = power
        self.precursor_densities = precursor_densities
        self.t_fuel = t_fuel
        self.I_135 = I_135
        self.Xe_135 = Xe_135
        self.Pm_149 = Pm_149
        self.Sm_149 = Sm_149
        self.U_235 = U_235
        self.U_238 = U_238
        self.rho_fuel_temp = rho_fuel_temp
        self.rho_fuel_gas_density = rho_fuel_gas_density
        self.rho_mes_gas_density = rho_mes_gas_density
        self.rho_mer_gas_density = rho_mer_gas_density
        self.rho_Xe_135 = rho_Xe_135
        self.rho_Sm_149 = rho_Sm_149
        self.drum_angle = drum_angle
        self.rho_con_drum = rho_con_drum
        self.j4a_enthalpy = j4a_enthalpy
        self.j5a_enthalpy = j5a_enthalpy
        self.j4b_enthalpy = j4b_enthalpy
        self.j5b_enthalpy = j5b_enthalpy
        self.j6_enthalpy = j6_enthalpy
        self.j7_enthalpy = j7_enthalpy
        self.j8_enthalpy = j8_enthalpy
        self.mix_enthalpy = mix_enthalpy
        self.j10_enthalpy = j10_enthalpy
        self.j11_enthalpy = j11_enthalpy
        self.j12_enthalpy = j12_enthalpy
        self.demand_pressure = demand_pressure
        self.demand_temperature = demand_temperature


    def to_array(self):
        return numpy.concatenate((numpy.array([self.power]),
                                  self.precursor_densities,
                                  numpy.array([self.t_fuel,
                                               self.I_135,
                                               self.Xe_135,
                                               self.Pm_149,
                                               self.Sm_149,
                                               self.U_235,
                                               self.U_238,
                                               self.rho_fuel_temp,
                                               self.rho_fuel_gas_density,
                                               self.rho_mes_gas_density,
                                               self.rho_mer_gas_density,
                                               self.rho_Xe_135,
                                               self.rho_Sm_149,
                                               self.drum_angle,
                                               self.rho_con_drum,
                                               self.j4a_enthalpy,
                                               self.j5a_enthalpy,
                                               self.j4b_enthalpy,
                                               self.j5b_enthalpy,
                                               self.j6_enthalpy,
                                               self.j7_enthalpy,
                                               self.j8_enthalpy,
                                               self.mix_enthalpy,
                                               self.j10_enthalpy,
                                               self.j11_enthalpy,
                                               self.j12_enthalpy,
                                               self.demand_pressure,
                                               self.demand_temperature])), axis=0)


    @staticmethod
    def from_array(state_array: numpy.ndarray):
        return State(power=state_array[StateComponent.Power],
                     precursor_densities=state_array[StateComponent.PrecursorDensity1:StateComponent.TFuel],
                     t_fuel=state_array[StateComponent.TFuel],
                     I_135=state_array[StateComponent.I135],
                     Xe_135=state_array[StateComponent.Xe135],
                     Pm_149=state_array[StateComponent.Pm149],
                     Sm_149=state_array[StateComponent.Sm149],
                     U_235=state_array[StateComponent.U235],
                     U_238=state_array[StateComponent.U238],
                     rho_fuel_temp=state_array[StateComponent.RhoFuelTemp],
                     rho_fuel_gas_density=state_array[StateComponent.RhoFuelGasDensity],
                     rho_mes_gas_density=state_array[StateComponent.RhoMESGasDensity],
                     rho_mer_gas_density=state_array[StateComponent.RhoMERGasDensity],
                     rho_Xe_135=state_array[StateComponent.RhoXe135],
                     rho_Sm_149=state_array[StateComponent.RhoSm149],
                     drum_angle=state_array[StateComponent.DrumAngle],
                     rho_con_drum=state_array[StateComponent.RhoConDrum],
                     j4a_enthalpy=state_array[StateComponent.J4AEnthalpy],
                     j5a_enthalpy=state_array[StateComponent.J5AEnthalpy],
                     j4b_enthalpy=state_array[StateComponent.J4BEnthalpy],
                     j5b_enthalpy=state_array[StateComponent.J5BEnthalpy],
                     j6_enthalpy=state_array[StateComponent.J6Enthalpy],
                     j7_enthalpy=state_array[StateComponent.J7Enthalpy],
                     j8_enthalpy=state_array[StateComponent.J8Enthalpy],
                     mix_enthalpy=state_array[StateComponent.MixEnthalpy],
                     j10_enthalpy=state_array[StateComponent.J10Enthalpy],
                     j11_enthalpy=state_array[StateComponent.J11Enthalpy],
                     j12_enthalpy=state_array[StateComponent.J12Enthalpy],
                     demand_pressure=state_array[StateComponent.DemandPressure],
                     demand_temperature=state_array[StateComponent.DemandTemperature])
