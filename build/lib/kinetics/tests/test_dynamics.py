"""Unittests for the inhour module
"""

import numpy

from ntptransient import reactor, reactivity
from ntptransient.config import Configs


class TestDynamics:
    cfg = Configs.Test
    POWER_INITIAL = cfg['reactor']['physical']['power_initial']
    BETA_VECTOR = cfg['reactor']['physical']['beta_vector']
    PRECURSOR_CONSTANTS = cfg['reactor']['physical']['precursor_constants']
    BETA = cfg['reactor']['physical']['beta']
    GEN_TIME = cfg['reactor']['physical']['gen_time']
    HEAT_COEFF = cfg['reactor']['thermo']['heat_coeff']
    MASS_H2 = cfg['reactor']['thermo']['mass_h2']
    HEAT_CAP_H2 = cfg['reactor']['thermo']['heat_cap_h2']
    MASS_FLOW_INITIAL = cfg['reactor']['thermo']['mass_flow_initial']
    MASS_FUEL = cfg['reactor']['thermo']['mass_fuel']
    HEAT_CAP_FUEL = cfg['reactor']['thermo']['heat_cap_fuel']
    TEMP_IN = cfg['reactor']['thermo']['temp_in']
    DRUM_ANGLE_INITIAL = cfg['control']['drum_angle_initial']
    DRUM_CONTROL_RULE = cfg['control']['rule']
    MASS_FLOW_RULE = cfg['flow']['rule']
    T_MAX = cfg['time']['t_max']
    NUM_ITERS = cfg['time']['num_iters']

    PRECURSOR_DENSITY_INITIAL = BETA_VECTOR / (PRECURSOR_CONSTANTS * GEN_TIME) * POWER_INITIAL
    TEMP_H2_FUEL_INITIAL = TEMP_IN + (POWER_INITIAL / (2 * MASS_FLOW_INITIAL * HEAT_CAP_H2))
    TEMP_FUEL_INITIAL = TEMP_IN + (1 / (2 * MASS_FLOW_INITIAL * HEAT_CAP_H2) + (1 / HEAT_COEFF)) * POWER_INITIAL

    RHO_FUEL_TEMP_INTIAL = reactivity.reactivity(reactivity.alpha_fuel(TEMP_FUEL_INITIAL), TEMP_FUEL_INITIAL, TEMP_FUEL_INITIAL)
    RHO_H2_FUEL_TEMP_INITIAL = reactivity.reactivity(reactivity.alpha_coolant(TEMP_H2_FUEL_INITIAL), TEMP_H2_FUEL_INITIAL, TEMP_H2_FUEL_INITIAL)
    RHO_CON_DRUM_INITIAL = reactivity.reactivity(reactivity.alpha_drum(DRUM_ANGLE_INITIAL), DRUM_ANGLE_INITIAL, DRUM_ANGLE_INITIAL)

    def test_total_neutron_deriv(self):
        res = reactor.total_neutron_deriv(beta=self.BETA,
                                          gen_time=self.GEN_TIME,
                                          power=self.POWER_INITIAL,
                                          precursor_constants=self.PRECURSOR_CONSTANTS,
                                          precursor_density=self.PRECURSOR_DENSITY_INITIAL,
                                          rho_fuel_temp=self.RHO_FUEL_TEMP_INTIAL,
                                          rho_h2_fuel_temp=self.RHO_H2_FUEL_TEMP_INITIAL,
                                          rho_con_drum=self.RHO_CON_DRUM_INITIAL)
        numpy.testing.assert_almost_equal(actual=res, desired=51256.35008239746, decimal=5)
