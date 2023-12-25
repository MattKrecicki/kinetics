"""Solves the transient equations for NTPs

Updated on Mon Oct 19 2022 @author: Matt Krecicki

"""

import functools
import copy
import numpy as np
from scipy.integrate import odeint, quad
from CoolProp.CoolProp import PropsSI

from ntptransient import reactor, reactivity, turbopump, pipe, nozzle, depletion
from ntptransient.steadysolver import initTransient
from ntptransient.turbopump import pumpSideTransient, turbineSideTransient, pumpSideSteadyState
from ntptransient.utilities import load_static
from ntptransient.control import ControlRule
from ntptransient.solution import Solution
from ntptransient.state import State
from ntptransient.postprocess import export
from ntptransient.control import shapeFunctions


GRAVITY = 9.80665  # [meters/sec^2]


def intializeStateArray(inputs, ssOutputs, precursor_density_initial,
                        temp_fuel_initial):
    
    statei = State(float(ssOutputs.power), precursor_density_initial, \
        float(temp_fuel_initial), float(inputs.I135), float(inputs.Xe135),
        float(inputs.Pm149), float(inputs.Sm149), float(inputs.U235),
        float(inputs.U238), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, float(inputs.theta0),
        0.0, float(ssOutputs.j4aH), float(ssOutputs.j5aH), float(ssOutputs.j4bH),
        float(ssOutputs.j5bH), float(ssOutputs.j6H), float(ssOutputs.j7H),
        float(ssOutputs.j8H), float(ssOutputs.mixH), float(ssOutputs.j10H),
        float(ssOutputs.j11H), float(ssOutputs.j12H), float(inputs.Pci),
        float(inputs.Tci))
    
    return statei
         

def _pumpCircuit(h2props, h2propsPT, mdot, inputs, state):
    
    hinPump = h2propsPT.evaluate("h", temperature=inputs.Tt, pressure=inputs.Pt)
    
    j2P, j2H, pumpPower, head, pumpSpeed, specificSuction, pumpStress, pumpEff,\
        pumpVolFlow =\
            pumpSideTransient(h2props, hinPump, inputs.Pt, state.demand_pressure,
                              state.demand_temperature, inputs.headPickle,
                              inputs.Hdesign, mdot, inputs.noPumps,
                              inputs.pumpDesignEff, inputs.Dp, inputs.nsp,
                              inputs.pumpMaterial, inputs.volFlowDesign,
                              inputs.omegaDesign, pv=1e+4, rotorShapeFactor=0.3)
    
    turbinePower = pumpPower / inputs.turbEff
        
    return j2P, j2H, hinPump, pumpPower, pumpSpeed, head, specificSuction,\
        pumpStress, pumpEff, pumpVolFlow, turbinePower
        

def _nozzleReflectorCircuit(h2props, mdot, j3aP, j3aH, inputs,
                            state):
    
    mdotNozz = (mdot * inputs.chi) / inputs.nozzleNoChans
    
    mdotRef = (mdot * inputs.chi) / (inputs.noDrums * inputs.refNoChans)
    
    j4aP = pipe.solvePipe(h2props,
                          length=inputs.nozzleL,
                          radius=inputs.nozzleR, 
                          mass_flow=mdotNozz,
                          p_in=j3aP,
                          h_in=j3aH,
                          h_out=state.j4a_enthalpy,
                          u=inputs.uNoz)
    
    j5aP = pipe.solvePipe(h2props,
                          length=inputs.refL,
                          radius=inputs.refR, \
                          mass_flow=mdotRef,
                          p_in=j4aP,
                          h_in=state.j4a_enthalpy, 
                          h_out=state.j5a_enthalpy)
    
    return j4aP, j5aP


def _moderatorCircuit(h2props, mdot, j3bP, j3bH, j5aP, inputs, state):
    
    mdotMod = (mdot * (1 - inputs.chi)) / inputs.noModeratorElems
    
    j4bP = pipe.solvePipe(h2props, length=inputs.supplyL, mass_flow=mdotMod,
                          p_in=j3bP, h_in=j3bH, h_out=state.j4b_enthalpy,
                          flowA=inputs.supplyFlowA, hper=inputs.supplyHper,
                          u=inputs.uMEs)

    j5bP = pipe.solvePipe(h2props, length=inputs.returnL, mass_flow=mdotMod,
                          p_in=j4bP, h_in=state.j4b_enthalpy, 
                          h_out=state.j5b_enthalpy, flowA=inputs.returnFlowA,
                          hper=inputs.returnHper, u=inputs.uMEr)

    if j5aP >= j5bP:
        j6P = j5bP
    else:
        j6P = j5bP

    j7P, j8P = j6P, j6P
    
    return j4bP, j5bP, j6P, j7P, j8P


def _turbineCircuit(h2props, j8P, htank, shaftSpeed, pumpPower, pumpEff, mdot,
                    inputs, state):
    
    j9P, j9H, mdotTurbine, mdotBypass, volFlowRate, j9T, turbineStress, \
        tbcvK, tbcvTheta = turbineSideTransient(h2props, htank, \
                state.j8_enthalpy, j8P, mdot, pumpPower, inputs.noPumps,
                pumpEff, inputs.turbEff, inputs.turbPR, inputs.nst,
                inputs.Dt, inputs.turbineMaterial, shaftSpeed, inputs.tbcvA,
                inputs.tbcvB, rotorShapeFactor=0.3)

    tbcvDensity = \
        h2props.evaluate('r', pressure=(0.5*(j9P+j8P)),
                         temperature=state.j8_enthalpy)
    
    j10P = j9P
    
    return j9P, j9H, j9T, j10P, mdotTurbine, mdotBypass, tbcvDensity, tbcvK,\
        tbcvTheta, turbineStress
        
        
def _fuelCircuit(h2props, mdot, j10P, inputs, state):
    
    mdotFE = mdot / inputs.noFuelElems
    
    j11P = pipe.solvePipe(h2props, length=inputs.pipeL, radius=inputs.pipeR,\
        mass_flow=mdot, p_in=j10P, h_in=state.j10_enthalpy,
        h_out=state.j11_enthalpy)

    j12P = pipe.solvePipe(h2props, length=inputs.fuelL, mass_flow=mdotFE,\
        p_in=j11P, h_in=state.j11_enthalpy, h_out=state.j12_enthalpy,
        flowA=inputs.fuelFlowA, hper=inputs.fuelHper, u=inputs.uFE)
    
    return j11P, j12P


def _obtainStationTemperatures(h2props, j2P, j2H, j3aP, j3aH, j3bH, j4aP, j5aP,
                               j3bP, j4bP, j5bP, j6P, j7P, j8P, j10P, j11P,
                               j12P, state):
    
    j2T = h2props.evaluate("T", pressure=j2P, temperature=j2H)
    j3aT = h2props.evaluate("T", pressure=j3aP, temperature=j3aH)
    j4aT = h2props.evaluate("T", pressure=j4aP, temperature=state.j4a_enthalpy)
    j5aT = h2props.evaluate("T", pressure=j5aP, temperature=state.j5a_enthalpy)
    j3bT = h2props.evaluate("T", pressure=j3bP, temperature=j3bH)
    j4bT = h2props.evaluate("T", pressure=j4bP, temperature=state.j4b_enthalpy)
    j5bT = h2props.evaluate("T", pressure=j5bP, temperature=state.j5b_enthalpy)
    j6T = h2props.evaluate("T", pressure=j6P, temperature=state.j6_enthalpy)
    j7T = h2props.evaluate("T", pressure=j7P, temperature=state.j7_enthalpy)
    j8T = h2props.evaluate("T", pressure=j8P, temperature=state.j8_enthalpy)
    j10T = h2props.evaluate("T", pressure=j10P, temperature=state.j10_enthalpy)
    j11T = h2props.evaluate("T", pressure=j11P, temperature=state.j11_enthalpy)
    j12T = h2props.evaluate("T", pressure=j12P, temperature=state.j12_enthalpy)
    
    return j2T, j3aT, j4aT, j5aT, j3bT, j4bT, j5bT, j6T, j7T, j8T, j10T,\
        j11T, j12T


def _calcHeatTransferCoefficient(h2props, mdot, j11P, j12P, inputs, state):
    
    Deh = 4 * inputs.fuelFlowA / inputs.fuelHper
    
    j11_cp = h2props.evaluate("cp", pressure=j11P, temperature=state.j11_enthalpy)
    j11_my = h2props.evaluate("my", pressure=j11P, temperature=state.j11_enthalpy)
    j11_tc = h2props.evaluate("tc", pressure=j11P, temperature=state.j11_enthalpy)
    j11_pr = (j11_cp * j11_my) / j11_tc

    j12_cp = h2props.evaluate("cp", pressure=j12P, temperature=state.j12_enthalpy)
    j12_my = h2props.evaluate("my", pressure=j12P, temperature=state.j12_enthalpy)
    j12_tc = h2props.evaluate("tc", pressure=j12P, temperature=state.j12_enthalpy)
    j12_pr = (j12_cp * j12_my) / j12_tc

    avCp = 0.5 * (j11_cp + j12_cp)
    avMy = 0.5 * (j11_my + j12_my)
    avPr = 0.5 * (j11_pr + j12_pr)
    avTc = 0.5 * (j11_tc + j12_tc)
    avRe = ((Deh * mdot) / (avMy * inputs.fuelFlowA))

    ht = pipe.heatTransferCoefficient(avRe, avPr, avTc, Deh)
    
    return ht


def _calcGasDensities(h2props, j3bP, j4bP, j5bP, j11P, j12P, j3bH, state):
    
    j3b_density = h2props.evaluate("r", pressure=j3bP, temperature=j3bH)
    j4b_density = h2props.evaluate("r", pressure=j4bP, temperature=state.j4b_enthalpy)
    j5b_density = h2props.evaluate("r", pressure=j5bP, temperature=state.j5b_enthalpy)
    
    supplyDensity = 0.5*(j3b_density + j4b_density)
    returnDensity = 0.5*(j4b_density + j5b_density)
    
    j11Rho = h2props.evaluate("r", pressure=j11P, temperature=state.j11_enthalpy)
    j12Rho = h2props.evaluate("r", pressure=j12P, temperature=state.j12_enthalpy)
    fuelDensity = 0.5*(j11Rho + j12Rho)

    mes_avCp = 0.5 * (h2props.evaluate("cp", pressure=j3bP, temperature=j3bH) +
                      h2props.evaluate("cp", pressure=j4bP, temperature=state.j4b_enthalpy))

    mer_avCp = 0.5 * (h2props.evaluate("cp", pressure=j4bP, temperature=state.j4b_enthalpy) +
                      h2props.evaluate("cp", pressure=j5bP, temperature=state.j5b_enthalpy))
    
    avCpFuel = 0.5 * (h2props.evaluate("cp", pressure=j11P, temperature=state.j11_enthalpy) +
                      h2props.evaluate("cp", pressure=j12P, temperature=state.j12_enthalpy))
    
    return supplyDensity, returnDensity, fuelDensity, mes_avCp, mer_avCp, avCpFuel


def _calcEnthalpyDerivatives(h2props, inputs, state, mdot, pumpmdot, noPumps,
                             turbinemdot, j3aH, j3bH, j3aP, j3bP, j4aP, j4bP,
                             j5aP, j5bP, j6P, j9H, j9P, j10P, j11P, j12P,
                             prev_j3aP, prev_j3bP, prev_j4aP, prev_j4bP,
                             prev_j5aP, prev_j5bP, prev_j6P, prev_j9P,
                             prev_j10P, prev_j11P, prev_j12P, t, prev_t):
    
    dj4aHdt = pipe.enthalpy_out_deriv(h2props,
                                      length=inputs.nozzleL,
                                      radius=inputs.nozzleR,
                                      volMod = inputs.nozzleNoChans,
                                      mass_flow=mdot*inputs.chi,
                                      heat=state.power * inputs.pfNozzle,
                                      h_in=j3aH,
                                      h_out=state.j4a_enthalpy,
                                      p_in=j3aP,
                                      p_out=j4aP,
                                      prev_p_in=prev_j3aP,
                                      prev_p_out=prev_j4aP,
                                      t=t,
                                      prev_t=prev_t)

    dj5aHdt = pipe.enthalpy_out_deriv(h2props,
                                      length=inputs.refL,
                                      radius=inputs.refR,
                                      volMod=inputs.refNoChans,
                                      mass_flow=(mdot * inputs.chi),
                                      heat=(state.power * inputs.pfReflector),
                                      h_in=state.j4a_enthalpy,
                                      h_out=state.j5a_enthalpy,
                                      p_in=j4aP,
                                      p_out=j5aP,
                                      prev_p_in=prev_j4aP,
                                      prev_p_out=prev_j5aP,
                                      t=t,
                                      prev_t=prev_t)

    dj4bHdt = pipe.enthalpy_out_deriv(h2props,
                                      length=inputs.supplyL, 
                                      radius=(inputs.supplyHper / (2 * np.pi)),
                                      volMod=inputs.noModeratorElems,
                                      mass_flow=(mdot * (1 - inputs.chi)),
                                      heat=(state.power * inputs.pfSupply),
                                      h_in=j3bH,
                                      h_out=state.j4b_enthalpy,
                                      p_in=j3bP,
                                      p_out=j4bP,
                                      prev_p_in=prev_j3bP,
                                      prev_p_out=prev_j4bP,
                                      t=t,
                                      prev_t=prev_t)

    dj5bHdt = pipe.enthalpy_out_deriv(h2props,
                                      length=inputs.returnL,
                                      radius=(inputs.returnHper / (2 * np.pi)),
                                      volMod=inputs.noModeratorElems,
                                      mass_flow=(mdot * (1 - inputs.chi)),
                                      heat=(state.power * inputs.pfReturn),
                                      h_in=state.j4b_enthalpy,
                                      h_out=state.j5b_enthalpy,
                                      p_in=j4bP,
                                      p_out=j5bP,
                                      prev_p_in=prev_j4bP,
                                      prev_p_out=prev_j5bP,
                                      t=t,
                                      prev_t=prev_t)

    dj6Hdt = pipe.enthalpy_out_join_tee_deriv(h2props,
                                              length=inputs.pipeL,
                                              radius=inputs.pipeR,
                                              mdotA=mdot*inputs.chi,
                                              mdotB=mdot*(1-inputs.chi),
                                              hinA=state.j5a_enthalpy,
                                              hinB=state.j5b_enthalpy,
                                              hout=state.j6_enthalpy,
                                              pinA=j5aP,
                                              pinB=j5bP,
                                              pout=j6P,
                                              prev_pinA=prev_j5aP,
                                              prev_pinB=prev_j5bP,
                                              prev_pout=prev_j6P,
                                              t=t,
                                              prev_t=prev_t)

    dj7Hdt = dj6Hdt

    dj8Hdt = dj6Hdt
    
    #A = turbine outlet, junction 9
    bypass_mdot = (pumpmdot - turbinemdot)
    
    #B = tubrine bypass line, junction 7
    
    dmixHdt = \
        pipe.enthalpy_out_join_tee_deriv(h2props,
                                         length=inputs.pipeL,
                                         radius=inputs.pipeR,
                                         mdotA=turbinemdot,
                                         mdotB=bypass_mdot,
                                         hinA=j9H,
                                         hinB=state.j8_enthalpy,
                                         hout=state.mix_enthalpy,
                                         pinA=j9P,
                                         pinB=j9P,
                                         pout=j9P,
                                         prev_pinA=prev_j9P,
                                         prev_pinB=prev_j9P,
                                         prev_pout=prev_j9P,
                                         t=t,
                                         prev_t=prev_t)
    
    #hmix = 3479674.3347098688
    #h10 = 6959348.6694197375
    #P9 = 4809879.670162864 / 4809854.563910066
    #P10 = 4809879.670162864 / 4809854.563910066
    
    dj10Hdt = pipe.enthalpy_out_join_tee_deriv(h2props,
                                               length=inputs.pipeL,
                                               radius=inputs.pipeR,
                                               mdotA=pumpmdot,
                                               mdotB=pumpmdot,
                                               hinA=state.mix_enthalpy,
                                               hinB=state.mix_enthalpy,
                                               hout=state.j10_enthalpy,
                                               pinA=j9P,
                                               pinB=j9P,
                                               pout=j10P,
                                               prev_pinA=prev_j9P,
                                               prev_pinB=prev_j9P,
                                               prev_pout=prev_j10P,
                                               t=t,
                                               prev_t=prev_t)
    
    
    dj11Hdt = pipe.enthalpy_out_deriv(h2props,
                                      length=inputs.pipeL,
                                      radius=inputs.pipeR,
                                      mass_flow=mdot,
                                      heat=(state.power * inputs.pfPipe),
                                      h_in=state.j10_enthalpy,
                                      h_out=state.j11_enthalpy,
                                      p_in=j10P,
                                      p_out=j11P,
                                      prev_p_in=prev_j10P,
                                      prev_p_out=prev_j11P,
                                      t=t,
                                      prev_t=prev_t)

    dj12Hdt = pipe.enthalpy_out_deriv(h2props, 
                                      length=inputs.fuelL,
                                      mass_flow=mdot,
                                      heat=(state.power * inputs.pfFuel),
                                      h_in=state.j11_enthalpy,
                                      h_out=state.j12_enthalpy,
                                      p_in=j11P,
                                      p_out=j12P,
                                      prev_p_in=prev_j11P,
                                      prev_p_out=prev_j12P,
                                      t=t,
                                      prev_t=prev_t,
                                      hper=inputs.fuelHper,
                                      volMod=inputs.noFuelElems)


    return dj4aHdt, dj5aHdt, dj4bHdt, dj5bHdt, dj6Hdt, dj7Hdt, dj8Hdt,\
        dmixHdt, dj10Hdt, dj11Hdt, dj12Hdt


def _calcReactivityFeedback(inputs, state, t, prev_t, j4bT, j5bT, j11T, j11P,
                            j12T, j4bP, j5bP, j12P, dj12Hdt, dj4bHdt, dj5bHdt,
                            htFuel, avCpFuel, mes_avCp, mer_avCp, fuelDens,
                            supplyDens, returnDens):
    
    # Compute derivatives for integrated quantities----------------------------
    dndt = reactor.total_neutron_deriv(beta=inputs.beta.sum(),
                                       gen_time=inputs.genTime,
                                       power=state.power,
                                       precursor_constants=inputs.Lambda,
                                       precursor_density=state.precursor_densities,
                                       rho_fuel_temp=state.rho_fuel_temp,
                                       rho_fuel_gas_density=state.rho_fuel_gas_density,
                                       rho_mes_gas_density=state.rho_mes_gas_density,
                                       rho_mer_gas_density=state.rho_mer_gas_density,
                                       rho_Xe_135=state.rho_Xe_135,
                                       rho_Sm_149=state.rho_Sm_149,
                                       rho_con_drum=state.rho_con_drum)
    
    # Caluclate precusor conc. time derv--------------------------------------- 
    dcdt = reactor.delay_neutron_deriv(beta_vector=inputs.beta,
                                       gen_time=inputs.genTime,
                                       power=state.power,
                                       precursor_constants=inputs.Lambda,
                                       precursor_density=state.precursor_densities)
    
    powerFE = (state.power * inputs.pfFuel)
    
    dTfdt = reactor.fuel_temp_deriv(power=powerFE,
                                    mass_fuel=inputs.massFuel,
                                    heat_cap_fuel=inputs.cpFuel,
                                    heat_coeff=htFuel,
                                    temp_fuel=state.t_fuel,
                                    temp_h2_fuel=0.5*(j11T + j12T))
    
    
    # Caluclate fission product reactivity deriv-------------------------------
    dI135dt, dXe135dt, dPm149dt, dSm149dt, dU235dt, dU238dt = \
        depletion.dndt([state.I_135, state.Xe_135, state.Pm_149, state.Sm_149,
                        state.U_235, state.U_238], t=t, power=state.power,
                       volume=inputs.fuelVol, Tc=j12T, Pc=j12P)
    
    # Calculate fuel temperature reactivity derviative ------------------------
    
    dRhoTfdt = reactor.temp_fuel_reactivity_deriv(power=powerFE,
                                                  mass_fuel=inputs.massFuel,
                                                  heat_cap_fuel=inputs.cpFuel,
                                                  heat_coeff=htFuel,
                                                  temp_fuel=state.t_fuel,
                                                  temp_h2_fuel=0.5*(j11T + j12T))
    
    
    # Calculate fuel gas density reactivity derviative ------------------------
    
    alphaFEdens =\
        reactivity.alpha_fuel_gas_density(fuel_gas_density=fuelDens,
                                          gas_temperature=0.5*(j11T + j12T),
                                          gas_pressure=0.5*(j11P + j12P)/1e+6)
    
    dRhoFEdt = (1 / avCpFuel) * dj12Hdt 
    
    dRhoFuelDensdt = alphaFEdens * dRhoFEdt
    
    
    # Calculate supply gas density reactivity derviative ----------------------
    
    alphaSupplyDens = reactivity.alpha_mes_gas_density(mes_gas_density=supplyDens,
                                                       gas_temperature=j4bT,
                                                       gas_pressure=j4bP/1e+6)
    dSupplyDensitydt = (1 / mes_avCp) * dj4bHdt
    
    dRhoSupplyDensdt = alphaSupplyDens * dSupplyDensitydt
    
    # Calculate return gas density reactivity derviative ----------------------

    alphaReturnDensity = reactivity.alpha_mer_gas_density(mer_gas_density=returnDens,
                                                          gas_temperature=j5bT,
                                                          gas_pressure=j5bP/1e+6)
    
    dReturnDensitydt = (1 / mer_avCp) * dj5bHdt
                     
    dRhoReturnDensdt = alphaReturnDensity * dReturnDensitydt

    dRhoXedt = reactivity.alpha_xenon(Xe_135=state.Xe_135) * dXe135dt

    dRhoSmdt = reactivity.alpha_samarium(Sm_149=state.Sm_149) * dSm149dt
    
    return dndt, dcdt, dTfdt, dI135dt, dXe135dt, dPm149dt, dSm149dt, dU235dt,\
        dU238dt, dRhoTfdt, dRhoFuelDensdt, dRhoSupplyDensdt, dRhoReturnDensdt,\
        dRhoXedt, dRhoSmdt


def state_deriv_array(state_array: np.ndarray, t, inputs, powerTarget,
                      state_history):    
    """Function to compute transient equations for NTPs"""

    # Access current state
    state = State.from_array(state_array)
    
    # Access previous state
    prev_t, prev_state = state_history['prev']

    # Load non-integrated quantites from previous iteration
    prev_non_integrated = state_history['selected_data']
    
    # Check to see if non-integrated quantities have been recorded, else use initial condition
    prev_chamber_pressure = prev_non_integrated[-1][15] if prev_non_integrated else inputs.Pci

    #split out h2props for easy of use
    h2props, h2propsPT = inputs.H2_PH, inputs.H2_PT

    # Compute non-integrated quantities----------------------------------------
    mdot = \
        nozzle.massflow(h2propsPT, inputs.At, state.demand_pressure,
                 state.demand_temperature)
    
    # Solve Pump Circuit-------------------------------------------------------
        
    j2P, j2H, hinPump, pumpPower, pumpSpeed, head, specificSuction, pumpStress,\
        pumpEff, pumpVolFlow, turbinePower = \
            _pumpCircuit(h2props, h2propsPT, mdot, inputs, state)
    
    j3aH, j3aP = j2H, j2P
    j3bH, j3bP = j2H, j2P
    
    
    #Solve Nozzle and Reflector Circuit---------------------------------------
    
    j4aP, j5aP = \
        _nozzleReflectorCircuit(h2props, mdot, j3aP, j3aH, inputs, state)
    
    
    #Solve Moderator Circuit---------------------------------------------------

    j4bP, j5bP, j6P, j7P, j8P = \
        _moderatorCircuit(h2props, mdot, j3bP, j3bH, j5aP, inputs, state)

    
    #Solve Turbine Circuit-----------------------------------------------------

    j9P, j9H, j9T, j10P, mdotTurbine, mdotBypass, tbcvDensity, tbcvK,\
        tbcvTheta, turbineStress =\
            _turbineCircuit(h2props, j8P, hinPump, pumpSpeed, pumpPower,
                            pumpEff, mdot, inputs, state)
    turbinePower = pumpPower/(pumpEff * inputs.turbEff)
    

    #Solve Pipe and Fuel Element Circuit---------------------------------------
    
    j11P, j12P = _fuelCircuit(h2props, mdot, j10P, inputs, state)
    

    # Calculate Station Temperatures-------------------------------------------
    
    j2T, j3aT, j4aT, j5aT, j3bT, j4bT, j5bT, j6T, j7T, j8T, j10T,\
        j11T, j12T = _obtainStationTemperatures(h2props, j2P, j2H, j3aP, j3aH,\
            j3bH, j4aP, j5aP, j3bP, j4bP, j5bP, j6P, j7P, j8P, j10P, j11P,
            j12P, state)
            
    # Calculate fuel heat transfer coefficient---------------------------------
    
    htfuel = _calcHeatTransferCoefficient(h2props, mdot, j11P, j12P, inputs,
                                         state)
    
    
    # Track densities and average specific heat values for gas density reactivity insertion
    
    supplyDensity, returnDensity, fuelDensity, mes_avCp, mer_avCp, avCpFuel =\
        _calcGasDensities(h2props, j3bP, j4bP, j5bP, j11P, j12P, j3bH, state)
    
    
    # Track non-integrated quantities------------------------------------------
    
    state_history['selected_data'].append((t, mdot, mdotTurbine, mdotBypass,\
        j3aP, j3aT, j4aP, j4aT, j5aP, j5aT, j3bP, j3bT, j4bP, j4bT, j5bP, j5bT,
        j6P, j6T, j7P, j7T, j8P, j8T, j9P, j9T, j10P, j10T, j11P, j11T, j12P,
        j12T, pumpSpeed, tbcvTheta, pumpPower, turbinePower, specificSuction,
        pumpStress, pumpEff, turbineStress))


    # Access 'previous' non-integrated quantities------------------------------
    
    prev_t, prev_mass_flow, prev_turbine_mass_flow, prev_bypass_mass_flow, \
        prev_j3aP, prev_j3aT, prev_j4aP, prev_j4aT, prev_j5aP, prev_j5aT, \
        prev_j3bP, prev_j3bT, prev_j4bP, prev_j4bT, prev_j5bP, prev_j5bT, \
        prev_j6P, prev_j6T, prev_j7P, prev_j7T, prev_j8P, prev_j8T, prev_j9P, \
        prev_j9T, prev_j10P, prev_j10T, prev_j11P, prev_j11T, prev_j12P, \
        prev_j12T, prev_pump_speed, prev_tbcv_position, prev_pump_p, \
        prev_turbine_p, prev_specificSuction, prev_pumpStress, prev_pumpEff,\
        prev_turbineStress = state_history['selected_data'][-1]


    #Calculate change in entalpy time derivatives---------- --------------------
    
    dj4aHdt, dj5aHdt, dj4bHdt, dj5bHdt, dj6Hdt, dj7Hdt, dj8Hdt, dmixHdt, \
        dj10Hdt, dj11Hdt, dj12Hdt =\
            _calcEnthalpyDerivatives(h2props, inputs, state, mdot, mdot/inputs.noPumps,\
                inputs.noPumps, mdotTurbine, j3aH, j3bH, j3aP, j3bP, j4aP,
                j4bP, j5aP, j5bP, j6P, j9H, j9P, j10P, j11P, j12P, prev_j3aP,
                prev_j3bP, prev_j4aP, prev_j4bP, prev_j5aP, prev_j5bP,
                prev_j6P, prev_j9P, prev_j10P, prev_j11P, prev_j12P, t, prev_t)
    
    
    # Calculate Reactivity Feedbacks-------------------------------------------
    
    dndt, dcdt, dTfdt, dI135dt, dXe135dt, dPm149dt, dSm149dt, dU235dt, dU238dt,\
        dRhoTfdt, dRhoFuelDensdt, dRhoSupplyDensdt, dRhoReturnDensdt, dRhoXedt,\
        dRhoSmdt = \
            _calcReactivityFeedback(inputs, state, t, prev_t, j4bT, j5bT, j11T,\
                j11P, j12T, j4bP, j5bP, j12P, dj12Hdt, dj4bHdt, dj5bHdt, 
                htfuel, avCpFuel, mes_avCp, mer_avCp, fuelDensity, supplyDensity,
                returnDensity)
                
    
    # Calculate reactivity control via MIT-SNL laws----------------------------
    
    drhodtFeedback = dRhoTfdt + dRhoFuelDensdt + dRhoSupplyDensdt + \
        dRhoReturnDensdt + dRhoXedt + dRhoSmdt

    control_signal, drho_con_drum_dt_automatic, ddrum_angle_dt_automatic = \
        reactor.auto_con_drum_reactivity_deriv(drum_angle=state.drum_angle,
                                               beta=inputs.beta.sum(),
                                               beta_vector=inputs.beta,
                                               gen_time=inputs.genTime,
                                               power=state.power,
                                               desired_power=powerTarget,
                                               precursor_constants=inputs.Lambda,
                                               precursor_density=state.precursor_densities,
                                               rho_fuel_temp=state.rho_fuel_temp,
                                               rho_fuel_gas_density=state.rho_fuel_gas_density,
                                               rho_mes_gas_density=state.rho_mes_gas_density,
                                               rho_mer_gas_density=state.rho_mer_gas_density,
                                               rho_Xe_135=state.rho_Xe_135,
                                               rho_Sm_149=state.rho_Sm_149,
                                               rho_con_drum=state.rho_con_drum,
                                               time_ramp=inputs.dtRamp,
                                               t=t,
                                               rho_feedback_deriv=drhodtFeedback)
    
    dpdt = inputs.pressureControlRule.evaluate(t)

    dtdt = inputs.temperatureControlRule.evaluate(t)


    # Update time derivatives--------------------------------------------------
    
    state_deriv = State(dndt, dcdt, dTfdt, dI135dt, dXe135dt, dPm149dt, \
        dSm149dt, dU235dt, dU238dt, dRhoTfdt, dRhoFuelDensdt, dRhoSupplyDensdt,
        dRhoReturnDensdt, dRhoXedt, dRhoSmdt, ddrum_angle_dt_automatic,
        drho_con_drum_dt_automatic, dj4aHdt, dj5aHdt, dj4bHdt, dj5bHdt, dj6Hdt,
        dj7Hdt, dj8Hdt, dmixHdt, dj10Hdt, dj11Hdt, dj12Hdt, dpdt, dtdt)


    # Set state history--------------------------------------------------------
    
    state_history['prev'] = (t, state)

    return state_deriv.to_array()


def _obtainStartupStatePoints(inputs):
    
    #setup solution array
    startupArray = {}
        
    #get inital conditions
    Tci, Pci = inputs.Tci, inputs.Pci
    
    #expand out control rules for evaluation
    TcControl, PcControl = inputs.temperatureControlRule, \
        inputs.pressureControlRule
        
    #get all unique time points for evaluation
    Ts_Tc = np.round(TcControl.tends, 3)
    Ts_Pc = np.round(PcControl.tends, 3)
    
    Ts = np.unique(np.append(Ts_Tc, Ts_Pc))
    
    #Perform integration for chamber conditions time points
    Tstart = 0.0
    
    Tct, Pct, time = [Tci], [Pci], [0.0]
    
    for Ti in Ts:
        
        Tend = Ti
        
        #chamber temperature integral
        To = TcControl.evaluateIntegral(Tstart, Tend) + Tci
        Tct.append(Tct)
        
        #chamber pressure integral
        Po = PcControl.evaluateIntegral(Tstart, Tend) + Pci
        Pct.append(Po)
        
        time.append(Ti)
        
        #update parameters
        Tci = To
        Pci = Po
        T0 = Ti

    
    #append results to dictionary
    startupArray["Tc"] = Tct
    startupArray["Pc"] = Pct
    startupArray["t"] = Ts
    
    return startupArray


def solve(inputs, outputfile="output.hdf5", writeoutput=True, verbose=True, 
          tolTc=1e-10, tolPc=1e-10, method="BFGS", maxiter=100):
    """Solving differential equations to calculate time dependent preformance
    of NTP reactor"""
    
    #correct chmaber conditions so ss solver get correct inital conditions-----
    
    #expand out hydrogen properties--------------------------------------------
    h2props = inputs.H2_PH
    
    #get target chamber condition array for target power controller 
    #startupArray = _obtainStartupStatePoints(inputs)
    
    
    # compute target power-----------------------------------------------------
    inputs.Pc, inputs.Tc = inputs.Pcf, inputs.Tcf
    powerTarget = \
        copy.deepcopy(initTransient(inputs, tolTc, tolPc, method, maxiter,
                                    False).power)
    
    #solve inital conditions---------------------------------------------------
    inputs.Pc, inputs.Tc = inputs.Pci, inputs.Tci
    ssOutputs = initTransient(inputs, tolTc, tolPc, method, maxiter, False)      
    inputs.H = ssOutputs.pumpHead
    

    # Initialize precursor concentrations--------------------------------------
    precursor_density_initial = inputs.beta / (inputs.Lambda * inputs.genTime)\
        * float(ssOutputs.power)

    # Output in temperatures---------------------------------------------------
    j11Ti = h2props.evaluate("T", pressure=ssOutputs.j11P,
                             temperature=ssOutputs.j11H)
    j12Ti = h2props.evaluate("T", pressure=ssOutputs.j12P,
                             temperature=ssOutputs.j12H)

    # Initialize fuel temperature----------------------------------------------
    j11_myi = h2props.evaluate("my", pressure=ssOutputs.j11P,
                               temperature=ssOutputs.j11H) 
    j11_cpi = h2props.evaluate("cp", pressure=ssOutputs.j11P,
                               temperature=ssOutputs.j11H) 
    j11_tci = h2props.evaluate("tc", pressure=ssOutputs.j11P,
                               temperature=ssOutputs.j11H) 
    j11_pri = (j11_cpi * j11_myi) / j11_tci

    j12_myi = h2props.evaluate("my", pressure=ssOutputs.j12P,
                               temperature=ssOutputs.j12H) 
    j12_cpi = h2props.evaluate("cp", pressure=ssOutputs.j12P,
                               temperature=ssOutputs.j12H) 
    j12_tci = h2props.evaluate("tc", pressure=ssOutputs.j12P,
                               temperature=ssOutputs.j12H) 
    j12_pri = (j12_cpi * j12_myi) / j12_tci

    De = 4 * inputs.fuelFlowA / inputs.fuelHper
    avMyi = 0.5 * (j11_myi + j12_myi)
    avPri = 0.5 * (j11_pri + j12_pri)
    avTci = 0.5 * (j11_tci + j12_tci)
    avRei = ((De * ssOutputs.mdot) / (avMyi * inputs.fuelFlowA))
    
    heat_transfer = pipe.heatTransferCoefficient(avRei, avPri, avTci, De)

    temp_fuel_initial = (0.5 * (j11Ti + j12Ti)) + \
                        ((ssOutputs.power * inputs.pfFuel) / heat_transfer)
                        
    # Builds a state array of all initial conditions to pass into odeint
    statei = \
        intializeStateArray(inputs, ssOutputs, precursor_density_initial,
                            temp_fuel_initial)

    # Compute time intervals for odeint integrator
    t = np.linspace(0.0, inputs.tmax, inputs.noIters)
    
    # ODE Integrator configuration
    interval = (inputs.tmax) / inputs.noIters

    # Create state history to hold "previous" time step
    # Setup storage array for non-integrated quantities
    # (massflow, pressure) --> 2 quantities to track at each "selected" time step
    state_history = {'prev': (0.0, statei),
                     'selected_data': []}
    
    # Partialize the state derivative function for signature compatibility with odeint
    deriv_func = functools.partial(state_deriv_array, inputs=inputs,\
        state_history=state_history, powerTarget=powerTarget)
        
    # Compute result using odeint integrator
    res = odeint(deriv_func, statei.to_array(), t)
    
    #writes results to hdf5 file if requested
    if writeoutput:
        export(outputfile, res, np.array(state_history['selected_data']),\
           np.arange(0.0, inputs.tmax, interval))
        
    # Returns result to Solution class for each state array for each time iteration
    res = Solution(array=res, t=np.arange(0.0, inputs.tmax, interval),
                   non_integrated_data=np.array(state_history['selected_data']))
    
    return res

