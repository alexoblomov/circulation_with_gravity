"""
Purpose : simulate ode model

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from parameters import *


def solve_Vsa(Vsa, Vsa0, Csa, Csa_l, Rs_u, Rs_l, Csa_u, Pext_l, Psa_u, Psv_u, 
              Psv_l, rho, g, H, Q):
    """ODE to solve for systemic arterial volume. Arguments are named below
    For more detailed information about units see the parameters.py file

    Args:
        Vsa (float): systemic arterial volume
        Vsa0 (float): reserve arterial volume
        Csa (float): systemic arterial compliance
        Csa_l (float): systemic upper arterial compliance
        Rs_u (float): systemic upper venous resistance
        Rs_l (float): systemic lower venous resistance
        Csa_u (float): systemic upper arterial compliance
        Pext_l (float): lower pressure external to vascular system
        Psa_u (float): pressure in upper systemic arteries
        Psv_u (float): pressure in upper systemic veins
        Psv_l (float): pressure in lower systemic veins
        rho (float): density of blood
        g (float): value of gravity
        H (float): total height
        Q (float): cardiac output

    Returns:
        float: dVsa/dt[t]
    """
    dydt = Q - Vsa/Csa*(1/Rs_u + 1/Rs_l) - (
            ((-Vsa0 + Csa_l * Pext_l - Csa_l*rho*g*H)/Csa) - Psv_u) / Rs_u - (
            ((-Vsa0 + Csa_l * Pext_l - Csa_u*rho*g*H)/Csa
            ) - Psv_l) / Rs_l # dVsa/dt
    
    return dydt

def solve_Vsv(Vsa, Vsa0, Csa, Csa_l, Rs_u, Rs_l, Csa_u, Pext_l, Psa_u, Psv_u, 
              Psv_l, rho, g, H, Q):
    """ODE to solve for systemic venous volume. Arguments are named below
    For more detailed information about units see the parameters.py file
    dVsv/dt = -dVsa/dt

    Args:
        Vsa (float): systemic arterial volume
        Vsa0 (float): reserve arterial volume
        Csa (float): systemic arterial compliance
        Csa_l (float): systemic upper arterial compliance
        Rs_u (float): systemic upper venous resistance
        Rs_l (float): systemic lower venous resistance
        Csa_u (float): systemic upper arterial compliance
        Pext_l (float): lower pressure external to vascular system
        Psa_u (float): pressure in upper systemic arteries
        Psv_u (float): pressure in upper systemic veins
        Psv_l (float): pressure in lower systemic veins
        rho (float): density of blood
        g (float): value of gravity
        H (float): total height
        Q (float): cardiac output

    Returns:
        float: dVsv/dt[t]
    """
    dydt = Q - Vsa/Csa*(1/Rs_u + 1/Rs_l) - (
            ((-Vsa0 + Csa_l * Pext_l - Csa_l*rho*g*H)/Csa) - Psv_u) / Rs_u - (
            ((-Vsa0 + Csa_l * Pext_l - Csa_u*rho*g*H)/Csa
            ) - Psv_l) / Rs_l # dVsv/dt
    return dydt

# time steps
# forward euler time step
h = 0.001
n_timesteps = 100
num_seconds = 60
T = np.linspace(0, num_seconds, num=n_timesteps)
# n_t = len(T)

# TODO : change to vary gravity as fn of time
# breakpoint()
g = g_earth*np.ones(n_timesteps)

# may want to make 2D in the future as we
# vary gravity and/or other parameters

Q = np.zeros(n_timesteps)

F = np.zeros(n_timesteps)
P_ra = np.zeros(n_timesteps)
Psv_u = np.zeros(n_timesteps)
Psv_l = np.zeros(n_timesteps)
Psa_u = np.zeros(n_timesteps)
Psa_l = np.zeros(n_timesteps)
Pext_l = np.zeros(n_timesteps)

# TODO consider moving to parameters.py
P_thorax = (- 4 * 1333) * np.ones(n_timesteps)

# TODO vary normally around the nominal values of the relative heights
Hu = Hu_patient
Hl = Hl_patient
H = Hu + Hl


Q = C_RVD * F*(P_ra - P_thorax)

# assume venous vol is roughly 70% of total BV.
# TODO: make this a normal random variable centered around 70%
pct_venous = 0.7
pct_arterial = 1 - pct_venous
init_Vsa = Vtotal * pct_arterial
init_Vsv = Vtotal * pct_venous
Vsa = np.zeros(n_timesteps+1)
Vsa[0] = init_Vsa
Vsv = np.zeros(n_timesteps+1)
Vsv[0] = init_Vsv

Vsa_s = np.zeros(n_timesteps+1)
Vsv_s = np.zeros(n_timesteps+1)

# reseve volumes
Vsa0 = np.zeros(n_timesteps)
Vsv0 = np.zeros(n_timesteps)

# for solving for vsa/vsv with odeint
t_interval = np.linspace(0, h, num=1)

for t in range(n_timesteps):
    # Case I:
    dP_RA = P_ra[t] - P_thorax[t]
    if P_thorax[t] <= - dP_RA:
        # eq 20 => six cases for Pra:
        # eq 20 gives Pra
        if Pext_l[t] < rho * g[t] * Hu:
            if P_ra[t] <= Pext_l[t]:
                P_ra[t] = (Vsv[t] - Vsv0[t] -
                           Csv_l*(rho*g[t]*Hl) +
                           Cra * P_thorax[t]) / Cra
            elif P_ra[t] >= Pext_l[t]:
                P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l *
                           (rho * g[t] * Hl - Pext_l[t]) + Cra * P_thorax[t]
                           ) / (Cra + Csv_l)
            elif rho * g[t] * Hu <= P_ra:
                P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l * (rho * g * Hl - Pext_l[t])
                           + Csv_u * rho * g * Hu + Cra *
                           P_thorax(t)) / (Cra + Csv_u + Csv_l)
        else:
            if P_ra[t] <= rho * g[t] * Hu:
                P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l * (rho * g[t] * Hl
                                                       - Pext_l[t])
                           + Cra * P_thorax[t] - Csv_l*Pext_l[t]) / Cra
            elif P_ra[t] >= rho * g[t] * Hu:
                P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l * rho * g[t] * Hl +
                           Csv_u * rho * g[t] * Hu + Cra * P_thorax[t]
                           ) / (Cra + Csv_u)
            elif Pext_l[t] <= P_ra[t]:
                P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l * (rho * g[t] * Hl -
                           Pext_l[t]) + Csv_u * rho * g[t] * Hu + Cra * P_thorax
                           ) / (Cra + Csv_u + Csv_l)

        Psv_u[t] = max(0, P_ra[t] - rho * g[t] * Hu)  # eq 15
        Psv_l[t] = max(P_ra[t], Pext_l[t]) + rho*g[t]*Hl  # eq 18

        Psa_u[t] = (Vsa[t] - Vsa0[t] + Csa_l*Pext_l[t] - Csa_l * rho * g[t]*H)/(
                    Csa)  # eq 12
        Psa_l[t] = (Vsa[t] - Vsa0[t] + Csa_l*Pext_l[t] - Csa_u * rho * g[t]*H)/(
                    Csa)  # eq 13

        # h is euler iteration timestep
        Vsa[t+1] = Vsa[t] + h * solve_Vsa(Vsa0[t], Vsa[t], Csa, Csa_l, Rs_u,
                                          Rs_l, Csa_u, Pext_l[t], Psa_u[t],
                                          Psv_u[t], Psv_l[t], rho, g[t], H,
                                          Q[t])

        Vsv[t+1] = Vsv[t] + h * solve_Vsv(Vsa0[t], Vsa[t], Csa, Csa_l, Rs_u,
                                          Rs_l, Csa_u, Pext_l[t], Psa_u[t],
                                          Psv_u[t], Psv_l[t], rho, g[t], H,
                                          Q[t])

        print("vsa+vsv ", Vsa[t+1]+Vsv[t+1])
    

    # # Case II
    # elif P_thorax[t] > - dP_RA and P_thorax[t] < rho * g[t] * Hu - dP_RA:z

    #     continue
    # # Case III
    # elif P_thorax[t] >= rho * g[t] * Hu - dP_RA:
    #     continue


