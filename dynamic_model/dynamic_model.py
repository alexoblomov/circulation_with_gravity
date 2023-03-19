"""
Purpose : simulate ode model

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp
from parameters import *

from volume_odes import *

# time steps
# forward euler time step
h = 0.0001
n_timesteps = 100
num_seconds = 60
T = np.linspace(0, num_seconds, num=n_timesteps)
# n_t = len(T)

# TODO : change to vary gravity as fn of time
# breakpoint()
g_range = g_earth*np.ones(n_timesteps)

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
    g = g_range[t]
    P_thorax = P_thorax[t]
    # Case I:
    dP_RA = P_ra[t] - P_thorax
    if P_thorax <= - dP_RA:
        # eq 20 => six cases for Pra:
        # eq 20 gives Pra
        if Pext_l[t] < rho * g * Hu:
            if P_ra[t] <= Pext_l[t]:
                P_ra[t] = (Vsv[t] - Vsv0[t] -
                           Csv_l*(rho*g*Hl) +
                           Cra * P_thorax) / Cra
            elif P_ra[t] >= Pext_l[t]:
                P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l *
                           (rho * g * Hl - Pext_l[t]) + Cra * P_thorax
                           ) / (Cra + Csv_l)
            elif rho * g * Hu <= P_ra:
                P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l * (rho * g * Hl - Pext_l[t])
                           + Csv_u * rho * g * Hu + Cra *
                           P_thorax(t)) / (Cra + Csv_u + Csv_l)
        else:
            if P_ra[t] <= rho * g * Hu:
                P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l * (rho * g * Hl
                                                       - Pext_l[t])
                           + Cra * P_thorax - Csv_l*Pext_l[t]) / Cra
            elif P_ra[t] >= rho * g * Hu:
                P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l * rho * g * Hl +
                           Csv_u * rho * g * Hu + Cra * P_thorax
                           ) / (Cra + Csv_u)
            elif Pext_l[t] <= P_ra[t]:
                P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l * (rho * g * Hl -
                           Pext_l[t]) + Csv_u * rho * g * Hu + Cra * P_thorax
                           ) / (Cra + Csv_u + Csv_l)

        Psv_u[t] = max(0, P_ra[t] - rho * g * Hu)  # eq 15
        Psv_l[t] = max(P_ra[t], Pext_l[t]) + rho*g*Hl  # eq 18

        Psa_u[t] = (Vsa[t] - Vsa0[t] + Csa_l*Pext_l[t] - Csa_l * rho * g*H)/(
                    Csa)  # eq 12
        Psa_l[t] = (Vsa[t] - Vsa0[t] + Csa_l*Pext_l[t] - Csa_u * rho * g*H)/(
                    Csa)  # eq 13

        # h is euler iteration timestep
        Vsa[t+1] = Vsa[t] + h * solve_Vsa(Vsa[t], Vsa0[t],Csa, Csa_l, Rs_u,
                                          Rs_l, Csa_u, Pext_l[t], Psa_u[t],
                                          Psv_u[t], Psv_l[t], rho, g, H,
                                          Q[t])

        Vsv[t+1] = Vsv[t] + h * solve_Vsv(Vsa[t], Vsa0[t], Csa, Csa_l, Rs_u,
                                          Rs_l, Csa_u, Pext_l[t], Psa_u[t],
                                          Psv_u[t], Psv_l[t], rho, g, H,
                                          Q[t])
        

        print("vsa+vsv ", Vsa[t+1]+Vsv[t+1])
        timvector = np.linspace(0,1,num=2)
        t_interval = np.linspace(0, h, num=2)
        print("t_interval ", t_interval, " timvector ", timvector)
        # Vsa_s[t+1] = solve_ivp(solve_ivp_Vsa, (0,1), [Vsa[t]],
        #                                                     args=(Vsa0[t],
        #                                                         Csa, Csa_l,
        #                                                         Rs_u, Rs_l,
        #                                                         Csa_u, Pext_l[t],
        #                                                         Psa_u[t], Psv_u[t], 
        #                                                         Psv_l[t], rho, g,
        #                                                         H, Q[t]))
        # print("Vsa odeint ",Vsa_s[t+1] )
        breakpoint()

    # # Case II
    # elif P_thorax > - dP_RA and P_thorax < rho * g * Hu - dP_RA:z

    #     continue
    # # Case III
    # elif P_thorax >= rho * g * Hu - dP_RA:
    #     continue


