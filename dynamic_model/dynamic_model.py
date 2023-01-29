"""
Purpose : simulate ode model

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from parameters import *


T = np.linspace(0,100,1)
n_t = len(T)

# TODO : change to vary gravity as fn of time
g = g_earth*np.ones(T)

# may want to make 2D in the future as we
# vary gravity and/or other parameters

Q = np.zeros(T)

# does F need to be a fn of time?
F = np.zeros(T)
P_ra = np.zeros(T)
Psv_u = np.zeros(T)
Psv_l = np.zeros(T)
Psa_u = np.zeros(T)
Psa_l = np.zeros(T)
Pext_l = np.zeros(T)

# TODO consider moving to parameters.py
P_thorax = (- 4 * 1333) * np.ones(T)

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
Vsa = np.zeros(T) ; Vsa[0] = init_Vsa
Vsv = np.zeros(T) ; Vsv[0] = init_Vsv

# reseve volumes
Vsa0 = np.zeros(T)
Vsv0 = np.zeros(T)

for t in T:
    # Case I:
    dP_RA = P_ra[t] - P_thorax[t]
    if P_thorax[j] <= - dP_RA:
        # eq 20 => six cases for Pra:
        # eq 20 gives Pra
        if Pext_l[t] < rho* g[t] * Hu:
            # /!\circular logic.
            if P_ra[t] <= Pext_l[t]:
                P_ra[t] = (Vsv[t] - Vsv0[t] -
                            Csv_l*(rho*g[t]*Hl) +
                            Cra * P_thorax[t]) / Cra
            elif P_ra[t] >= Pext_l[t]:
                P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l * 
                            (rho * g[t] * Hl - Pext_l[t]) + Cra * P_thorax[t]
                            ) / (Cra + Csv_l)
            elif rho * g[t] * Hu <= P_ra:
                P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l  * (rho * g * Hl 
                            - Pext_l[t]) + Csv_u * rho * g *Hu + Cra * 
                            P_thorax(t)) / (Cra + Csv_u + Csv_l)
        else:
            if P_ra[t] <= rho * g[t] * Hu:
                P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l * (rho * g[t] * Hl 
                            - Pext_l[t]) + Cra * P_thorax[t] 
                            - Csv_l*Pext_l[t]) / Cra
            elif P_ra[t] >= rho * g[t] * Hu:
                P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l * rho * g[t] * Hl + 
                           Csv_u * rho * g[t] * Hu + Cra * P_thorax[t]
                          ) / (Cra + Csv_u)
            elif Pext_l[t] <= P_ra[t]:
                P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l * (rho * g[t] * Hl -
                           Pext_l[t]) + Csv_u * rho * g[t] * Hu + Cra * P_thorax
                          ) / (Cra + Csv_u + Csv_l)
    
 
                

        Psv_u[t] = max(0, P_ra[t] - rho * g[t] * Hu) # eq 15
        Psv_l[t] = max(P_ra,Pext_l) +rho*g[t]*Hl # eq 18

        Psa_u[t] = (Vsa[t] - Vsa0[0] + Csa_l*Pext_l[t] - Csa_l* rho * g[t]*H)/(
                    Csa) # eq 12
        Psa_l[t] = (Vsa[t] - Vsa0[0] + Csa_l*Pext_l[t] - Csa_u* rho * g[t]*H)/(
                    Csa) # eq 13

    # Case II
    elif P_thorax[t] > - dP_RA and P_thorax[t] < rho * g[t] * Hu - dP_RA:

        continue
    # Case III
    elif P_thorax[t] >= rho * g[t] * Hu - dP_RA:
        continue


# def volume_odes(x,t):
#     Qrv = x[0]
#     Qlv = x[1]
#     Qsa_l = x[2]
#     Qsa_    u = x[3]
#     Qsv_l = x[4]
#     Qsv_u = x[5]
#
#     dVstroke_dt = Qrv - Qlv
#     dVsa_udt = Qlv - Qsa_l
#     dVsa_ldt = Qsa_u - Qsv_l
#     dVsv_ldt = Qsa_l - Qsv_u
#     dVsv_udt = Qsv_l -Qrv
#
#     return [dVstroke_dt, dVsa_udt, dVsa_ldt, dVsv_ldt, dVsv_udt]

# initial conditions - wrong need flows not volume ??
# oVstroke = 70 #ml
# artery_vein_split = 0.5
# oVsa_u = artery_vein_split * Hu_factor * Vtotal
# oVsv_u = artery_vein_split * Hu_factor * Vtotal
# oVsa_l = artery_vein_split * Hl_factor * Vtotal
# oVsv_l = artery_vein_split * Hl_factor * Vtotal

# x0 = [oVstroke, oVsa_u, oVsa_l, oVsv_l, oVsv_u]

Qrv = 500
Qlv = 500
Qsa_l = 500
Qsa_u = 500
Qsv_l = 500
Qsv_u = 500
x0 = 6*[Qrv]
test = volume_odes(x0, 0)

# timestep
t = np.linspace(0,3*60)

# SOLVE
x = odeint(volume_odes, x0, t)
breakpoint()
# plot
