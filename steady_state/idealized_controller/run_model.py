import numpy as np
import matplotlib.pyplot as plt
from parameters import *
from matplotlib import rc
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'


# Fix P_thorax to a single value
P_thorax = -4 * 1333
P_RA = P_thorax + dP_RA

# Discretization
G = np.linspace(g_earth, 10 * 980, 100)

# Initialize array to store G tolerance values
G_tolerance = np.zeros(len(G))

# Initialize Vd_total vector
Vd_total_vec = np.empty(len(G))
F_vec = np.empty(len(G))

for k in range(len(G)):
    if P_thorax <= -dP_RA:
        Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * (dP_RA) \
                    - (Tp * Gs + Csa_u + Csa_l) * Psa_u_star \
                    - (Tp * Gs_l + Csa_l) * rho * G[k] * Hu_patient - (Csa_l + Csv_l) * rho * G[k] * (-Hl_patient)
        
        Q = Gs * Psa_u_star + rho * G[k] * Hu_patient / Rs_l
        F = Q / (C_RVD * (P_RA - P_thorax))
        
        # print(f'G is  {G[k]} Vd_total is  {Vd_total}')
    elif P_thorax > -dP_RA and P_thorax < rho * G[k] * Hu_patient - dP_RA:
        Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA \
                    - (Tp * Gs + Csa_u + Csa_l) * Psa_u_star \
                    - (Tp * Gs_l + Csa_l) * (rho * G[k] * Hu_patient) \
                    - Cs_l * rho * G[k] * (-Hl_patient) \
                    - (Csv_l - Tp * Gs_l) * (P_thorax + dP_RA)
        
        Psv_l = - rho * G[k] * Hl_patient + P_thorax + dP_RA
        Psv_u = 0
        Psa_l = Psa_u_star + rho * G[k] * (Hu_patient - Hl_patient)
        
        Qs_u = Psa_u_star / Rs_u # eq 61
        #eq 62
        Qs_l = (Psa_u_star + rho * G[k] * Hu_patient \
                - P_thorax-dP_RA) * 1/Rs_l
        #Qs_l = (Psa_l - Psv_l) / Rs_l
        Q = Qs_u + Qs_l
        
        F = Q / (C_RVD * (dP_RA))

    elif P_thorax >= rho * G[k] * Hu_patient - dP_RA:
        Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA \
                    - (Tp * Gs + Csa_l + Csa_u) * Psa_u_star \
                    - (Tp * Gs + Csa_l - Csv_u) * rho * G[k] * Hu_patient \
                    - (Csa_l + Csv_l) * rho * G[k] * (-Hl_patient) \
                    - (Csv_l + Csv_u - Tp * Gs) * (P_thorax + dP_RA)
        
        Psv_l = P_thorax + dP_RA + rho * G[k] * (- Hl_patient)
        Psv_u = P_thorax + dP_RA - rho * G[k] * Hu_patient
        Psa_l = Psa_u_star + rho * G[k] * (Hu_patient- Hl_patient)
        Qs_u = (Psa_u_star - Psv_u) / Rs_u
        Qs_l = (Psa_l - Psv_l) / Rs_l
        Q = Qs_u + Qs_l

        F = Q / (C_RVD * (dP_RA))

    if Vd_total >= 0:
        Vd_total_vec[k] = Vd_total
        F_vec[k] = F
    else:
        Vd_total_vec[k] = float('nan')
        F_vec[k] = float('nan')

# Find G tolerance value
zero_indices = np.where(np.isnan(Vd_total_vec))[0]
if len(zero_indices) > 0:
    G_tolerance = G[zero_indices[0]] / 100 / 9.8
else:
    G_tolerance = G[np.argmin(np.abs(Vd_total_vec))] / 100 / 9.8

print(f'g tolerance is {G_tolerance}')

#conversions:
G = G / 100 / (g_earth / 100)
sol_F_Vtotal_G = F_vec * 60
sol_Vd_Vtotal_G = Vd_total_vec / 1000

breakpoint()
plt.figure()


plt.plot(G, G_tolerance)
plt.save("g_vs_g_tolerance.png")