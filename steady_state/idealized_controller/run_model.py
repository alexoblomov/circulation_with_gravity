import numpy as np
import matplotlib.pyplot as plt
from parameters import *
from matplotlib import rc
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'
# Fix P_thorax to a single value
P_thorax = -4 * 1333

# Discretization
G = np.linspace(g_earth, 10 * 980, 10)

# Initialize array to store G tolerance values
G_tolerance = np.zeros(len(G))

# Initialize Vd_total vector
Vd_total_vec = np.empty(len(G))

for k in range(len(G)):
    if P_thorax <= -dP_RA:
        Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * (dP_RA) \
                    - (Tp * Gs + Csa_u + Csa_l) * Psa_u_star \
                    - (Tp * Gs_l + Csa_l) * rho * G[k] * Hu_patient - (Csa_l + Csv_l) * rho * G[k] * (-Hl_patient)
        print(f'G is  {G[k]} Vd_total is  {Vd_total}')
    elif P_thorax > -dP_RA and P_thorax < rho * G[k] * Hu_patient - dP_RA:
        Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA \
                    - (Tp * Gs + Csa_u + Csa_l) * Psa_u_star \
                    - (Tp * Gs_l + Csa_l) * (rho * G[k] * Hu_patient) \
                    - Cs_l * rho * G[k] * (-Hl_patient) \
                    - (Csv_l - Tp * Gs_l) * (P_thorax + dP_RA)
    elif P_thorax >= rho * G[k] * Hu_patient - dP_RA:
        Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA \
                    - (Tp * Gs + Csa_l + Csa_u) * Psa_u_star \
                    - (Tp * Gs + Csa_l - Csv_u) * rho * G[k] * Hu_patient \
                    - (Csa_l + Csv_l) * rho * G[k] * (-Hl_patient) \
                    - (Csv_l + Csv_u - Tp * Gs) * (P_thorax + dP_RA)

    if Vd_total >= 0:
        Vd_total_vec[k] = Vd_total
    else:
        Vd_total_vec[k] = float('nan')

# Find G tolerance value
zero_indices = np.where(np.isnan(Vd_total_vec))[0]
if len(zero_indices) > 0:
    G_tolerance = G[zero_indices[0]] / 100 / 9.8
else:
    G_tolerance = G[np.argmin(np.abs(Vd_total_vec))] / 100 / 9.8

print(f'g tolerance is {G_tolerance}')
