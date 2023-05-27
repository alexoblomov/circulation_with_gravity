import numpy as np
import matplotlib.pyplot as plt
from parameters import *

# Fix P_thorax to a single value
P_thorax = -4 * 1333

# Discretization
dx = 100
Hu_values = np.linspace(0.5 * Hu_patient, 3 * Hu_patient, dx)
Hl_values = np.linspace(0.5 * Hl_patient, 1 * Hl_patient, dx)

Csa_u = Csa_u
Csa_l = Csa_l

# Initialize arrays to store G tolerance values
G_tolerance = np.zeros((len(Hu_values), len(Hl_values)))

# Loop over Hu and Hl values
for i in range(len(Hu_values)):
    for j in range(len(Hl_values)):
        G = np.linspace(g_earth, 10 * 980, 3000)

        # Initialize Vd_total vector
        Vd_total_vec = np.empty(len(G))

        for k in range(len(G)):
            if P_thorax <= -dP_RA:
                Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * (dP_RA) \
                           - (Tp * Gs + Csa_u + Csa_l) * Psa_u_star \
                           - (Tp * Gs_l + Csa_l) * rho * G[k] * Hu_values[i] - (Csa_l + Csv_l) * rho * G[k] * (-Hl_values[j])
            elif P_thorax > -dP_RA and P_thorax < rho * G[k] * Hu_values[i] - dP_RA:
                Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA \
                           - (Tp * Gs + Csa_u + Csa_l) * Psa_u_star \
                           - (Tp * Gs_l + Csa_l) * (rho * G[k] * Hu_values[i]) \
                           - Cs_l * rho * G[k] * (-Hl_values[j]) \
                           - (Csv_l - Tp * Gs_l) * (P_thorax + dP_RA)
            elif P_thorax >= rho * G[k] * Hu_values[i] - dP_RA:
                Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA \
                           - (Tp * Gs + Csa_l + Csa_u) * Psa_u_star \
                           - (Tp * Gs_l + Csa_l - Csv_u) * rho * G[k] * Hu_values[i] \
                           - (Csa_l + Csv_l) * rho * G[k] * (-Hl_values[j]) \
                           - (Csv_l + Csv_u - Tp * Gs) * (P_thorax + dP_RA)

            if Vd_total >= 0:
                Vd_total_vec[k] = Vd_total
            else:
                Vd_total_vec[k] = float('nan')

        # Find G tolerance value
        zero_indices = np.where(Vd_total_vec <= 0)[0]
        if len(zero_indices) > 0:
            G_tolerance[i, j] = G[zero_indices[0]] / 1000
        else:
            G_tolerance[i, j] = G[np.argmin(np.abs(Vd_total_vec))] / 1000

# Plotting the heatmap
fig, ax = plt.subplots()
heatmap = ax.imshow(G_tolerance, cmap='jet', aspect='auto', origin='lower')

ax.set_xlabel('Hl (cm)')
ax.set_ylabel('Hu (cm)')
ax.set_title('G Tolerance Heatmap')
plt.colorbar(heatmap, label='G Multiples')
plt.grid(True)
plt.show()
