import numpy as np
import matplotlib.pyplot as plt
from parameters import *
from parameters import *
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'
# Fix P_thorax to a single value
P_thorax = -4 * 1333

# Discretization
dx = 100
Csa_u = np.linspace(0.01 * Csa_u, 8 * Csa_u, dx)
Csa_l = np.linspace(0.01 * Csa_l, 8* Csa_l, dx)

Hu = Hu_patient
Hl = Hl_patient

# Initialize arrays to store G tolerance values
G_tolerance = np.zeros((len(Csa_u), len(Csa_l)))

# Loop over Csa_u and Csa_l values
for i in range(len(Csa_u)):
    for j in range(len(Csa_l)):
        G = np.linspace(g_earth, 10 * 980, 3000)

        # Initialize Vd_total vector
        Vd_total_vec = np.empty(len(G))

        for k in range(len(G)):
            if P_thorax <= -dP_RA:
                Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * (dP_RA) \
                           - (Tp * Gs + Csa_u[i]+Csa_l[j]) * Psa_u_star \
                           - (Tp * Gs_l + Csa_l[j]) * rho * G[k] * Hu - (Csa_l[j] +Csv_l) * rho * G[k] * (-Hl)
            elif P_thorax > -dP_RA and P_thorax < rho * G[k] * Hu - dP_RA:
                Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA \
                           - (Tp * Gs + Csa_u[i]+Csa_l[j]) * Psa_u_star \
                           - (Tp * Gs_l + Csa_l[j]) * (rho * G[k] * Hu) \
                           - Cs_l * rho * G[k] * (-Hl) \
                           - (Csv_l - Tp * Gs_l) * (P_thorax + dP_RA)
            elif P_thorax >= rho * G[k] * Hu - dP_RA:
                Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA \
                           - (Tp * Gs + Csa_l[j] + Csa_u[i]) * Psa_u_star \
                           - (Tp * Gs_l + Csa_l[j] - Csv_u) * rho * G[k] * Hu \
                           - (Csa_l[j] + Csv_l) * rho * G[k] * (-Hl) \
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
heatmap = ax.imshow(G_tolerance, cmap='inferno', aspect='auto', origin='lower')

x_tick_indices = np.linspace(0, len(Csa_l) - 1, 5, dtype=int)
y_tick_indices = np.linspace(0, len(Csa_u) - 1, 5, dtype=int)
ax.set_xticks(x_tick_indices)
ax.set_yticks(y_tick_indices)
ax.set_xticklabels([f'{Csa_l[i]:.3f}' for i in x_tick_indices])
ax.set_yticklabels([f'{Csa_u[i]:.3f}' for i in y_tick_indices])

ax.set_xlabel(r'$C_{\mathrm{sa}}^{\mathrm{l}}$ $\mathrm{(ml/mmHg)}$')
ax.set_ylabel(r'$C_{\mathrm{sa}}^{\mathrm{u}}$ $\mathrm{(ml/mmHg)}$')
ax.set_title(r'$\mathrm{+Gz}$ $\mathrm{Tolerance}$ $\mathrm{Varying}$ $\mathrm{Arterial}$ $\mathrm{Compliances}$')
plt.colorbar(heatmap, label=r'$g$ $\mathrm{Multiples}$')
plt.grid(False)
plt.savefig('vary_Csa_gtol_heatmap')

# Plotting the heatmap
fig, ax = plt.subplots()
heatmap = ax.imshow(G_tolerance, cmap='jet', aspect='auto', origin='lower')

ax.set_xlabel('Csa_l (ml/mmHg)')
ax.set_ylabel('Csa_u (ml/mmHg)')
ax.set_title('G Tolerance Heatmap')
ax.tick_params(axis='both', labelsize=8)  # Set tick label font size
plt.colorbar(heatmap, label='G Multiples')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.grid(False)
#plt.show()
