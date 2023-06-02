import numpy as np
import matplotlib.pyplot as plt
from parameters import *
from matplotlib import rc
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'
# Fix P_thorax to a single value
P_thorax = -4 * 1333

# Discretization
dx = 90
Hu_values = np.linspace(20, 75, dx)
Hl_values = np.linspace(-25, -55, dx)

Csa_u = Csa_u
Csa_l = Csa_l

# Initialize array to store G tolerance values
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
plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
heatmap = plt.imshow(G_tolerance, cmap=color_map, aspect='auto', origin='lower')

x_tick_indices = np.linspace(0, len(Hl_values) - 1, 5, dtype=int)
y_tick_indices = np.linspace(0, len(Hu_values) - 1, 5, dtype=int)
x_tick_labels = [f'{Hl_values[i]:.2f}' for i in x_tick_indices]
y_tick_labels = [f'{Hu_values[i]:.2f}' for i in y_tick_indices]

plt.xticks(x_tick_indices, [float(label) for label in x_tick_labels])
plt.yticks(y_tick_indices, [float(label) for label in y_tick_labels])

plt.xlabel(r'$H_{\mathrm{l}}$ $\mathrm{(cm)}$')
plt.ylabel(r'$H_{\mathrm{u}}$ $\mathrm{(cm)}$')
plt.title(r'$\mathrm{+Gz}$ $\mathrm{Tolerance}$ $\mathrm{Varying}$ $\mathrm{Height}$ $\mathrm{Ratio}$')
plt.tick_params(axis='both')  # Set tick label font size
cbar = plt.colorbar(heatmap, label=r"$g$ $\mathrm{multiples}$")


plt.grid(False)
plt.savefig('figures/varyH_gtol', bbox_inches='tight', dpi=300)