"""
Purpose: simulate model with varying values of G and plot output variables
of interest
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc
from parameters import *
mpl.rcParams['mathtext.fontset'] = 'cm'

# Enable LaTeX rendering
mpl.rcParams['text.usetex'] = True

# Set the font family to serif
mpl.rcParams['font.family'] = 'serif'



G = np.linspace(g_earth,3*980, 3000)

P_thorax = np.linspace(- 4 * 1333, 24 * 1333, 8)
# Define shades of red and yellow
red_shades = np.linspace(1, 0, len(P_thorax))  # From 1 (bright red) to 0 (dark red)
yellow_shades = np.linspace(0.8, 0, len(P_thorax))  # From 1 (bright yellow) to 0 (dark yellow)

# Convert shades to colors in the colormap
cmap = plt.cm.get_cmap(color_map)  # Get the colormap
line_colors = [cmap(shade) for shade in np.concatenate([red_shades, yellow_shades])]

P_RA = P_thorax + dP_RA

cases = np.empty((len(P_thorax),len(G)), dtype=np.uint8)

Vd_total_vec = np.empty(len(G))
Q_vec = np.empty(len(G))
F_vec = np.empty(len(G))
Ppa_vec = np.empty(len(G))

sol_Vd_Pthorax_G = np.empty((len(P_thorax),len(G)))
sol_Q_Pthorax_G = np.empty((len(P_thorax),len(G)))
sol_F_Pthorax_G = np.empty((len(P_thorax),len(G)))
sol_Ppa_Pthorax_G = np.empty((len(P_thorax),len(G)))

Hu = Hu_patient
Hl = Hl_patient

for j in range(len(P_thorax)):
    dP_RA = P_RA[j] - P_thorax[j]
    for i in range(len(G)):
        if P_thorax[j] <= - dP_RA:
            Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * (dP_RA) \
                       - (Tp * Gs + Csa) * Psa_u_star \
                       - (Tp * Gs_l + Csa_l) \
                       * rho * G[i] * Hu - Cs_l * rho * G[i] * (- Hl)
            Q = Gs * Psa_u_star + rho * G[i] * Hu / Rs_l
            F = Q / (C_RVD * (P_RA[j] - P_thorax[j]))
            Ppv = P_thorax[j] + (C_RVD / C_LVD) * dP_RA # eq 44
            Ppa = Ppv + Q * Rp
            Ppa_also = P_thorax[j] + (C_RVD / C_LVD) * dP_RA + Q*Rp
            # assert (nabs(Ppa - Ppa_also) < 0.01)
            cases[j, i] = 1
        elif P_thorax[j] > - dP_RA and P_thorax[j] < rho * G[i] * Hu - dP_RA:
            Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA \
                      - (Tp * Gs + Csa)* Psa_u_star \
                      - (Tp * Gs_l + Csa_l) * (rho * G[i] * Hu) \
                      - Cs_l * rho * G[i] * (- Hl) \
                      - (Csv_l - Tp * Gs_l) * (P_thorax[j] + dP_RA)
            Psv_l = - rho * G[i] * Hl + P_thorax[j] + dP_RA
            Psv_u = 0
            Psa_l = Psa_u_star + rho * G[i] * (Hu - Hl)
            
            Qs_u = Psa_u_star / Rs_u # eq 61
            #eq 62
            Qs_l = (Psa_u_star + rho * G[i] * Hu \
                    - P_thorax[j]-dP_RA) * 1/Rs_l
            #Qs_l = (Psa_l - Psv_l) / Rs_l
            Q = Qs_u + Qs_l
            F = Q / (C_RVD * (dP_RA))
            F_also = (Gs * Psa_u_star +
                      (1/Rs_l)* (rho*G[i]*Hu-P_thorax[j] - dP_RA)) / (
                      C_RVD*dP_RA)

            Ppv = P_thorax[j] + (C_RVD / C_LVD) * dP_RA # eq 66
            Ppa = Ppv + Q * Rp
            Ppa_also = (C_RVD / C_LVD) * dP_RA + Rp*(Gs*Psa_u_star + 1/Rs_l *
                        (rho * G[i] *Hu - P_thorax[j] - dP_RA))
            cases[j, i] = 2
        elif P_thorax[j] >= rho * G[i] * Hu - dP_RA:
            # print("entered case III")
            Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA - \
                    (Tp*Gs + Csa_l+Csa_u)*Psa_u_star \
                    - (Tp*Gs + Csa_l - Csv_u) * rho * G[i] * Hu \
                    - (Csa_l+Csv_l) * rho * G[i] * (-Hl) \
                    - (Csv_l+Csv_u - Tp*Gs)* (P_thorax[j] + dP_RA)

            Psv_l = P_thorax[j] + dP_RA + rho * G[i] * (- Hl)
            Psv_u = P_thorax[j] + dP_RA - rho * G[i] * Hu
            Psa_l = Psa_u_star + rho * G[i] * (Hu - Hl)
            Qs_u = (Psa_u_star - Psv_u) / Rs_u
            Qs_u_also = (Psa_u_star + rho*G[i]*Hu - P_thorax[j] - dP_RA)/Rs_u
            Qs_l = (Psa_l - Psv_l) / Rs_l
            Qs_l_also = (Psa_u_star + rho*G[i]*Hu - P_thorax[j] - dP_RA)/Rs_l
            Q = Qs_u + Qs_l
            F = Q / (C_RVD * (dP_RA))
            Ppv = P_thorax[j] + (C_RVD / C_LVD) * dP_RA
            Ppa = Ppv + Q * Rp
            cases[j, i] = 3
        if Vd_total >= 0:
            Vd_total_vec[i] = Vd_total
            Q_vec[i] = Q
            F_vec[i] = F
            Ppa_vec[i] = Ppa
        else:
            Vd_total_vec[i] = float('nan')
            Q_vec[i] = float('nan')
            F_vec[i] = float('nan')
            Ppa_vec[i] = float('nan')
    sol_Vd_Pthorax_G[j,:] = Vd_total_vec
    sol_Q_Pthorax_G[j,:] = Q_vec
    sol_F_Pthorax_G[j,:] = F_vec
    sol_Ppa_Pthorax_G[j,:] = Ppa_vec

#conversions:
G = G / 100 / (g_earth / 100)
sol_Q_Pthorax_G = sol_Q_Pthorax_G * 60 / 1000
sol_F_Pthorax_G = sol_F_Pthorax_G * 60
sol_Ppa_Pthorax_G = sol_Ppa_Pthorax_G / 1333
sol_Vd_Pthorax_G = sol_Vd_Pthorax_G / 1000

P_thorax_string = r'$P_\mathrm{thorax}$'
mmHg = r'$\mathrm{mmHg}$'
# Plotting G against Vd for different P_thorax values
fig, ax = plt.subplots()
for i in range(len(P_thorax)):
    ax.plot(G, sol_Vd_Pthorax_G[i, :], label=fr'{P_thorax_string} = {P_thorax[i]/1333} {mmHg}', color=line_colors[i])
ax.set_xlabel(r'$+\mathrm{Gz}$ $(g$ $\mathrm{multiples})$')
ax.set_ylabel(r"$V_{\mathrm{total}}^0$ $(\mathrm{L}$)")
ax.set_title(r'$\mathrm{Reserve}$ $\mathrm{Volume}$ $\mathrm{v.}$ $\mathrm{Acceleration}$ $\mathrm{Varying}$ $\mathrm{Intrathoracic}$ $\mathrm{Pressures}$')
ax.legend()
ax.set_xlim(1, 2.3)
ax.tick_params(axis='both', labelsize=10)  # Set tick label font size
plt.rc('font', family='serif')  # Set font family to use LaTeX font
plt.grid(True)
plt.savefig('figures/varyPthorax_Vd_G', bbox_inches='tight', dpi=300)

# Plotting G against F for different P_thorax values
plt.figure()
for i in range(len(P_thorax)):
    plt.plot(G, sol_F_Pthorax_G[i, :], label=fr'{P_thorax_string} = {P_thorax[i]/1333} {mmHg}', color=line_colors[i])
plt.xlabel(r'$+\mathrm{Gz}$ $(g$ $\mathrm{multiples})$')
plt.ylabel(r'$F$ $\mathrm{(beats/min)}$')
plt.title(r'$\mathrm{Heart}$ $\mathrm{Rate}$ $\mathrm{v.}$ $\mathrm{Acceleration}$ $\mathrm{Varying}$ $\mathrm{Intrathoracic}$ $\mathrm{Pressures}$')
plt.legend()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.xlim(1, 2.3)
plt.grid(True)
#plt.show()
plt.savefig('figures/varyPthorax_F_G', bbox_inches='tight', dpi=300)


plt.figure()
plt.title(r"$+\mathrm{Gz}$ $\mathrm{Tolerance}$ $\mathrm{Varying}$ $\mathrm{Intrathoracic}$ $\mathrm{Pressure}$")
plt.xlabel(r"$P_\mathrm{thorax}$ $\mathrm{(mmHg)}$")
plt.ylabel(r"$+\mathrm{Gz}$ $\mathrm{Tolerance}$ $(g$ $\mathrm{multiples})$")
for i in range(len(P_thorax)):
    min_Vd_total = np.nanmin(sol_Vd_Pthorax_G[i, :])
    if np.isnan(min_Vd_total):
        g_intercept = G[np.isnan(sol_Vd_Pthorax_G[i, :])][-1]
    else:
        min_Vd_total_index = np.nanargmin(sol_Vd_Pthorax_G[i, :])
        g_intercept = G[min_Vd_total_index]
    plt.plot(P_thorax[i] / 1333, g_intercept, 'o-', markersize=6, color=line_colors[i])
# Enable LaTeX rendering
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.grid(True)
plt.savefig('figures/varyPthorax_gtol_plot', bbox_inches='tight', dpi=300)
