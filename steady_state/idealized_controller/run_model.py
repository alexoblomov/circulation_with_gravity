"""
Purpose: simulate model with varying values of G and plot case I (or I or III)
for the purpose of comparison with the centrifuge data
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from parameters import *
# mpl.rcParams['mathtext.fontset'] = 'cm'
from matplotlib.lines import Line2D


G = np.linspace(g_earth, 4*980, 4)

P_thorax = -4 * 1333
P_RA = P_thorax + dP_RA

cases = np.empty((len(Vtotal),len(G)), dtype=np.uint8)

Vd_total_vec = np.empty(len(G))
Q_vec = np.empty(len(G))
F_vec = np.empty(len(G))
Ppa_vec = np.empty(len(G))


sol_F_Vtotal_G = np.empty((len(Vtotal),len(G)))
# sol_Ppa_Vtotal_G = np.empty((len(Vtotal),len(G)))

Hu = Hu_patient
Hl = Hl_patient
dP_RA = P_RA - P_thorax
for j in range(len(Vtotal)):
    for i in range(len(G)):
        if P_thorax <= - dP_RA:
            Vd_total = Vtotal[j] - Cp * (C_RVD / C_LVD) * (dP_RA) \
                       - (Tp * Gs + Csa) * Psa_u_star \
                       - (Tp * Gs_l + Csa_l) \
                       * rho * G[i] * Hu - Cs_l * rho * G[i] * (- Hl)
            Q = Gs * Psa_u_star + rho * G[i] * Hu / Rs_l
            F = Q / (C_RVD * (P_RA - P_thorax))
            Ppv = P_thorax + (C_RVD / C_LVD) * dP_RA # eq 44
            Ppa = Ppv + Q * Rp
            Ppa_also = P_thorax + (C_RVD / C_LVD) * dP_RA + Q*Rp
        
            cases[j, i] = 1
        elif P_thorax > - dP_RA and P_thorax < rho * G[i] * Hu - dP_RA:
            Vd_total = Vtotal[j] - Cp * (C_RVD / C_LVD) * dP_RA \
                      - (Tp * Gs + Csa)* Psa_u_star \
                      - (Tp * Gs_l + Csa_l) * (rho * G[i] * Hu) \
                      - Cs_l * rho * G[i] * (- Hl) \
                      - (Csv_l - Tp * Gs_l) * (P_thorax + dP_RA)
            Psv_l = - rho * G[i] * Hl + P_thorax + dP_RA
            Psv_u = 0
            Psa_l = Psa_u_star + rho * G[i] * (Hu - Hl)
            # Qs_u = Psa_u_star / Rs_u
            # Qs_l = (Psa_l - Psv_l) / Rs_l
            Qs_u = Psa_u_star * 1 / Rs_u # eq 61
            #eq 62
            Qs_l = (Psa_u_star + rho * G[i] * Hu \
                    - P_thorax-dP_RA) * 1/Rs_l
            Q = Qs_u + Qs_l
            F = Q / (C_RVD * (dP_RA))
            F_also = (Gs * Psa_u_star +
                      (1/Rs_l)* (rho*G[i]*Hu-P_thorax - dP_RA)) / (
                      C_RVD*dP_RA)
          

            Ppv = P_thorax + (C_RVD / C_LVD) * dP_RA # eq 66
            Ppa = Ppv + Q * Rp
            Ppa_also = (C_RVD / C_LVD) * dP_RA + Rp*(Gs*Psa_u_star + 1/Rs_l *
                        (rho * G[i] *Hu - P_thorax - dP_RA))
         
            cases[j, i] = 2
        elif P_thorax[j] >= rho * G[i] * Hu - dP_RA:
            
            Vd_total = Vtotal[j] - Cp * (C_RVD / C_LVD) * dP_RA - \
                    (Tp*Gs + Csa_l+Csa_u)*Psa_u_star \
                    - (Tp*Gs + Csa_l - Csv_u) * rho * G[i] * Hu \
                    - (Csa_l+Csv_l) * rho * G[i] * (-Hl) \
                    - (Csv_l+Csv_u - Tp*Gs)* (P_thorax + dP_RA)
            Psv_l = P_thorax + dP_RA + rho * G[i] * (- Hl)
            Psv_u = P_thorax + dP_RA - rho * G[i] * Hu
            Psa_l = Psa_u_star + rho * G[i] * (Hu - Hl)
            Qs_u = (Psa_u_star - Psv_u) / Rs_u
            Qs_u_also = (Psa_u_star + rho*G[i]*Hu - P_thorax - dP_RA)/Rs_u   
            Qs_l = (Psa_l - Psv_l) / Rs_l
            Qs_l_also = (Psa_u_star + rho*G[i]*Hu - P_thorax - dP_RA)/Rs_l          
            Q = Qs_u + Qs_l
            F = Q / (C_RVD * (dP_RA))
            Ppv = P_thorax + (C_RVD / C_LVD) * dP_RA
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
    sol_Vd_Vtotal_G[j,:] = Vd_total_vec
    sol_Q_Vtotal_G[j,:] = Q_vec
    sol_F_Vtotal_G[j,:] = F_vec
    sol_Ppa_Vtotal_G[j,:] = Ppa_vec

#conversions:
G = G / 100 / (g_earth / 100)
sol_Q_Vtotal_G = sol_Q_Vtotal_G * 60 / 1000
sol_F_Vtotal_G = sol_F_Vtotal_G * 60
sol_Ppa_Vtotal_G = sol_Ppa_Vtotal_G / 1333
sol_Vd_Vtotal_G = sol_Vd_Vtotal_G / 1000

Vtotal_string = r'$V_\mathrm{total}$'
L = r'$\mathrm{L}$'

Vtotal_titles = [Vtotal_string + " $=$ " + str(vtot/1000) + r" $\mathrm{L}$" for vtot in Vtotal]
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
# Define shades of red and yellow
red_shades = np.linspace(1, 0, len(Vtotal) + 1)  # From 1 (bright red) to 0 (dark red)
yellow_shades = np.linspace(0.8, 0, len(Vtotal) + 1)  # From 1 (bright yellow) to 0 (dark yellow)

# Convert shades to colors in the colormap
cmap = plt.get_cmap(color_map)  # Get the colormap
line_colors = [cmap(shade) for shade in np.concatenate([red_shades, yellow_shades])]
line_style = ['--', ':', '-']


fig, ax = plt.subplots() 
for n, plt_title in enumerate(Vtotal_titles):
    case_1_indices = np.where(cases[n,:] == 1)
    case_2_indices = np.where(cases[n,:] == 2)
    case_3_indices = np.where(cases[n,:] == 3)
    ax.plot(G[case_1_indices], sol_Vd_Vtotal_G[n, case_1_indices[0]], color=line_colors[len(Vtotal) - n], linestyle = line_style[0])
    ax.plot(G[case_2_indices], sol_Vd_Vtotal_G[n, case_2_indices[0]], color=line_colors[len(Vtotal) - n], linestyle = line_style[1])
    ax.plot(G[case_3_indices], sol_Vd_Vtotal_G[n, case_3_indices[0]], color=line_colors[len(Vtotal) - n], linestyle = line_style[2], label = Vtotal_string + r" $=$ " + fr'{Vtotal[n]/1000}' + L)
    ax.set_title(r'$\mathrm{Reserve}$ $\mathrm{Volume}$ $\mathrm{v.}$ $\mathrm{Acceleration}$ $\mathrm{Varying}$ $\mathrm{Total}$ $\mathrm{Volume}$')
    # ax.get_legend().remove()
    ax.set_xlabel(r"$g$ $\mathrm{multiple}$")
    ax.set_ylabel(r"$V_{\mathrm{total}}^0$ $(\mathrm{L})$")
    
color_legend = ax.legend(ncol=2, loc = 'upper right')
legend_handles = [
    Line2D([], [], linestyle='dashed', color='maroon', label=r'$\mathrm{Case}$ $1$'),
    Line2D([], [], linestyle='dotted', color='maroon', label=r'$\mathrm{Case}$ $2$'),
    Line2D([], [], linestyle='solid', color='maroon', label=r'$\mathrm{Case}$ $3$')
]
linestyle_legend = ax.legend(handles = legend_handles, loc = 'center right')
ax.add_artist(color_legend)
plt.grid(True)
plt.savefig("varyVtotal_V0_G", bbox_inches='tight', dpi=300)

plt.figure()
plt.title(r"$+\mathrm{Gz}$ $\mathrm{Tolerance}$ $\mathrm{Varying}$ $\mathrm{Total}$ $\mathrm{Volume}$")
plt.xlabel(r"$V_\mathrm{total}$ $\mathrm{(L)}$")
plt.ylabel(r"$+\mathrm{Gz}$ $\mathrm{Tolerance}$ $(g$ $\mathrm{multiple})$")

for i in range(len(Vtotal)):
    min_Vd_total = np.nanmin(sol_Vd_Vtotal_G[i, :])
    if np.isnan(min_Vd_total):
        g_intercept = G[np.isnan(sol_Vd_Vtotal_G[i, :])][-1]
    else:
        min_Vd_total_index = np.nanargmin(sol_Vd_Vtotal_G[i, :])
        g_intercept = G[min_Vd_total_index]
    plt.plot(Vtotal[i] / 1000, g_intercept, 'o-', markersize=6, color=line_colors[len(Vtotal) - i])

plt.grid(True)

# Specify the file path for saving the figure
plt.savefig("varyVtotal_gtol_plot.png", bbox_inches='tight', dpi=300)

