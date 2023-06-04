"""
Purpose: simulate model with varying values of G and plot output variables
of interest
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from parameters import *
mpl.rcParams['mathtext.fontset'] = 'cm'

G = np.linspace(g_earth,7*980, 3000)

Vtotal = np.linspace(3, 8, 6)
Vtotal = Vtotal*1000
P_thorax = -4 * 1333
P_RA = P_thorax + dP_RA

cases = np.empty((len(Vtotal),len(G)), dtype=np.uint8)

Vd_total_vec = np.empty(len(G))
Q_vec = np.empty(len(G))
F_vec = np.empty(len(G))
Ppa_vec = np.empty(len(G))

sol_Vd_Vtotal_G = np.empty((len(Vtotal),len(G)))
sol_Q_Vtotal_G = np.empty((len(Vtotal),len(G)))
sol_F_Vtotal_G = np.empty((len(Vtotal),len(G)))
sol_Ppa_Vtotal_G = np.empty((len(Vtotal),len(G)))

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
            # print("entered case III")
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
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
# Define shades of red and yellow
red_shades = np.linspace(1, 0, len(Vtotal))  # From 1 (bright red) to 0 (dark red)
yellow_shades = np.linspace(1, 0, len(Vtotal))  # From 1 (bright yellow) to 0 (dark yellow)

# Convert shades to colors in the colormap
cmap = plt.cm.get_cmap(color_map)  # Get the colormap
line_colors = [cmap(shade) for shade in np.concatenate([red_shades, yellow_shades])]

plt.figure()
plt.title(r'$\mathrm{Heart}$ $\mathrm{Rate}$ $\mathrm{v.}$ $\mathrm{Acceleration}$ $\mathrm{Varying}$ $\mathrm{Total}$ $\mathrm{Volume}$')
plt.xlabel(r"$g$ $\mathrm{multiples}$")
plt.ylabel(r"$F$ $(\mathrm{beats/min})$")
for n, plt_title in enumerate(Vtotal_titles):
    # Create a copy of the data to avoid overwriting
    F_values = np.copy(sol_F_Vtotal_G[n, :])
    plt.plot(G, F_values, color=line_colors[n])
    
plt.legend(Vtotal_titles)
plt.grid(True)
plt.savefig("figures/varyVtotal_F_G", bbox_inches='tight', dpi=300)
plt.figure(figsize=(15, 12))
plt.subplots_adjust(hspace=0.5)
plt.suptitle(r"$\mathrm{Reserve}$ $\mathrm{volume}$ v. g", fontsize=18, y=0.95)



plt.figure()
for n, plt_title in enumerate(Vtotal_titles):
    # add a new subplot iteratively
    # filter df and plot ticker on the new subplot axis
    idx_case_1 = cases == 1
    idx_case_2 = cases == 2
    idx_case_3 = cases == 3
    plt.plot(G, sol_Vd_Vtotal_G[n, :], color=line_colors[n])
    plt.title(r'$\mathrm{Reserve}$ $\mathrm{Volume}$ $\mathrm{v.}$ $\mathrm{Acceleration}$ $\mathrm{Varying}$ $\mathrm{Total}$ $\mathrm{Volume}$')
    # ax.get_legend().remove()
    plt.xlabel(r"$g$ $\mathrm{multiples}$")
    plt.ylabel(r"$V_{\mathrm{total}}^0$ $(\mathrm{L})$")
    
plt.legend(Vtotal_titles)
plt.grid(True)
plt.savefig("figures/varyVtotal_V0_G", bbox_inches='tight', dpi=300)

plt.figure()
plt.title(r"$+\mathrm{Gz}$ $\mathrm{Tolerance}$ $\mathrm{Varying}$ $\mathrm{Total}$ $\mathrm{Volume}$")
plt.xlabel(r"$V_\mathrm{total}$ $\mathrm{(L)}$")
plt.ylabel(r"$+\mathrm{Gz}$ $\mathrm{Tolerance}$ $(g$ $\mathrm{multiples})$")

for i in range(len(Vtotal)):
    min_Vd_total = np.nanmin(sol_Vd_Vtotal_G[i, :])
    if np.isnan(min_Vd_total):
        g_intercept = G[np.isnan(sol_Vd_Vtotal_G[i, :])][-1]
    else:
        min_Vd_total_index = np.nanargmin(sol_Vd_Vtotal_G[i, :])
        g_intercept = G[min_Vd_total_index]
    plt.plot(Vtotal[i] / 1000, g_intercept, 'o-', markersize=6, color=line_colors[i])
# Enable LaTeX rendering
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.grid(True)
# Specify the file path for saving the figure
plt.savefig("figures/varyVtotal_gtol_plot.png", bbox_inches='tight', dpi=300)
