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
from matplotlib.lines import Line2D
from quadratic import find_zeros
# Enable LaTeX rendering
mpl.rcParams['text.usetex'] = True

# Set the font family to serif
mpl.rcParams['font.family'] = 'serif'

G = np.linspace(g_earth,7*980, 300)

P_thorax = np.linspace(- 4 * 1333, 24 * 1333, 8)
# Define shades of red and yellow
red_shades = np.linspace(1, 0, len(P_thorax) + 1)  # From 1 (bright red) to 0 (dark red)
yellow_shades = np.linspace(0.8, 0, len(P_thorax) + 1)  # From 1 (bright yellow) to 0 (dark yellow)

# Convert shades to colors in the colormap
cmap = plt.cm.get_cmap(color_map)  # Get the colormap
line_colors = [cmap(shade) for shade in np.concatenate([red_shades, yellow_shades])]

P_RA = P_thorax + dP_RA

cases = np.empty((len(P_thorax),len(G)), dtype=np.uint8)

Vd_total_vec = np.empty(len(G))
Q_vec = np.empty(len(G))
F_vec = np.empty(len(G))
Ppa_vec = np.empty(len(G))
Gs_l_vec = np.empty(len(G))

sol_Vd_Pthorax_G = np.empty((len(P_thorax),len(G)))
sol_Q_Pthorax_G = np.empty((len(P_thorax),len(G)))
sol_F_Pthorax_G = np.empty((len(P_thorax),len(G)))
sol_Ppa_Pthorax_G = np.empty((len(P_thorax),len(G)))
sol_Gsl_Pthorax_G = np.empty((len(P_thorax),len(G)))

Hu = Hu_patient
Hl = Hl_patient
#compute Normal Reserve Volume at each P_thorax value:
Vd_total_normal_vec = []
for j in range(len(P_thorax)):
    if P_thorax[j] <= -dP_RA:
        Vd_total_normal = Vtotal - Cp * (C_RVD / C_LVD) * (dP_RA) \
                       - (Tp * (Gs_l_normal +Gs_u) + Csa) * Psa_u_star \
                       - (Tp * Gs_l_normal + Csa_l) \
                       * rho * g_earth * Hu - Cs_l * rho * g_earth * (- Hl)
    elif P_thorax[j] > - dP_RA and P_thorax[j] < rho * g_earth * Hu - dP_RA:
        Vd_total_normal = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA \
                      - (Tp * (Gs_l_normal + Gs_u) + Csa)* Psa_u_star \
                      - (Tp * Gs_l_normal + Csa_l) * (rho * g_earth * Hu) \
                      - Cs_l * rho * g_earth * (- Hl) \
                      - (Csv_l - Tp * Gs_l_normal) * (P_thorax[j] + dP_RA)
    elif P_thorax[j] >= rho * g_earth * Hu - dP_RA:
            # print("entered case III")
        Vd_total_normal = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA - \
                    (Tp*(Gs_l_normal +Gs_u) + Csa)*Psa_u_star \
                    - (Tp*(Gs_l_normal +Gs_u) + Csa_l - Csv_u) * rho * g_earth * Hu \
                    - (Cs_l) * rho * g_earth * (-Hl) \
                    - (Csv_l+Csv_u - Tp*(Gs_l_normal +Gs_u))* (P_thorax[j] + dP_RA)  
    Vd_total_normal_vec.append(Vd_total_normal)
for j in range(len(P_thorax)):
    dP_RA = P_RA[j] - P_thorax[j]
    for i in range(len(G)):#case 1
        if P_thorax[j] <= - dP_RA:
            a = (Tp * Gs_l_normal / Vd_total_normal_vec[j]**2) * (Psa_u_star + rho * G[i] * Hu)
            b = 1
            c = -Vtotal + Cp*(C_RVD/C_LVD)*dP_RA + Tp * Gs_u * Psa_u_star + Csa*(Psa_u_star) \
                + Csa_l*rho*G[i]*Hu + Cs_l * rho *G[i]*(-Hl)
            Vd_total = find_zeros(a,b,c)
            Gs_l = Gs_l_normal*(Vd_total**2/Vd_total_normal_vec[j]**2)
            Rs_l = 1 / Gs_l
            Q = (Gs_u + Gs_l) * Psa_u_star + rho * G[i] * Hu / Rs_l
            F = Q / (C_RVD * (P_RA[j] - P_thorax[j]))
            Ppv = P_thorax[j] + (C_RVD / C_LVD) * dP_RA # eq 44
            Ppa = Ppv + Q * Rp
            Ppa_also = P_thorax[j] + (C_RVD / C_LVD) * dP_RA + Q*Rp
            # assert (nabs(Ppa - Ppa_also) < 0.01)
            cases[j, i] = 1
        elif P_thorax[j] > - dP_RA and P_thorax[j] < rho * G[i] * Hu - dP_RA:#case 2
            a = Tp * (Gs_l_normal/Vd_total_normal_vec[j]**2)*(rho*G[i]*Hu + Psa_u_star - \
                (dP_RA + P_thorax[j]))
            b = 1
            c = -Vtotal + Cp*(C_RVD/C_LVD)*dP_RA + Tp * Gs_u * Psa_u_star + Csa*(Psa_u_star) + \
                Csa_l*rho*G[i]*Hu + Cs_l * rho *G[i]*(-Hl) + Csv_l * (dP_RA + P_thorax[j])
            Vd_total = find_zeros(a,b,c)
            Gs_l = Gs_l_normal * (Vd_total**2/Vd_total_normal_vec[j]**2)
            Rs_l = 1 / Gs_l
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
            F_also = ((Gs_u + Gs_l) * Psa_u_star +
                      (1/Rs_l)* (rho*G[i]*Hu-P_thorax[j] - dP_RA)) / (
                      C_RVD*dP_RA)

            Ppv = P_thorax[j] + (C_RVD / C_LVD) * dP_RA # eq 66
            Ppa = Ppv + Q * Rp
            Ppa_also = (C_RVD / C_LVD) * dP_RA + Rp*(Gs*Psa_u_star + 1/Rs_l *
                        (rho * G[i] *Hu - P_thorax[j] - dP_RA))
            cases[j, i] = 2
        elif P_thorax[j] >= rho * G[i] * Hu - dP_RA:# print("entered case III")
            
            a = (Gs_l_normal / Vd_total_normal_vec[j]**2) * Tp * (Psa_u_star + rho * G[i]* Hu - P_RA[j])
            b = 1
            c = -Vtotal + Cp*(C_RVD/C_LVD)*dP_RA + Tp * Gs_u * Psa_u_star + Csa*(Psa_u_star) + \
                Tp * Gs_u * rho * G[i] * Hu + Csa_l * rho * G[i] * Hu - Csv_u * rho * G[i] * \
                    (Hu) + Cs_l * rho * G[i] *(-Hl) + (Csv_l + Csv_u)* P_RA[j] - Tp * Gs_u * P_RA[j]
            Vd_total = find_zeros(a,b,c)
            Gs_l = mu * Gs_l_normal * (Vd_total**2/Vd_total_normal_vec[j]**2)
            Rs_l = 1 / Gs_l
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
            Gs_l_vec[i] = Gs_l
        else:
            Vd_total_vec[i] = float('nan')
            Q_vec[i] = float('nan')
            F_vec[i] = float('nan')
            Ppa_vec[i] = float('nan')
            Gs_l_vec[i] = float('nan')
    sol_Vd_Pthorax_G[j,:] = Vd_total_vec
    sol_Q_Pthorax_G[j,:] = Q_vec
    sol_F_Pthorax_G[j,:] = F_vec
    sol_Ppa_Pthorax_G[j,:] = Ppa_vec
    sol_Gsl_Pthorax_G[j,:] = Gs_l_vec
#conversions:
G = G / 100 / (g_earth / 100)
sol_Q_Pthorax_G = sol_Q_Pthorax_G * 60 / 1000
sol_F_Pthorax_G = sol_F_Pthorax_G * 60
sol_Ppa_Pthorax_G = sol_Ppa_Pthorax_G / 1333
sol_Vd_Pthorax_G = sol_Vd_Pthorax_G / 1000

#finding indices for each case:
line_style = ['--', ':', '-']
P_thorax_string = r'$P_\mathrm{thorax}$'
mmHg = r'$\mathrm{mmHg}$'
P_thorax_titles = [fr'{P_thorax_string} = {pt/1333} {mmHg}' for pt in P_thorax]
# Plotting G against Vd for different P_thorax values
fig, ax = plt.subplots()

for i in range(len(P_thorax)):
    #finding indices for each case:
    case_1_indices = np.where(cases[i,:] == 1)
    case_2_indices = np.where(cases[i,:] == 2)
    case_3_indices = np.where(cases[i,:] == 3)
    ax.plot(G[case_1_indices], sol_Vd_Pthorax_G[i, case_1_indices[0]],  color=line_colors[len(P_thorax) - i], linestyle = line_style[0])
    ax.plot(G[case_2_indices], sol_Vd_Pthorax_G[i, case_2_indices[0]],  color=line_colors[len(P_thorax) - i], linestyle = line_style[1])
    ax.plot(G[case_3_indices], sol_Vd_Pthorax_G[i, case_3_indices[0]], color=line_colors[len(P_thorax) - i], linestyle = line_style[2], label = fr'{P_thorax_string} = {P_thorax[i]/1333} {mmHg}')

color_legend = ax.legend(ncol=2, loc = 'upper right')
legend_handles = [
    Line2D([], [], linestyle='dashed', color='maroon', label=r'$\mathrm{Case}$ $1$'),
    Line2D([], [], linestyle='dotted', color='maroon', label=r'$\mathrm{Case}$ $2$'),
    Line2D([], [], linestyle='solid', color='maroon', label=r'$\mathrm{Case}$ $3$')
]
linestyle_legend = ax.legend(handles = legend_handles, loc = 'center right')
ax.add_artist(color_legend)
ax.set_xlim(1,4)

ax.set_xlabel(r'$+\mathrm{Gz}$ $(g$ $\mathrm{multiple})$')
ax.set_ylabel(r"$V_{\mathrm{total}}^0$ $(\mathrm{L}$)")
ax.set_title(r'$\mathrm{Reserve}$ $\mathrm{Volume}$ $\mathrm{v.}$ $\mathrm{Acceleration}$ $\mathrm{Varying}$ $\mathrm{Intrathoracic}$ $\mathrm{Pressures}$')
ax.tick_params(axis='both', labelsize=10)  # Set tick label font size
plt.rc('font', family='serif')  # Set font family to use LaTeX font
plt.grid(True)
plt.savefig('steady_state/figures/lowerGs_controller/varyPthorax_Vd_G', bbox_inches='tight', dpi=300)

# Plotting G against F for different P_thorax values
fig, ax = plt.subplots()
for i in range(len(P_thorax)):
    case_1_indices = np.where(cases[i,:] == 1)
    case_2_indices = np.where(cases[i,:] == 2)
    case_3_indices = np.where(cases[i,:] == 3)
    ax.plot(G[case_1_indices], sol_F_Pthorax_G[i, case_1_indices[0]], color=line_colors[len(P_thorax) - i], linestyle = line_style[0])
    ax.plot(G[case_2_indices], sol_F_Pthorax_G[i, case_2_indices[0]], color=line_colors[len(P_thorax) - i], linestyle = line_style[1])
    ax.plot(G[case_3_indices], sol_F_Pthorax_G[i, case_3_indices[0]], label=fr'{P_thorax_string} = {P_thorax[i]/1333} {mmHg}', color=line_colors[len(P_thorax) - i], linestyle = line_style[2])
color_legend = ax.legend(ncols = 2, loc = 'upper right')
linestyle_legend = ax.legend(handles = legend_handles, loc = 'center right')
ax.add_artist(color_legend)
ax.set_xlabel(r'$+\mathrm{Gz}$ $(g$ $\mathrm{multiple})$')
ax.set_ylabel(r'$F$ $\mathrm{(beats/min)}$')
ax.set_title(r'$\mathrm{Heart}$ $\mathrm{Rate}$ $\mathrm{v.}$ $\mathrm{Acceleration}$ $\mathrm{Varying}$ $\mathrm{Intrathoracic}$ $\mathrm{Pressures}$')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.grid(True)
plt.savefig('steady_state/figures/lowerGs_controller/varyPthorax_F_G', bbox_inches='tight', dpi=300)

# Plotting G against GsL for different P_thorax values
fig, ax = plt.subplots()
for i in range(len(P_thorax)):
    case_1_indices = np.where(cases[i,:] == 1)
    case_2_indices = np.where(cases[i,:] == 2)
    case_3_indices = np.where(cases[i,:] == 3)
    ax.plot(G[case_1_indices], sol_Gsl_Pthorax_G[i, case_1_indices[0]], color=line_colors[len(P_thorax) - i], linestyle = line_style[0])
    ax.plot(G[case_2_indices], sol_Gsl_Pthorax_G[i, case_2_indices[0]], color=line_colors[len(P_thorax) - i], linestyle = line_style[1])
    ax.plot(G[case_3_indices], sol_Gsl_Pthorax_G[i, case_3_indices[0]], label=fr'{P_thorax_string} = {P_thorax[i]/1333} {mmHg}', color=line_colors[len(P_thorax) - i], linestyle = line_style[2])
color_legend = ax.legend(ncols = 2, loc = 'upper right')
linestyle_legend = ax.legend(handles = legend_handles, loc = 'center right')
ax.add_artist(color_legend)
ax.set_xlabel(r'$+\mathrm{Gz}$ $(g$ $\mathrm{multiple})$')
ax.set_ylabel(r'$G_s^{\mathrm{l}}$ $\mathrm{((dynes/cm2)/(cm3/s))}$')
ax.set_title(r'$\mathrm{Lower}$ $\mathrm{Systemic}$ $\mathrm{Conductance}$ $\mathrm{v.}$ $\mathrm{Acceleration}$ $\mathrm{Varying}$ $\mathrm{Intrathoracic}$ $\mathrm{Pressures}$')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.grid(True)
plt.savefig('steady_state/figures/lowerGs_controller/varyPthorax_Gsl_G', bbox_inches='tight', dpi=300)




