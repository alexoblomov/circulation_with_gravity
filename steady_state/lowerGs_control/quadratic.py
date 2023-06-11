import cmath
from parameters import * 
import numpy as np

def find_zeros(a, b, c):
    roots = np.roots([a, b, c])  # Find all roots

    # Filter out the real roots
    real_roots = roots[np.isreal(roots)].real

    # Select the larger root if there are multiple real roots
    if len(real_roots) > 1:
        larger_root = max(real_roots)
        return larger_root
    elif len(real_roots) == 1:
        return real_roots[0]  # Return the only real root
    else:
        return float('nan')  # No real roots found


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc
from parameters import *
from quadratic import find_zeros
mpl.rcParams['mathtext.fontset'] = 'cm'

# Enable LaTeX rendering
mpl.rcParams['text.usetex'] = True

# Set the font family to serif
mpl.rcParams['font.family'] = 'serif'



G = np.linspace(g_earth,3*980, 300)

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
    for i in range(len(G)):
        if P_thorax[j] <= - dP_RA:
            a = (Tp * Gs_l_normal / Vd_total_normal_vec[j]**2) * (Psa_u_star + rho * G[i] * Hu)
            b = 1
            c = -Vtotal + Cp*(C_RVD/C_LVD)*dP_RA + Tp * Gs_u * Psa_u_star + Csa*(Psa_u_star) \
                + Csa_l*rho*G[i]*Hu + Cs_l * rho *G[i]*(-Hl)
            Vd_total = find_zeros(a,b,c)
            cases[j, i] = 1
        elif P_thorax[j] > - dP_RA and P_thorax[j] < rho * G[i] * Hu - dP_RA:
            a = Tp * (Gs_l_normal/Vd_total_normal_vec[j]**2)*(rho*G[i]*Hu + Psa_u_star - \
                (dP_RA + P_thorax[j]))
            b = 1
            c = -Vtotal + Cp*(C_RVD/C_LVD)*dP_RA + Tp * Gs_u * Psa_u_star + Csa*(Psa_u_star) + \
                Csa_l*rho*G[i]*Hu + Cs_l * rho *G[i]*(-Hl) + Csv_l * (dP_RA + P_thorax[j])
            Vd_total = find_zeros(a,b,c)

            cases[j, i] = 2
        elif P_thorax[j] >= rho * G[i] * Hu - dP_RA:
            # print("entered case III")
            a = (Gs_l_normal / Vd_total_normal_vec[j]**2) * Tp * (Psa_u_star + rho * G[i]* Hu - P_RA[j])
            b = 1
            c = -Vtotal + Cp*(C_RVD/C_LVD)*dP_RA + Tp * Gs_u * Psa_u_star + Csa*(Psa_u_star) + \
                Tp * Gs_u * rho * G[i] * Hu + Csa_l * rho * G[i] * Hu - Csv_u * rho * G[i] * \
                    (Hu) + Cs_l * rho * G[i] *(-Hl) + (Csv_l + Csv_u)* P_RA[j] - Tp * Gs_u * P_RA[j]
            Vd_total = find_zeros(a,b,c)
            cases[j, i] = 3

        print(Vd_total)