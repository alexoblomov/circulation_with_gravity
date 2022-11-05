"""
This parameter study of varying heights is conducted for case I only
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_theme()
from matplotlib.ticker import FormatStrFormatter

from parameters import *

# Pressures:
P_thorax = -4*1333
P_RA = P_thorax + dP_RA

# for a given hypothetical maximum HR
max_F = 180 # maximum allowed heart rate

# WE NEED TO SOLVE WHOLE SYSTEM THEN CALCULATE INTERCEPT. PLACEHOLDER CODE
# for a given g
dx = 100
Hu_range = np.linspace(0.5*Hu_patient, 1.8*Hu_patient, dx)
Hl_range = np.linspace(0.5*Hl_patient, 1.5*Hl_patient, dx)

gtol_vs_Hu_vs_Hl= np.zeros((len(Hu_range), len(Hl_range)))
F = max_F
for i in range(dx):
    for j in range(dx):
            # Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * (dP_RA) - ...
            #           (Tp * Gs + Csa) * Psa_u_star - (Tp * Gs_l + Csa_l) ...
            #           * rho * g_earth * Hu - Cs_l * rho * g_earth * (- Hl)
            # Setting Vd_total to zero to find g tol

            Hu = Hu_range[i]
            Hl = Hl_range[j]
            g_tol = (Vtotal - Cp*(C_RVD/C_LVD)*(dP_RA) -
                                (Tp*Gs+Csa)*Psa_u_star)/((Tp*Gs_l+ Csa_l)
                                                         *rho*Hu-Cs_l*rho*Hl)
            gtol_vs_Hu_vs_Hl[i,j] = g_tol
            Q = ((1 / Rs_u) + (1 / Rs_l)) * Psa_u_star + rho * g_tol * Hu / Rs_l
            Ppv = P_thorax + (C_LVD / C_RVD) * dP_RA
            Ppa = Ppv + Q * Rp
        # cases[j, i] = 1

# conversions:
gtol_vs_Hu_vs_Hl = gtol_vs_Hu_vs_Hl / 100 / (g_earth / 100)
Hu_range = np.around(Hu_range, 2)
Hl_range = np.around(Hu_range, 2)
# heatmap plot

num_ticks = 10
# the index of the position of yticks
yticks_ = np.linspace(0, dx - 1, num_ticks, dtype=int)
# the content of labels of these yticks
yticklabels_ = [str(Hu_range[idx]) for idx in yticks_]

xticks_ = np.linspace(0, dx - 1, num_ticks, dtype=int)
# the content of labels of these yticks
xticklabels_ = [str(Hl_range[idx]) for idx in xticks_]

fig, ax = plt.subplots()
im = ax.imshow(gtol_vs_Hu_vs_Hl)
fig.colorbar(im, ax=ax)
# ax = sns.heatmap(gtol_vs_CsaU_vs_CsaL)

# Show all ticks and label them with the respective list entries
ax.set_xticks(xticks_, labels=xticklabels_)
ax.set_yticks(yticks_, labels=yticklabels_)

# optional if you want to put ticks on the top
# ax.xaxis.tick_top()
# ax.xaxis.set_label_position('top')
# ax.tick_params(axis='x', which='major', pad=18, direction='out', rotation=45)  # move the tick labels

# plt.setp( ax.xaxis.get_majorticklabels(), rotation=45, ha="left", rotation_mode="anchor")
# plt.setp( ax.xaxis.get_majorticklabels(), rotation=45, ha="right" )

ax.grid(False)
# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
fig.tight_layout(pad=3)

plt.ylabel('heart to head height', fontsize = 12)
plt.xlabel('heart to seat height', fontsize = 12)
plt.savefig('heights_gtol_height_factor.png')
