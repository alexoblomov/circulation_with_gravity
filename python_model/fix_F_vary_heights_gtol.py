import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_theme()
from matplotlib.ticker import FormatStrFormatter
Psa_u_star = 100 * 1333
Psa_u = Psa_u_star
height = 167.64

Hu_karen = 32
Hl_karen = -42
# Hr = -Hu/Hl #height_ratio
# Hu_factor = 1-Hr
# Hl_factor = Hr
rho = 1
g_earth = 980

# Resistance and compliance constants

# systemic resistance (mmHg/(liters/s))
Rs = (16.49) * 1333 / (1000 / 60)

Gs = 1 / Rs
Gs_u = (1 / 3) * Gs
Gs_l = (2 / 3) * Gs
Rs_l = 1 / Gs_l
Rs_u = 1 / Gs_u
Rp = (1.61 * 1333) / (1000 / 60)

# right-ventricular diastolic compliance (liters/mmHg)
C_RVD = (0.035 / 1333) * 1000
# left
C_LVD = (0.00583 / 1333) * 1000

# other compliances
Csa_l = (2 / 3) * (0.00175 / 1333) * 1000
Csa_u = (1 / 3) * (0.00175 / 1333) * 1000
Csv_l = (2 / 3) * (0.09 / 1333) * 1000
Csv_u = (1 / 3) * (0.09 / 1333) * 1000
Cpa = (0.00412 / 1333) * 1000
Cpv = (0.01 / 1333) * 1000
Cp = Cpa + Cpv
Vtotal = 3.7 * 1000

# aggregate constants
Csa = Csa_l + Csa_u
Cs_l = Csa_l + Csv_l

Gs = 1 / Rs_u + 1 / Rs_l
Gs_l = 1 / Rs_l

# time constants
Ts = Csa_u * Rs_u
Tp = Rp * Cpa

# P_thorax = np.linspace(- 4 * 1333,3 * 1333,8)
# P_thorax = np.linspace(- 4 * 1333, 3, 1)

# Pressures:
dP_RA = 2 * 1333
P_thorax = -4*1333
P_RA = P_thorax + dP_RA

# for a given hypothetical maximum HR
max_F = 180 # maximum allowed heart rate

# WE NEED TO SOLVE WHOLE SYSTEM THEN CALCULATE INTERCEPT. PLACEHOLDER CODE
# for a given g
dx = 10
Hu_range = np.linspace(0.5*Hu_karen, 1.8*Hu_karen, dx)
Hl_range = np.linspace(0.5*Hl_karen, 1.5*Hl_karen, dx)

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

# heatmap plot

num_ticks = 10
# the index of the position of yticks
yticks_ = np.linspace(0, dx - 1, num_ticks, dtype=np.int)
# the content of labels of these yticks
yticklabels_ = [str(Hu_range[idx]) for idx in yticks_]

xticks_ = np.linspace(0, dx - 1, num_ticks, dtype=np.int)
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

ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.grid(False)
# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
fig.tight_layout(pad=3)

plt.ylabel('heart to head height', fontsize = 12)
plt.xlabel('heart to seat height', fontsize = 12)
plt.savefig('heights_gtol.png')
