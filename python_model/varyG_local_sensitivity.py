import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_theme()

Psa_u_star = 100 * 1333
max_F = 180 # maximum allowed heart rate
Psa_u = Psa_u_star
dP_RA = 2 * 1333
height = 167.64

Hu_karen = 32
Hl_karen = -42
# Hr = -Hu/Hl #height_ratio
# Hu_factor = 1-Hr
# Hl_factor = Hr
rho = 1

g_earth = 980

G = np.linspace(980,980 * 10,100)

Rs = (16.49) * 1333 / (1000 / 60)

Gs = 1 / Rs
Gs_u = (1 / 3) * Gs
Gs_l = (2 / 3) * Gs
Rs_l = 1 / Gs_l
Rs_u = 1 / Gs_u
Rp = (1.61 * 1333) / (1000 / 60)

C_RVD = (0.035 / 1333) * 1000

C_LVD = (0.00583 / 1333) * 1000

Csa_l = (2 / 3) * (0.00175 / 1333) * 1000
Csa_u = (1 / 3) * (0.00175 / 1333) * 1000
Csv_l = (2 / 3) * (0.09 / 1333) * 1000
Csv_u = (1 / 3) * (0.09 / 1333) * 1000
Cs_l = Csa_l + Csv_l
Cpa = (0.00412 / 1333) * 1000
Cpv = (0.01 / 1333) * 1000
Cp = Cpa + Cpv
Vtotal = 3.7 * 1000

Csa = Csa_l + Csa_u
Gs = 1 / Rs_u + 1 / Rs_l
Gs_l = 1 / Rs_l
Ts = Csa_u * Rs_u

Tp = Rp * Cpa
Csa = Csa_u + Csa_l
# P_thorax = np.linspace(- 4 * 1333,3 * 1333,8)
# P_thorax = np.linspace(- 4 * 1333, 3, 1)
P_thorax = -4*1333

#P_thorax = -4*1333;
P_RA = P_thorax + dP_RA

# for a given g
dx = 10
Hu_range = np.linspace(0.5*Hu_karen, 1.8*Hu_karen, dx)
Hl_range = np.linspace(0.5*Hl_karen, 1.5*Hl_karen, dx)

gtol_vs_Hu_vs_Hl= np.zeros((len(Hu_range), len(Hl_range)))
# WE NEED TO SOLVE WHOLE SYSTEM THEN CALCULATE INTERCEPT. PLACEHOLDER CODE
for g in range(len(G)):
    for i in range(dx):
        for j in range(dx):
            if P_thorax <= - dP_RA:
                Q = ((1 / Rs_u) + (1 / Rs_l)) * Psa_u_star + rho * G[g] * Hu / Rs_l
                F = Q / (C_RVD * (P_RA - P_thorax))
                if F < max_F:
                    Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * (dP_RA) - ...
                              (Tp * Gs + Csa) * Psa_u_star - (Tp * Gs_l + Csa_l) ...
                              * rho * G[i] * Hu - Cs_l * rho * G[g] * (- Hl)
                    # fix before continuing
                    break
                    # Setting Vd_total to zero to find g tol
                    Hu = Hu_karen[i]
                    Hl = Hl_karen[j]
                    gtol_vs_Hu_vs_Hl[i,j] = (Vtotal - Cp*(C_RVD/C_LVD)*(dP_RA) -
                                        (Tp*Gs+Csa)*Psa_u_star)/((Tp*Gs_l+ Csa_l)
                                                                 *rho*Hu-Cs_l*rho*Hl)
                # Q = ((1 / Rs_u) + (1 / Rs_l)) * Psa_u_star + rho * G[i] * Hu / Rs_l

                # Ppv = P_thorax[j] + (C_LVD / C_RVD) * dP_RA
                # Ppa = Ppv + Q * Rp
                # cases[j, i] = 1
            # elif P_thorax[j] > - dP_RA and P_thorax[j] < rho * G[i] * Hu - dP_RA:
            #         Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA - (Tp * Gs + Csa) * Psa_u_star - (Tp * Gs_l + Csa_l) * rho * G[i] * Hu - Cs_l * rho * G[i] * (- Hl) - (Csv_l - Tp * Gs_l) * (P_thorax[j] + dP_RA)
            #         Psv_l = - rho * G[i] * Hl + P_thorax[j] + dP_RA
            #         Psv_u = 0
            #         Psa_l = Psa_u_star + rho * G[i] * (Hu - Hl)
            #         Qs_u = Psa_u_star / Rs_u
            #         Qs_l = (Psa_l - Psv_l) / Rs_l
            #         Q = Qs_u + Qs_l
            #         F = Q / (C_RVD * (dP_RA))
            #         Ppv = P_thorax[j] + (C_LVD / C_RVD) * dP_RA
            #         Ppa = Ppv + Q * Rp
            #         cases[j, i] = 2
            # else:
            #     if P_thorax[j] >= rho * G[i] * Hu - dP_RA:
            #         Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA - (Tp * Gs + Csa) * Psa_u_star - (Tp * Gs + Csa_l - Csv_u) * rho * G[i] * Hu - Cs_l * rho * G[i] * (- Hl) - (Csv_l - Tp * Gs) * (P_thorax[j] + dP_RA)
            #         Psv_l = P_thorax[j] + dP_RA + rho * G[i] * (- Hl)
            #         Psv_u = P_thorax[j] + dP_RA - rho * G[i] * Hu
            #         Psa_l = Psa_u_star + rho * G[i] * (Hu - Hl)
            #         Qs_u = (Psa_u_star - Psv_u) / Rs_u
            #         Qs_l = (Psa_l - Psv_l) / Rs_l
            #         Q = Qs_u + Qs_l
            #         F = Q / (C_RVD * (dP_RA))
            #         Ppv = P_thorax[j] + (C_LVD / C_RVD) * dP_RA
            #         Ppa = Ppv + Q * Rp
            #         cases[j, i] = 3
        #     if Vd_total > 0:
        #         Vd_total_vec[i] = Vd_total
        #         Q_vec[i] = Q
        #         F_vec[i] = F
        #         Ppa_vec[i] = Ppa
        #     else:
        #         Vd_total_vec[i] = np.nan
        #         Q_vec[i] = np.nan
        #         F_vec[i] = np.nan
        #         Ppa_vec[i] = np.nan
        # sol_Vd_Pthorax_G[j,:] = Vd_total_vec
        # sol_Q_Pthorax_G[j,:] = Q_vec
        # sol_F_Pthorax_G[j,:] = F_vec
        # sol_Ppa_Pthorax_G[j,:] = Ppa_vec

#conversions:
# G = G / 100 / (g_earth / 100)

### PLOTS ###
ax = sns.heatmap(uniform_data)
