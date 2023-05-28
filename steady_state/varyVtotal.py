"""
Purpose: simulate model with varying values of G and plot output variables
of interest
Nb case III is only actualized when P_thorax is roughly 70x its nominal value
"""
import numpy as np
import matplotlib.pyplot as plt
from parameters import *

G = np.linspace(g_earth,7*980, 3000)

Vtotal = np.linspace(3, 10, 8)
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

plt.figure(figsize=(15, 12))
plt.subplots_adjust(hspace=0.5)
plt.suptitle("Reserve volume vs g", fontsize=18, y=0.95)
Vtotal_titles = ["Total volume " + str(vtot/1000) + " L" for vtot in Vtotal]
# loop through the length of tickers and keep track of index
for n, plt_title in enumerate(Vtotal_titles):
    # add a new subplot iteratively
    ax = plt.subplot(4, 2, n + 1)

    # filter df and plot ticker on the new subplot axis
    idx_case_1 = cases == 1
    idx_case_2 = cases == 2
    idx_case_3 = cases == 3
    ax.plot(G,sol_Vd_Vtotal_G[n,:])
    # ax.plot(G,sol_Vd_Pthorax_G[n,idx_case_2],'g')
    # ax.plot(G,sol_Vd_Pthorax_G[n,idx_case_3],'b')
    # chart formatting
    ax.set_title(plt_title)
    # ax.get_legend().remove()
    ax.set_xlabel("g multiple")
    ax.set_ylabel("VT0")
# plt.show()
plt.savefig("varyVtot_VT0_vs_g_V0_w_height_factor.png")

plt.figure()
for n, plt_title in enumerate(Vtotal_titles):
    # add a new subplot iteratively
    # filter df and plot ticker on the new subplot axis
    idx_case_1 = cases == 1
    idx_case_2 = cases == 2
    idx_case_3 = cases == 3
    plt.plot(G, sol_Vd_Vtotal_G[n, :])
    plt.title('Reserve Volume v. Acceleration')
    # ax.get_legend().remove()
    plt.xlabel("g multiple")
    plt.ylabel("VT0")
plt.legend(Vtotal_titles)
plt.grid(True)
plt.show()
plt.savefig("varyVT_VT0_vs_g_V0")

plt.figure(figsize=(12, 8))
plt.title("G Tolerance vs total volume", fontsize=18)
plt.xlabel("Vtotal (L)", fontsize=14)
plt.ylabel("G Tolerance ", fontsize=14)

for i in range(len(Vtotal)):
    min_Vd_total = np.nanmin(sol_Vd_Vtotal_G[i, :])
    if np.isnan(min_Vd_total):
        g_intercept = G[np.isnan(sol_Vd_Vtotal_G[i, :])][-1]
    else:
        min_Vd_total_index = np.nanargmin(sol_Vd_Vtotal_G[i, :])
        g_intercept = G[min_Vd_total_index]
    plt.plot(Vtotal[i] / 1000, g_intercept, 'ro', markersize=6)

plt.grid(True)
plt.show()
