"""
Purpose: simulate model with varying values of G and plot output variables
of interest
Nb case III is only actualized when P_thorax is roughly 70x its nominal value
"""


import numpy as np
import matplotlib.pyplot as plt
from parameters import *

G = np.linspace(g_earth,1.7*980, 3)

P_thorax = np.linspace(- 4 * 1333,3 * 1333,8)

#P_thorax = -4*1333;
P_RA = P_thorax + dP_RA

cases = np.zeros((len(P_thorax),len(G)), dtype=np.uint8)

Vd_total_vec = np.zeros(len(G))
Q_vec = np.zeros(len(G))
F_vec = np.zeros(len(G))
Ppa_vec = np.zeros(len(G))

sol_Vd_Pthorax_G = np.zeros((len(P_thorax),len(G)))
sol_Q_Pthorax_G = np.zeros((len(P_thorax),len(G)))
sol_F_Pthorax_G = np.zeros((len(P_thorax),len(G)))
sol_Ppa_Pthorax_G = np.zeros((len(P_thorax),len(G)))

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
            Ppa_also = P_thorax[i] + (C_RVD / C_LVD) * dP_RA + Q*Rp
            # assert (nabs(Ppa - Ppa_also) < 0.01)
            cases[j, i] = 1
        elif P_thorax[j] > - dP_RA and P_thorax[j] < rho * G[i] * Hu - dP_RA:
            Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA \
                      - (Tp * Gs + Csa)* Psa_u_star \
                      - (Tp * Gs_l + Csa_l) * (rho * G[i] * Hu - Cs_l *
                                               rho * G[i] * (- Hl)) \
                      - (Csv_l - Tp * Gs_l) * (P_thorax[j] + dP_RA)
            Psv_l = - rho * G[i] * Hl + P_thorax[j] + dP_RA
            Psv_u = 0
            Psa_l = Psa_u_star + rho * G[i] * (Hu - Hl)
            # Qs_u = Psa_u_star / Rs_u
            # Qs_l = (Psa_l - Psv_l) / Rs_l
            Qs_u = Psa_u_star * 1 / Rs_u # eq 61
            #eq 62
            Qs_l = (Psa_u_star + rho * G[i] * Hu \
                    - P_thorax[j]-dP_RA) * 1/Rs_l
            Q = Qs_u + Qs_l
            F = Q / (C_RVD * (dP_RA))
            F_also = (Gs * Psa_u_star +
                      (1/Rs_l)* (rho*G[i]*Hu-P_thorax[i] - dP_RA)) / (
                      C_RVD*dP_RA)
            # breakpoint()
            # assert F == F_also, "F, equation 64 inconsistent"
            # (Pdb) F
            # 1.7248725549967021
            # (Pdb) F_also
            # 1.7494574795209836

            Ppv = P_thorax[j] + (C_RVD / C_LVD) * dP_RA # eq 66
            Ppa = Ppv + Q * Rp
            Ppa_also = (C_RVD / C_LVD) * dP_RA + Rp*(Gs*Psa_u_star + 1/Rs_l *
                        (rho * G[i] *Hu - P_thorax[i] - dP_RA))
            # (Pdb) Ppa
            # 30219.708890709935
            # (Pdb) Ppa_also
            # 29308.135553880078
            # assert Ppa == Ppa_also, "Ppa eq 67 inconsistent"
            cases[j, i] = 2
        elif P_thorax[j] >= rho * G[i] * Hu - dP_RA:
            # print("entered case III")
            Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA - \
                    (Tp*Gs + Csa_l+Csa_u)*Psa_u_star \
                    +(Tp*Gs + Csa_u - Csv_l) * rho * G[i] * Hu \
                    + (Csa_l+Csv_l) * rho * G[i] * (-Hl) \
                    + (Csv_l+Csv_u - Tp*Gs)* (P_thorax[j] - dP_RA)

            Psv_l = P_thorax[j] + dP_RA + rho * G[i] * (- Hl)
            Psv_u = P_thorax[j] + dP_RA - rho * G[i] * Hu
            Psa_l = Psa_u_star + rho * G[i] * (Hu - Hl)
            Qs_u = (Psa_u_star - Psv_u) / Rs_u
            Qs_u_also = (Psa_u_star + rho*G[i]*Hu - P_thorax[j] - dP_RA)/Rs_u
            # assert  Qs_u == Qs_u_also, "QsU eq 81 inconsistent"
            # (Pdb) Qs_u
            # 40.523557770364924
            # (Pdb) Qs_u_also
            # 40.523557770364924

            # breakpoint()
            Qs_l = (Psa_l - Psv_l) / Rs_l
            Qs_l_also = (Psa_u_star + rho*G[i]*Hu - P_thorax[j] - dP_RA)/Rs_l
            # Qs_l == Qs_l_also, "QsL eq 82 inconsistent"
            # (Pdb) Qs_l
            # 53.18716957360396
            # (Pdb) Qs_l_also
            # 53.18716957360396

            Q = Qs_u + Qs_l
            F = Q / (C_RVD * (dP_RA))
            Ppv = P_thorax[j] + (C_RVD / C_LVD) * dP_RA
            Ppa = Ppv + Q * Rp
            cases[j, i] = 3
        if Vd_total > 0:
            Vd_total_vec[i] = Vd_total
            Q_vec[i] = Q
            F_vec[i] = F
            Ppa_vec[i] = Ppa
        else:
            Vd_total_vec[i] = nnan
            Q_vec[i] = nnan
            F_vec[i] = nnan
            Ppa_vec[i] = nnan
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
### PLOT OPTIONS ###

##plot(x, y,myLineColorPref,'color', myLineColorVec,'LineWidth', myLineWidth) #buffer times in black
# myLineColorPref = 'k-'
#
# myLineColorVec = narray([0,0,0])
#
# myLineWidth = 1
# myLabelFontSize = 18
# alpha = 0.1
#
# beta = 1

#family of plots for Reserve Volume vs. G for different Pthorax values:

# plt.show()
# xticks([1*9.80 2*9.80 3*9.80 4*9.80 5*9.80 6*9.80 7*9.80 8*9.80 9*9.80 10*9.80])
# xticklabels({'g','2g','3g','4g','5g','6g','7g', '8g', '9g', '10g'})


plt.figure(figsize=(15, 12))
plt.subplots_adjust(hspace=0.5)
plt.suptitle("Reserve volume vs g", fontsize=18, y=0.95)
pthorax_titles = ["P thorax " + str(pth / 1333) + " mmHg" for pth in P_thorax]
# loop through the length of tickers and keep track of index
for n, plt_title in enumerate(pthorax_titles):
    # add a new subplot iteratively
    ax = plt.subplot(4, 2, n + 1)

    # filter df and plot ticker on the new subplot axis
    idx_case_1 = cases == 1
    idx_case_2 = cases == 2
    idx_case_3 = cases == 3
    ax.plot(G,sol_Vd_Pthorax_G[n,:])
    # ax.plot(G,sol_Vd_Pthorax_G[n,idx_case_2],'g')
    # ax.plot(G,sol_Vd_Pthorax_G[n,idx_case_3],'b')
    # chart formatting
    ax.set_title(plt_title)
    # ax.get_legend().remove()
    ax.set_xlabel("g multiple")
    ax.set_ylabel("VT0")
# plt.show()
plt.savefig("fixed_VT0_vs_g_V0_w_height_factor.png")

#################### code needs to be adapted from here on #####################
# h1 = plt.figure(101)
# clf(101)
# for i in narange(1,len(P_thorax)+1).reshape(-1):
#     plt.plot(G,sol_Q_Pthorax_G(i,:))
#     plt.xlabel('G','interpreter','latex')
#     plt.ylabel('Cardiac OuTput (L/min)','interpreter','latex')
#     hold('on')

# xticks([1*9.80 2*9.80 3*9.80 4*9.80 5*9.80 6*9.80 7*9.80 8*9.80 9*9.80 10*9.80])
# xticklabels({'g','2g','3g','4g','5g','6g','7g', '8g', '9g', '10g'})

# h2 = plt.figure(102)
# clf(102)
# hold('on')
# for i in narange(1,len(P_thorax)+1).reshape(-1):
#     myLineColor = myLineColorVec + alpha * (i - 1)
#     plt.plot(G,sol_F_Pthorax_G(i,:),myLineColorPref,'color',myLineColor,'LineWidth',myLineWidthPlot)
#
# hold('off')
# plt.xlabel('G')
# plt.ylabel('Heart Rate (per minute)')
# hold('on')
# xticks([1*9.80 2*9.80 3*9.80 4*9.80 5*9.80 6*9.80 7*9.80 8*9.80 9*9.80 10*9.80])
# xticklabels({'g','2g','3g','4g','5g','6g','7g', '8g', '9g', '10g'})
#
# h3 = plt.figure(103)
# clf(103)
# for i in narange(1,len(P_thorax)+1).reshape(-1):
#     plt.plot(G,sol_Ppa_Pthorax_G(i,:))
#     plt.xlabel('G','interpreter','latex')
#     plt.ylabel('Pulmonary Arterial Pressure (mmHg)','interpreter','latex')
#     hold('on')

# xticks([1*9.80 2*9.80 3*9.80 4*9.80 5*9.80 6*9.80 7*9.80 8*9.80 9*9.80 10*9.80])
# xticklabels({'g','2g','3g','4g','5g','6g','7g', '8g', '9g', '10g'})
