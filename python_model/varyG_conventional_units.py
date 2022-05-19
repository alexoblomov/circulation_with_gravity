import numpy as np
import matplotlib.pyplot as plt

Psa_u_star = 100 * 1333

Psa_u = Psa_u_star
dP_RA = 2 * 1333
height = 167.64

Hu = 0.5 * 32
Hl = - (42 / 2)
Hr = -Hu/Hl #height_ratio
Hu_factor = 1-Hr
Hl_factor = Hr
rho = 1

g_earth = 980

G = np.linspace(980,980 * 10,1000)

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
P_thorax = np.linspace(- 4 * 1333,3 * 1333,8)

#P_thorax = -4*1333;
P_RA = P_thorax + dP_RA
Vd_total_vec = np.zeros(len(G))
Q_vec = np.zeros(len(G))
F_vec = np.zeros(len(G))
Ppa_vec = np.zeros(len(G))

sol_Vd_Pthorax_G = np.zeros((len(P_thorax),len(G)))
sol_Q_Pthorax_G = np.zeros((len(P_thorax),len(G)))
sol_F_Pthorax_G = np.zeros((len(P_thorax),len(G)))
sol_Ppa_Pthorax_G = np.zeros((len(P_thorax),len(G)))
for j in range(len(P_thorax)):
    for i in range(len(G)):
        if P_thorax[j] <= - dP_RA:
            Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * (dP_RA) - (Tp * Gs + Csa) * Psa_u_star - (Tp * Gs_l + Csa_l) * rho * G[i] * Hu - Cs_l * rho * G[i] * (- Hl)
            Q = ((1 / Rs_u) + (1 / Rs_l)) * Psa_u_star + rho * G[i] * Hu / Rs_l
            F = Q / (C_RVD * (P_RA[j] - P_thorax[j]))
            Ppv = P_thorax[j] + (C_LVD / C_RVD) * dP_RA
            Ppa = Ppv + Q * Rp
        else:
            if P_thorax[j] > - dP_RA and P_thorax[j] < rho * G[i] * Hu - dP_RA:
                Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA - (Tp * Gs + Csa) * Psa_u_star - (Tp * Gs_l + Csa_l) * rho * G[i] * Hu - Cs_l * rho * G[i] * (- Hl) - (Csv_l - Tp * Gs_l) * (P_thorax[j] + dP_RA)
                Psv_l = - rho * G[i] * Hl + P_thorax[j] + dP_RA
                Psv_u = 0
                Psa_l = Psa_u_star + rho * G[i] * (Hu - Hl)
                Qs_u = Psa_u_star / Rs_u
                Qs_l = (Psa_l - Psv_l) / Rs_l
                Q = Qs_u + Qs_l
                F = Q / (C_RVD * (dP_RA))
                Ppv = P_thorax[j] + (C_LVD / C_RVD) * dP_RA
                Ppa = Ppv + Q * Rp
            else:
                if P_thorax[j] >= rho * G[i] * Hu - dP_RA:
                    Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA - (Tp * Gs + Csa) * Psa_u_star - (Tp * Gs + Csa_l - Csv_u) * rho * G[i] * Hu - Cs_l * rho * G[i] * (- Hl) - (Csv_l - Tp * Gs) * (P_thorax[j] + dP_RA)
                    Psv_l = P_thorax[j] + dP_RA + rho * G[i] * (- Hl)
                    Psv_u = P_thorax[j] + dP_RA - rho * G[i] * Hu
                    Psa_l = Psa_u_star + rho * G[i] * (Hu - Hl)
                    Qs_u = (Psa_u_star - Psv_u) / Rs_u
                    Qs_l = (Psa_l - Psv_l) / Rs_l
                    Q = Qs_u + Qs_l
                    F = Q / (C_RVD * (dP_RA))
                    Ppv = P_thorax[j] + (C_LVD / C_RVD) * dP_RA
                    Ppa = Ppv + Q * Rp
        if Vd_total > 0:
            Vd_total_vec[i] = Vd_total
            Q_vec[i] = Q
            F_vec[i] = F
            Ppa_vec[i] = Ppa
        else:
            Vd_total_vec[i] = np.nan
            Q_vec[i] = np.nan
            F_vec[i] = np.nan
            Ppa_vec[i] = np.nan
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
myLineColorPref = 'k-'

myLineColorVec = np.array([0,0,0])

myLineWidth = 1
myLabelFontSize = 18
alpha = 0.1

beta = 1
#label_P_thorax = num2str([-3:4]'); #in mm Hg, from -3 to 4

## saveas(h_overlay_nonDM, myPlotOverlay_nonDM, 'pdf')
#set(gca, 'FontSize', myLabelFontSize)
#print(myfig,figName,"-dpdf", "-S4016,2362")

#family of plots for Reserve Volume vs. G for different Pthorax values:

h = plt.figure(100)
plt.plot(G,sol_Vd_Pthorax_G[1,:])

# for i in np.arange(1,len(P_thorax)+1).reshape(-1):
#     myLineColor = myLineColorVec + alpha * (i - 1)
#     #myLineWidthPlot=myLineWidth+beta*(i-1); #set linewidth
#     myLineWidthPlot = myLineWidth
#     plt.plot(G,sol_Vd_Pthorax_G[i,:])

plt.xlabel('G')
plt.ylabel('Reserve Volume (L)')
# set(gca,'FontSize',myLabelFontSize)
# saveas(h,'QvsV0','pdf')
# saveas(h,'QvsV0','png')
plt.show()
# xticks([1*9.80 2*9.80 3*9.80 4*9.80 5*9.80 6*9.80 7*9.80 8*9.80 9*9.80 10*9.80])
# xticklabels({'g','2g','3g','4g','5g','6g','7g', '8g', '9g', '10g'})

# h1 = plt.figure(101)
# clf(101)
# for i in np.arange(1,len(P_thorax)+1).reshape(-1):
#     plt.plot(G,sol_Q_Pthorax_G(i,:))
#     plt.xlabel('G','interpreter','latex')
#     plt.ylabel('Cardiac Output (L/min)','interpreter','latex')
#     hold('on')

# xticks([1*9.80 2*9.80 3*9.80 4*9.80 5*9.80 6*9.80 7*9.80 8*9.80 9*9.80 10*9.80])
# xticklabels({'g','2g','3g','4g','5g','6g','7g', '8g', '9g', '10g'})

# h2 = plt.figure(102)
# clf(102)
# hold('on')
# for i in np.arange(1,len(P_thorax)+1).reshape(-1):
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
# for i in np.arange(1,len(P_thorax)+1).reshape(-1):
#     plt.plot(G,sol_Ppa_Pthorax_G(i,:))
#     plt.xlabel('G','interpreter','latex')
#     plt.ylabel('Pulmonary Arterial Pressure (mmHg)','interpreter','latex')
#     hold('on')

# xticks([1*9.80 2*9.80 3*9.80 4*9.80 5*9.80 6*9.80 7*9.80 8*9.80 9*9.80 10*9.80])
# xticklabels({'g','2g','3g','4g','5g','6g','7g', '8g', '9g', '10g'})
