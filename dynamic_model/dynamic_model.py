"""
Purpose : simulate ode model

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import odeint, solve_ivp
from parameters import *
from volume_odes import *
from controller import (get_linear_heart_rate, get_exp_heart_rate,
                        get_lower_peripheral_resistance)

control_type = "linear" # linear or exp -- which controller we call


# time steps
# forward euler time step
n_timesteps = 1000
n_seconds = 10
T = np.linspace(0.01, n_seconds, num=n_timesteps)
h = n_seconds/n_timesteps


# INPUTS 
# g_range = g_earth*np.ones(n_timesteps)
g_range = g_earth*np.zeros(n_timesteps)
Pext_l = Pext * np.ones(n_timesteps)
P_thorax = (- 4 * 1333) * np.ones(n_timesteps)


# TODO vary normally around the nominal values of the relative heights
Hu = Hu_patient
Hl = Hl_patient
H = Hu + Hl


# may want to make the following 2D in the future as we
# vary gravity and/or other parameters

# OUTPUTS
Q = np.zeros(n_timesteps)
P_ra = np.zeros(n_timesteps)
Psv_u = np.zeros(n_timesteps)
Psv_l = np.zeros(n_timesteps)
Psa_u = np.zeros(n_timesteps)
Psa_l = np.zeros(n_timesteps)
F = np.zeros(n_timesteps) # set by controller inside loop

# TODO consider moving to parameters.py
dP_RA = np.zeros(n_timesteps)
dP_RA[0] = init_dP_RA




# assume venous vol is roughly 70% of total BV.
# TODO: make this a normal random variable centered around 70%
pct_venous = 0.7
pct_arterial = 1 - pct_venous
init_Vsa = Vtotal * pct_arterial
init_Vsv = Vtotal * pct_venous
Vsa = np.zeros(n_timesteps+1)
Vsa[0] = init_Vsa
Vsv = np.zeros(n_timesteps+1)
Vsv[0] = init_Vsv

# reseve volumes
Vsa0 = np.ones(n_timesteps) * (1-pct_venous )* VT0_steady_state_cntrl
Vsv0 = np.ones(n_timesteps) * pct_venous * VT0_steady_state_cntrl


for t in range(n_timesteps):
    g = g_range[t]
    # eq 20 => six cases for Pra:
    # eq 20 gives Pra
    if Pext_l[t] < rho * g * Hu:
        if P_ra[t] <= Pext_l[t]:
            P_ra[t] = (Vsv[t] - Vsv0[t] -
                        Csv_l*(rho*g*Hl) +
                        Cra * P_thorax[t]) / Cra
        elif P_ra[t] >= Pext_l[t]:
            P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l *
                        (rho * g * Hl - Pext_l[t]) + Cra * P_thorax[t]
                        ) / (Cra + Csv_l)
        elif rho * g * Hu <= P_ra:
            P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l * (rho * g * Hl - Pext_l[t])
                        + Csv_u * rho * g * Hu + Cra *
                        P_thorax[t]) / (Cra + Csv_u + Csv_l)
        else:
            raise ValueError
    else:
        if P_ra[t] <= rho * g * Hu:
            P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l * (rho * g * Hl
                                                    - Pext_l[t])
                        + Cra * P_thorax[t]- Csv_l*Pext_l[t]) / Cra
        elif P_ra[t] >= rho * g * Hu:
            P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l * rho * g * Hl +
                        Csv_u * rho * g * Hu + Cra * P_thorax[t]
                        ) / (Cra + Csv_u)
        elif Pext_l[t] <= P_ra[t]:
            P_ra[t] = (Vsv[t] - Vsv0[t] - Csv_l * (rho * g * Hl -
                        Pext_l[t]) + Csv_u * rho * g * Hu + Cra * P_thorax[t]
                        ) / (Cra + Csv_u + Csv_l)
        else:
            raise ValueError

    Psv_u[t] = max(0, P_ra[t] - rho * g * Hu)  # eq 15
    Psv_l[t] = max(P_ra[t], Pext_l[t]) + rho*g*Hl  # eq 18

    Psa_u[t] = (Vsa[t] - Vsa0[t] + Csa_l*Pext_l[t] - Csa_l * rho * g*H)/(
                Csa)  # eq 12
    Psa_l[t] = (Vsa[t] - Vsa0[t] + Csa_l*Pext_l[t] - Csa_u * rho * g*H)/(
                Csa)  # eq 13
    

    dPsa = (Psa_u[t] - Psa_u_star)
    if control_type == "linear":
        F[t] = get_linear_heart_rate(dPsa, F_star, F_min, Psa_u_star, P_sa_u_min)
        # fit curve through max points and starred points doesn't work. maybe simple point slope rewrite
        # F[t] = get_linear_heart_rate(dPsa, F_star, F_max, Psa_u_star, P_sa_u_max) # not work
    elif control_type == "exp":
        F[t] = get_exp_heart_rate(dPsa, Psa_u_star, F_star)
    
    print("F ", F[t])
    Q[t] = C_RVD * F[t]*(P_ra[t] - P_thorax[t])

    # Rs_l = get_lower_peripheral_resistance(dPsa, Psa_u_star, Rs_l_star,
    #                                        Rs_l_min, P_sa_u_min)
    # Q[t] = 500
    # h is euler iteration timestep
    Vsa[t+1] = Vsa[t] + h * solve_Vsa(Vsa[t], Vsa0[t],Csa, Csa_l, Rs_u,
                                        Rs_l, Csa_u, Pext_l[t], Psa_u[t],
                                        Psv_u[t], Psv_l[t], rho, g, H,
                                        Q[t])

    Vsv[t+1] = Vsv[t] + h * solve_Vsv(Vsa[t], Vsa0[t], Csa, Csa_l, Rs_u,
                                        Rs_l, Csa_u, Pext_l[t], Psa_u[t],
                                        Psv_u[t], Psv_l[t], rho, g, H,
                                        Q[t])
    

    print("vsa+vsv ", Vsa[t+1]+Vsv[t+1])
    # # Case II
    # elif P_thorax[t]> - dP_RA[t] and P_thorax[t]< rho * g * Hu - dP_RA[t]:z

    #     continue
    # # Case III
    # elif P_thorax[t]>= rho * g * Hu - dP_RA[t]:
    #     continue
# breakpoint()


Q = 60/1000 * Q # convert cm3/s to L/min
dynes_2_mmhg = 1/1333
fig = plt.figure(constrained_layout=False)

gs = gridspec.GridSpec(3, 3, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
# identical to ax1 = plt.subplot(gs.new_subplotspec((0, 0), colspan=3))
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])
ax5 = fig.add_subplot(gs[1, 2])
ax6 = fig.add_subplot(gs[0, 2])
ax7 = fig.add_subplot(gs[2, 0])
ax8 = fig.add_subplot(gs[2, 1])
ax9 = fig.add_subplot(gs[2, 2])

start = 0
ax1.plot(T[start:], Vsa[start:-1], label='Vsa')
ax1.set_ylabel("ml")
# ax1.set_title("Vsa")
ax1.legend()
ax2.plot(T[start:], Vsv[start:-1], label='Vsv')
ax2.legend()

# ax2.set_title('Vsv')
ax3.plot(T[start:], Psa_u[start:]*dynes_2_mmhg, label='Psa_u')
# ax3.set_title("Psa_u")
ax3.set_ylabel("mmHg")
ax3.legend()
ax4.plot(T[start:], Psa_l[start:]*dynes_2_mmhg, label ='Psa_l')
# ax4.set_title("Psa_l")
ax4.legend()
ax5.plot(T[start:], P_ra[start:]*dynes_2_mmhg, label='Pra')
ax5.plot(T, rho*g_range*Hu_patient*dynes_2_mmhg, label='rho g Hu')

# ax5.set_title("Pra")
ax5.legend()

ax6.plot(T[start:], Q[start:], label='Q')
ax6.set_ylabel("L/min")
ax6.legend()

ax7.plot(T[start:], F[start:]*60, label='F')
ax7.set_ylabel("B/min")
ax7.legend()

ax8.plot(T, g_range, label='g')
ax8.set_ylabel("x G")
ax8.legend()

ax9.plot(T, rho*g_range*Hu_patient*dynes_2_mmhg, label='rho g Hu')
ax9.set_ylabel("mmHg")
ax9.legend()

fig.tight_layout()

title = control_type + '_control_dynamic_model.png'
plt.savefig(title)

p_range = np.linspace(0, P_sa_u_max)
# f_range = np.array([get_exp_heart_rate(p, F_star, F_min, Psa_u_star, 
#                                    P_sa_u_min) for p in p_range])
# f_range = np.array([get_linear_heart_rate(p, F_star, F_min, Psa_u_star, 
#                                    P_sa_u_min) for p in p_range])
# plt.figure()
# plt.plot(p_range, f_range)
# plt.savefig("P_vs_F_controller.png")