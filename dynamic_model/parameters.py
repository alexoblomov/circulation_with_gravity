"""
Karen's parameters for the model

pressures originally in mmHg, converted to dynes/cm2
heights in cm
volumes originally in L, converted to cm^3
g in m/s2
compliances originally in L/mmHg and converted to cm^3/dynes/cm2
resistances in (dynes/cm2)/(cm^3/s)
Ts, the time constant is in seconds
"""
Psa_u_star = 100 * 1333

Psa_u = Psa_u_star
dP_RA = 2 * 1333
height = 167.64

Hu_patient = 32
Hl_patient = 42
lumped_height = Hu_patient + Hl_patient

Hu_factor = Hu_patient /lumped_height
Hl_factor = Hl_patient/lumped_height

rho = 1

g_earth = 980

Rs = (16.49) * 1333 / (1000 / 60)

Gs = 1 / Rs
Gs_u = Hu_factor * Gs
Gs_l = Hl_factor * Gs
Rs_l = 1 / Gs_l
Rs_u = 1 / Gs_u
Rp = (1.61 * 1333) / (1000 / 60)

C_RVD = (0.035 / 1333) * 1000

C_LVD = (0.00583 / 1333) * 1000

Csa_l = Hl_factor * (0.00175 / 1333) * 1000
Csa_u = Hu_factor * (0.00175 / 1333) * 1000
Csv_l = Hl_factor * (0.09 / 1333) * 1000
Csv_u = Hu_factor * (0.09 / 1333) * 1000
Cs_l = Csa_l + Csv_l
Cpa = (0.00412 / 1333) * 1000
Cpv = (0.01 / 1333) * 1000
Cp = Cpa + Cpv
#  RANDOM GUESS. # TODO: change to real value
Cra = (0.01 / 1333) * 1000
Vtotal = 3.7 * 1000

Csa = Csa_l + Csa_u
# Gs = 1 / Rs_u + 1 / Rs_l
Ts = Csa_u * Rs_u

Tp = Rp * Cpa
Csa = Csa_u + Csa_l

F_patient = 60 / 60 # beats per second

# at G_earth
# VT0_steady_state_cntrl = 1.6785631 * 1000
VT0_steady_state_cntrl = 0.01 * 1000

# else need to simulate and input VTO array vs G