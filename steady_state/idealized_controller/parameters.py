"""

Karen's parameters for the model

pressures originally in mmHg, converted to dynes/cm2
1 dyne = 1 g⋅cm/s2
heights in cm
volumes originally in L, converted to cm^3
g in cm/s2
compliances originally in L/mmHg and converted to cm^3/dynes/cm2
resistances in (dynes/cm2)/(cm^3/s) - force/area / flow (vol/time)
Ts, the time constant is in seconds
density is in g/cm3
"""
color_map = 'inferno'
Psa_u_star = 100 * 1333

Psa_u = Psa_u_star
dP_RA = 2 * 1333

# 1/2 factor to average over compartment height (bernoulli's principle - 
# fluid dynamics)

Hu_patient = 0.5*32 # (heart 2 eyeball)
Hl_patient = -0.5*42 # (heart 2 seat)
lumped_height = Hu_patient + (-Hl_patient)

Hu_factor = 1/3
Hl_factor = 2/3

rho = 1

g_earth = 980

# resistance of systemic arteries
Rs = (16.49) * 1333 / (1000 / 60)

Gs = 1 / Rs
Gs_u = Hu_factor * Gs
Gs_l = Hl_factor * Gs
Rs_l = 1 / Gs_l
Rs_u = 1 / Gs_u

# resistance of pulmary arteries
Rp = (1.61 * 1333) / (1000 / 60)

C_RVD = (0.035 / 1333) * 1000

C_LVD = (0.00583 / 1333) * 1000


# systemic compliances 

# literature compliance is probably aortic compliance / diastolic compliance
# and this may be the reason that subdividing works
# if we increase the pressure, how much do the vessels stretch in the upper body and lower body
# legs stretch more because there's more mileage of blood vessels there
# AND there's more blood pooling there because of gravity
# even if each individual vessel is the same level of stretchy, 
# we're averaging over a whole column of fluid (upper and )

#currently set systemic compliances- CSA so that the 0.00175L/mmHG value with Hi/Hu are used for compliances

Csa_l = Hl_factor * (0.00175 / 1333) * 1000 
Csa_u = Hu_factor * (0.00175 / 1333) * 1000

Csv_l = Hl_factor * (0.09 / 1333) * 1000 
Csv_u = Hu_factor * (0.09 / 1333) * 1000 

# Csa_l = (0.000141 / 1333) * 1000
# Csa_u = Csa_l


Cs_l = Csa_l + Csv_l
Csa = Csa_l + Csa_u


# pulmonary compliances
Cpa = (0.00412 / 1333) * 1000
Cpv = (0.01 / 1333) * 1000
Cp = Cpa + Cpv
Vtotal = 3.7 * 1000


# Gs = 1 / Rs_u + 1 / Rs_l
Ts = Csa_u * Rs_u
Tp = Rp * Cpa