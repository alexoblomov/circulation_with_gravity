"""
Purpose : simulate ode model

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from parameters import *


def volume_odes(x,t):
    Qrv = x[0]
    Qlv = x[1]
    Qsa_l = x[2]
    Qsa_u = x[3]
    Qsv_l = x[4]
    Qsv_u = x[5]

    dVstroke_dt = Qrv - Qlv
    dVsa_udt = Qlv - Qsa_l
    dVsa_ldt = Qsa_u - Qsv_l
    dVsv_ldt = Qsa_l - Qsv_u
    dVsv_udt = Qsv_l -Qrv

    return [dVstroke_dt, dVsa_udt, dVsa_ldt, dVsv_ldt, dVsv_udt]

# initial conditions - wrong need flows not volume ??
# oVstroke = 70 #ml
# artery_vein_split = 0.5
# oVsa_u = artery_vein_split * Hu_factor * Vtotal
# oVsv_u = artery_vein_split * Hu_factor * Vtotal
# oVsa_l = artery_vein_split * Hl_factor * Vtotal
# oVsv_l = artery_vein_split * Hl_factor * Vtotal

# x0 = [oVstroke, oVsa_u, oVsa_l, oVsv_l, oVsv_u]

Qrv = 500
Qlv = 500
Qsa_l = 500
Qsa_u = 500
Qsv_l = 500
Qsv_u = 500
x0 = 6*[Qrv]
test = volume_odes(x0, 0)

# timestep
t = np.linspace(0,3*60)

# SOLVE
x = odeint(volume_odes, x0, t)
breakpoint()
# plot
