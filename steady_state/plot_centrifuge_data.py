"""plot centrifuge data. Gz vs symptoms
    
    1 plot with no strain "2 or 5"

    50% of Gz is 2.2
    100% of Gz is roughly 3.5

    1 plot with hook only - model intrathoracic pressure

    (1 plot with legs +hook to motivate future model.
    maybe don't include because legs and hook because 
    reviewers will comment on it)

    us = "ultrasound"

    legs + hook, hr goes lower than hook only 
    and lower on average than legs only

    the data point at 180 is from a winged vehicle profile
    that has Gz as well and was done possibly before the
    pure Gz (so excitement)

    data point at 180 is mixed gz gx (boost)
    data point at 108 is reentry (almost all gz)
    
    data point at 100 HR is second (different) winged 
    vehicle profile - also from re-rentry
    

    feather is when the vehicle folds in half
    look up virgin galactic 
    winged vehicles 1 and 2 are different planes

    ultrasound on all data plot : relaxed , regular profile (non winged vehicle)

"""

import numpy as np
import matplotlib.pyplot as plt

gz_values = np.ones((1,11))
gz_values[0][6:8] *= 2.2
gz_values[0][8:11] *= 3.5

hr_values = np.array([[73, 77, 63, 55, 64, 78, 77, 140, 100, 112, 92]])
gz_relaxed_table = np.vstack((gz_values, hr_values))