"""code to plot data points from centrifuge and model prediction on the 
same plot
"""

import numpy as np
import matplotlib.pyplot as plt
from parameters import *
from matplotlib import rc
import matplotlib as mpl
from matplotlib.lines import Line2D

mpl.rcParams['mathtext.fontset'] = 'cm'

from run_model import plot_F_vs_G_case_I
from plot_centrifuge_data import plot_centrifuge_data

def main():
    fig, ax = plt.subplots()

    gz_list, hr_list = plot_centrifuge_data()
    G, sol_F_Vtotal_G = plot_F_vs_G_case_I()

    plt.plot(G, sol_F_Vtotal_G, label="model prediction")
    plt.plot(gz_list, hr_list, 'ro', label="data points")
    
    plt.axis((0, 4, 40, 200))
    plt.legend(loc="upper right")

    plt.savefig("combined_data.png")

if __name__ == "__main__":
    main()