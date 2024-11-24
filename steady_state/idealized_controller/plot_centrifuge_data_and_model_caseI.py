"""code to plot data points from centrifuge and model prediction on the 
same plot
"""

import numpy as np
import matplotlib.pyplot as plt
from parameters import *
from matplotlib import rc
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'

from run_model import plot_F_vs_G_case_I
from plot_centrifuge_data import plot_centrifuge_data

def main():
    fig, ax = plt.subplots()
    # plt.axis((0, 4, 40, 200))
    # plot_F_vs_G_case_I()
    gz_list, hr_list = plot_centrifuge_data()
    G, sol_F_Vtotal_G = plot_F_vs_G_case_I()
    plt.plot(G, sol_F_Vtotal_G)
    plt.plot(gz_list, hr_list, 'ro')
    
    # breakpoint()
    plt.axis((0, 4, 40, 200))

    plt.savefig("combined_data.png")

    # breakpoint()

if __name__ == "__main__":
    main()