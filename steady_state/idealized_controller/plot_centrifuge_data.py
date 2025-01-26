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

 #   ultrasound on all data plot : relaxed , regular profile (non winged vehicle)
 # when running the below participants- parameters file was adjusted appropriately for Hu, Hi, and V total
 # values have been returned back to the original participant for saving
 #
# Initial participant HR values 1Gz: 73,77,63,55,64,78; 2.2Gz: 77; 3.5Gz 100, 112,92
# participant #2 HR values 1Gz: 108, 114, 103, 120, 93, 96, 91, 95; 2.2Gz: 120; 3.5Gz 142
# participant #3 HR values 1Gz: 90, 92, 103, 120, 120, 93, 96, 91; 2.2Gz: 122; 3.5Gz: 142
# participant #4 HR values 1Gz, 98, 100, 96, 90, 98, 100, 100, 85; 2.2Gz:108, 3.5Gz: 152
# participant #5 HR values 1Gz: 82, 91, 90, 82, 82, 90, 80, 87; 2.2Gz: 94; 3.5Gz: 109
# participant #6 HR values 1Gz: 92, 89, 84, 83, 91, 75, 110, 76; 2.2Gz: 117; 3.5Gz: 139
# participant #7 HR values 1Gz: 98, 86, 110, 94; 2.2 Gz: 124, 122, 138; 3.5Gz: 138, 152, 166
# participant #8 HR values 1Gz: 90, 102, 110; 2.2Gz: 120; 3.5Gz: 125, 110, 110
# participant #9 HR values 1Gz: 125, 125, 130, 124; 2.2Gz: 135; 3.5Gz: 155, 155
"""

import numpy as np
import matplotlib.pyplot as plt

def plot_centrifuge_data():
    gz_values = np.ones((1,10))
    gz_values[0][7] *= 2.2
    gz_values[0][7:10] *= 3.5

    hr_values = np.array([[73, 77, 63, 55, 64, 78, 77, 100, 112, 92]])
    gz_relaxed_table = np.vstack((gz_values, hr_values))

    # as 1D arrays
    # can add in the two blue data points 62 and 65 at idle (gz=1)
    gz_list = np.ones((10))
    gz_list[6] *= 2.2
    gz_list[7:10] *= 3.5
    hr_list = np.array([73, 77, 63, 55, 64, 78, 77, 100, 112, 92])

    fig, ax = plt.subplots()
    x = ax.plot(gz_list, hr_list, 'ro')
    # plt.axis((0, 4, 40, 200))
    # plt.savefig("relaxed_gz_HR_centrifuge_data.png")

    # plot hook only 
    # 3 data points

    gz_hook = np.ones((3)) * 3.5
    hr_hook = np.array([142, 152, 126])

    # currently adds all data points to plot
    # may want to clear relaxed hr and have hook only on the plot
    # fig = plt.figure()
    # plt.axis((0, 4, 40, 200))
    # plt.plot(gz_hook, hr_hook, 'bo')
    # plt.axis((0, 4, 40, 200))
    # plt.savefig("hook_maneuver_gz_HR_centrifuge_data.png")
    return (gz_list, hr_list)

# if __name__ == "__main__":
#     plot_centrifuge_data()