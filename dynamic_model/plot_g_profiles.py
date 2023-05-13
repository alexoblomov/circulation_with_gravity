import csv
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from utils.io import import_g_profile

path = Path("dynamic_model/g_profiles")
fname = path / "Run_1_G.csv"

T, g_range = import_g_profile(fname)


plt.figure()
plt.plot(T, g_range)
plt.savefig("g_1_profile.png")
