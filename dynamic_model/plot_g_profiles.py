import csv
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

path = Path("dynamic_model/g_profiles")
# fname = path / "1_Gz.csv"
fname = path / "1_Gz.csv"
# with open(fname, newline='') as csvfile:
    # data = list(csv.reader(csvfile))

with open(fname, 'r') as f:
    data = f.readlines()
# df= pd.read_csv(fname)

# df_temp = df.iloc[:,0:2]
# data = df.values.tolist()
# oX = [line.replace('\n', '').split(',')[0] for line in data]
# oY = [float(line.replace('\n', '').split(',')[1]) for line in data]
# print(data)

df = pd.read_csv(fname, usecols=[0,1], header=0)
t= df['t'].to_numpy()
gz = df['Gz'].to_numpy()

plt.figure()
plt.plot(t, gz)
plt.savefig("g_1_profile.png")
print(df)
breakpoint()