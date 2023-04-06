import pandas as pd
import matplotlib.pyplot as plt


def import_g_profile(fname):

    df = pd.read_csv(fname, header=0)
    t= df['t'].to_numpy()
    gz = df['Gz'].to_numpy()
    return t, gz