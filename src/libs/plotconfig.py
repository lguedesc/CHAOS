import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from matplotlib.colors import BoundaryNorm, ListedColormap
import matplotlib.colors as mpl_col

def plot_params(tex, fontsize, lineweight):
    mpl.rcParams['agg.path.chunksize'] = 10000
    plt.rcParams.update({
        'font.family': 'serif',  # use serif/main font for text elements
        'text.usetex': tex,     # use inline math for ticks
        'axes.titlesize': 'medium',
        'font.size': fontsize,
        'axes.linewidth' : lineweight,
        'xtick.major.width': lineweight,
        'ytick.major.width': lineweight,
        'xtick.major.size': 1.75,
        'ytick.major.size': 1.75,
        'xtick.top': True,
        'ytick.right': True,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
    })

def convert_dir(dir):
    if os.name == 'nt':      # For Windows
        if (dir.find('/') > 0):
            dir.replace('/', '\\')
    else:                   # For Linux and MacOs
        if (dir.find('\\') > 0):
            dir.replace('\\', '/')
    return dir

def set_colormap(df):
    if (df.max() >= 8):
        colormap = ListedColormap(['#404040', '#FFEE00', '#00DC00', '#FF8000', '#9900FF', '#007BFF', '#FF0000', '#700000'])
    else:
        colormap = ListedColormap(['#404040', '#FFEE00', '#00DC00', '#FF8000', '#9900FF', '#007BFF', '#FF0000'])
    return colormap

# ====== IN PROGRESS FUNCTIONS ========= #
def handle_s_lyap(df, trans, columnName, dim):
    maxtime = trans*df[columnName].max()
    for i in range(dim):
        df["sLE["+str(i)+"]"].loc[df[columnName] < maxtime] = np.nan
    print(df[columnName].loc[:maxtime])