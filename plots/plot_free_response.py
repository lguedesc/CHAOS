import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import os
from libs import plotconfig as pltconf
from matplotlib.ticker import FormatStrFormatter
from math import remainder
import numpy as np

def max_absolute_value(v1, v2):
    if abs(v1) > abs(v2):
        return v1
    elif abs(v2) > abs(v1):
        return v2
    else:
        return v1
    
def format_label_range(maxabs):
    # Get nearest whole number related to maxabs
    roundednum = round(maxabs)
    # Check if it is smaller than the original number
    if roundednum < maxabs:
        roundednum = roundednum + 1
    return roundednum

pltconf.plot_params(True, 10, 0.2)

save = False
init_cond = "phi"

#=======================================================================#
# Load Data                                                             #
#=======================================================================#
filenum = 11
directory = "/Users/luaguedescosta/Desktop/CHAOS/data/FTimeSeries/out/"
savepath = f"{directory}figs/"
readpath = f"{directory}multidirect_hybrid_EH_ftimeseries({filenum}).csv"
df = pd.read_csv(readpath, delimiter = " ")
#=======================================================================#
# Figure Parameters                                                     #
#=======================================================================#
rows = 1; cols = 3
width = 15
dpi = pltconf.set_fig_quality(save = save, base_dpi = 200)
figsize = pltconf.figsize_in_cm(width, 0.15*width)
fig, axs = pltconf.makefig_and_axs(figsize, rows, cols, dpi = dpi, wspace = 0.02, hpad = 0.02)
#=======================================================================#
# Plot Data                                                             #
#=======================================================================#
lwd = 0.7

labels = [r"$\bar{x}$", r"$\bar{z}$", r"$\bar{\phi}$"]
colors = ['black', 'darkorange', 'red']
variables = ['x[0]', 'x[2]', 'x[4]']

for ax, lbl, clr, var in zip(axs, labels, colors, variables):
    ax.plot(df['Time'], df[var], color = clr, lw = lwd)
    ax.set_ylabel(lbl)
    ax.set_xlabel(r"$\tau$")
    
    # Find maximum absolute value
    maxval = df[var].max()
    minval = df[var].min()
    maxabs = max_absolute_value(maxval, minval)
     
    ax.set_ylim(-maxabs - maxabs*0.1, maxabs + maxabs*0.1)
    ax.set_xlim(df['Time'].min(), round(df['Time'].max()))
    ax.set_xticks([df['Time'].min(), round((df['Time'].max() - df['Time'].min())/2) , round(df['Time'].max())])
    #ax.set_yticks([-lblrange, lblrange])
    
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useMathText=True)

#axs[0].set_yticks([-4e-5, 0, 4e-5])
#========================================================================#
# Show and Save Figure                                                   #
#========================================================================#
pltconf.savefigure(savepath, f"free_response_{init_cond}", ".pdf", fig, save = save)
if save == False:
    plt.show()
