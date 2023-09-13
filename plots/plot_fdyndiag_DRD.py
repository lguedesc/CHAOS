# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt    
import matplotlib.colors as mpl_col
from matplotlib.ticker import FormatStrFormatter

import os
from libs import plotconfig as pltconf

pltconf.plot_params(False, 10, 0.2)

# =========================================================================== #
#                           Function to change axis labels 
# =========================================================================== #
'''
def customize_labels(ax, title1, custom = False):
    ax.set_title(title1, loc = 'left')
    #ax.set_title(title2, loc = 'center')
    if custom == '(a)':
        ax.set_xticks([0.01, 1, 2])
        ax.axes.xaxis.set_ticklabels([])
        ax.set_yticks([y2.min(), y2.max()])
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.set_ylabel(r'$\gamma$', labelpad = -10)
    elif custom == '(b)':
        ax.set_xticks([0.01, 1, 2])
        ax.axes.xaxis.set_ticklabels([])
        ax.set_yticks([y2.min(), y2.max()])
        ax.axes.yaxis.set_ticklabels([])
        #ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        #ax.set_ylabel(r'$\gamma$', labelpad = -10)
    elif custom == '(c)':
        ax.set_xlabel(r'$\Omega$')
        ax.set_xticks([0.01, 1, 2])
        ax.set_yticks([y2.min(), y2.max()])
        ax.set_ylabel(r'$\gamma$', labelpad = -10)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    elif custom == '(d)':
        ax.set_xlabel(r'$\Omega$')
        ax.set_xticks([0.01, 1, 2])
        ax.set_yticks([y2.min(), y2.max()])
        ax.axes.yaxis.set_ticklabels([])
        ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
'''
# =========================================================================== #
#                                    Read Data                                #
# =========================================================================== #
maxper = 6
save = False
folder = "data/FDynDiagram/"
system = "duffing"
number = 8
ext = ".pdf"

if number == 0:
    readpath = f"{folder}out/" + system + f"_fdyndiag.csv"; readpath = pltconf.convert_dir(readpath)
else:
    readpath = f"{folder}out/" + system + f"_fdyndiag({number}).csv"; readpath = pltconf.convert_dir(readpath)
savepath = f"{folder}figs"; savepath = pltconf.convert_dir(savepath)

raw_data = pd.read_csv(readpath, delimiter = " ")

x, y, z = pltconf.process_data(raw_data, 'CparY', 'CparX', 'Attractor')
# =========================================================================== #
#                           Define figure parameters                          #
# =========================================================================== #
figsize = pltconf.figsize_in_cm(15, 0.6*15)
raster = True
if save == True:
    dpi = 300
else:
    dpi = 100
# =========================================================================== #
#                               Create figure and axes
# =========================================================================== #
fig, axs = pltconf.makefig_and_axs(figsize, 1, 1, dpi, hspace = 0.1, wspace = 0.1)
# =========================================================================== #
#                               Plot Data  
# =========================================================================== #
pltconf.plot_attractor_map(fig, axs, x, y, z, maxper)
# =========================================================================== #
#                          Customize titles and labels                        #
# =========================================================================== #
# =========================================================================== #
#                                Save Figure                                  #
# =========================================================================== #
pltconf.savefigure(fig, save, savepath, system, ext)
plt.show()

