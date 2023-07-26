# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt    
import matplotlib.colors as mpl_col
from matplotlib.ticker import FormatStrFormatter

import os
from libs import plotconfig as pltconf

pltconf.plot_params(False, 10, 0.5)

# =========================================================================== #
#                           Function to change axis labels 
# =========================================================================== #
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
# =========================================================================== #
#                                    Read Data                                #
# =========================================================================== #
maxper = 6
save = False
system = "duffing"
ext = ".pdf"

readpath = "data/FDynDiagram/out/" + system + "_fdyndiag.csv"; readpath = pltconf.convert_dir(readpath)
savepath = "data/FDynDiagram/figs"; savepath = pltconf.convert_dir(savepath)

raw_data = pd.read_csv(readpath, delimiter = " ")

x1, y1, z1 = pltconf.process_data(raw_data, 'CparY', 'CparX', 'Attractor')
x2, y2, z2 = pltconf.process_data(raw_data, 'CparY', 'CparX', 'xMAX[0]')
x3, y3, z3 = pltconf.process_data(raw_data, 'CparY', 'CparX', 'xMIN[0]')
#x4, y4, z4 = pltconf.process_data(raw_data, 'CparY', 'CparX', 'ddxbRMS')
#x5, y5, z5 = pltconf.process_data(raw_data, 'CparY', 'CparX', 'TotalPout')
#x6, y6, z6 = pltconf.process_data(raw_data, 'CparY', 'CparX', 'xRMS[4]')
#x7, y7, z7 = pltconf.process_data(raw_data, 'CparY', 'CparX', 'xRMS[5]')
#x8, y8, z8 = pltconf.process_data(raw_data, 'CparY', 'CparX', 'xRMS[2]')
# =========================================================================== #
#                           Define figure parameters                          #
# =========================================================================== #
cm = 1/2.54
x_inches = 15     # [mm]*constant
y_inches = x_inches*(1.2)
raster = True
if save == True:
    dpi = 300
else:
    dpi = 96
# =========================================================================== #
#                               Create figure 
# =========================================================================== #
fig = plt.figure(1, figsize = (x_inches*cm,y_inches*cm), dpi = dpi)#, constrained_layout = True)
#fig.tight_layout()
#fig.set_constrained_layout_pads(hspace = 0, wspace = 0.1, w_pad = 0, h_pad = 0)


fig.subplots_adjust(top=1,
                    bottom=0.08,
                    left=0.063,
                    right=0.927,
                    hspace=0.1,
                    wspace=0.26)


lin = 4; col = 2
ax1 = fig.add_subplot(lin,col,1)
ax2 = fig.add_subplot(lin,col,2)
ax4 = fig.add_subplot(lin,col,3)
ax3 = fig.add_subplot(lin,col,4)
ax5 = fig.add_subplot(lin,col,5)
ax6 = fig.add_subplot(lin,col,6)
ax7 = fig.add_subplot(lin,col,7)
ax8 = fig.add_subplot(lin,col,8)
# =========================================================================== #
#                               Plot Data  
# =========================================================================== #
pltconf.plot_attractor_map(fig, ax1, x1, y1, z1, maxper)
pltconf.plot_neg_pos_map(ax2, x2, y2, z2)
pltconf.plot_neg_pos_map(ax3, x3, y3, z3)
#pltconf.plot_rainbow_map(ax4, x4, y4, z4)
#pltconf.plot_rainbow_map(fig, ax5, x5, y5, z5)
#pltconf.plot_rainbow_map(ax6, x6, y6, z6)
#pltconf.plot_rainbow_map(ax7, x7, y7, z7)
#pltconf.plot_rainbow_map(ax8, x8, y8, z8)
# =========================================================================== #
#                          Customize titles and labels                        #
# =========================================================================== #
customize_labels(ax1, r'(a) Dynamical Attractors', custom = '(a)')
customize_labels(ax2, r'(b) xMAX[0]', custom = '(b)')
customize_labels(ax3, r'(d) xMIN[0]',  custom = '(b)')
customize_labels(ax4, r'(c) ', custom = '(a)')
customize_labels(ax5, r'(e) TotalPout', custom = '(a)')
customize_labels(ax6, r'(f) xRMS[4]', custom = '(b)')
customize_labels(ax7, r'(g) xRMS[5]', custom = '(c)')
customize_labels(ax8, r'(h) ', custom = '(d)')
# =========================================================================== #
#                                Save Figure                                  #
# =========================================================================== #
if save == True:
    isExist = os.path.exists(savepath)
    if (isExist == False):
        os.makedirs(savepath)
    
    name = "/sample_" + system + "_dyndiag(serial)" + ext; name = pltconf.convert_dir(name)
    fig.savefig(savepath + name)
    plt.show()
else:
    plt.show()

