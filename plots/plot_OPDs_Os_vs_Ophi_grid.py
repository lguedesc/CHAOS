# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from libs import plotconfig as pltconf
from matplotlib.patches import Rectangle

pltconf.plot_params(True, 9, 0.2)

def mark_point_in_plot(ax, x, y, number, color, mode = '1'):
    if mode == '1':
        ax.scatter(x, y, s = 5, linewidth = 0.5, facecolors = 'none', edgecolors = color)
        ax.text(x + 0.11, y - 0.08 , number, color = color, fontsize=6)
    else:
        ax.scatter(x, y, s = 5, linewidth = 0.5, facecolors = 'none', edgecolors = color)
        ax.text(x + 0.06, y - 0.04, number, color = color, fontsize=6)
        
#=======================================================================#
# General Input
#=======================================================================#
save = True
maxper = 6

#power = 'PoutPZ_Avg'
#power = 'PoutEM_Avg'
power = 'TotalPout'
#=======================================================================#
# Make Figure
#=======================================================================#
rows = 4; cols = 4
if save == True:
    dpi = 600
else:
    dpi = 200
figsize = pltconf.figsize_in_cm(15, 0.27*rows*15)

fig, axs = plt.subplots(rows, cols , figsize=figsize, dpi = dpi, layout='constrained')

Omegas = [0.25, 0.5, 1.0, 1.5]
gammas = [0.1, 0.3, 0.5, 0.7]
letters = ['a', 'b', 'c', 'd']

limited_axs = [axs[1,3], axs[2,2], axs[2,3], axs[3,2], axs[3,3]]

for i in range(len(axs[:, 0])):        
    for ax, O in zip(axs[i, :], Omegas):
        #=======================================================================#
        # Read and Process Diagram Data
        #=======================================================================#
        path = "/Volumes/Luã SSD/Hybrid-Multidirectional Energy Harvester/Omega_s vs Omega_phi/"
        pathname = f"{path}pend_oscillator_EH_dyndiag_O={O:.2f}_g={gammas[i]:.1f}.csv"
        pathname = pltconf.convert_dir(pathname)

        raw_data = pltconf.read_data(pathname)
        x, y, z = pltconf.process_data(raw_data, 'CparY', 'CparX', power)
        z = z*1e3
        cbarmax = z.max()
        #=======================================================================#
        # Plot Dynamical Diagrams 
        #=======================================================================#
        if power == 'TotalPout' or power == 'PoutEM_Avg':
            if ax in limited_axs:
                factor = 0.2
            else:
                factor = 1.0
            
        Xmargin = 0.033   #0.022
        Ymargin = 0.032  #0.065
        # Full Attractor Diagram
        plot, cbar = pltconf.plot_rainbow_map(fig, ax, x, y, z, cbarmax = cbarmax*factor)
        ax.use_sticky_edges = False
        ax.margins(x=Xmargin, y=Ymargin)
        
        if power == 'TotalPout' or power == 'PoutEM_Avg':
            if ax in limited_axs:
                cbar.ax.set_title(rf"${cbarmax:.2f}$", fontsize = 6, loc = "left")
#'''
#========================================================================#
# Define Labels and ticklabels
#========================================================================#
for ax in axs.flat:
    ax.set_xticks([0.01, 1, 2])
    ax.set_yticks([0.01, 1, 2])
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))

for ax in axs[:,0]:
    ax.set_ylabel(r"$\Omega_{\phi}$")    

for ax in axs[-1,:]:
    ax.set_xlabel(r"$\Omega_{s}$")

for ax_row in axs[:, 1:]:
    for ax in ax_row:
        ax.yaxis.set_ticklabels([])
    
for ax_col in axs[:-1, :]:
    for ax in ax_col:
        ax.xaxis.set_ticklabels([])

for i in range(len(axs[:, 0])):        
    for ax, Os in zip(axs[i, :], Omegas):
        ax.set_title(rf"$\Omega = {Os}$", loc = 'center')

for ax, ltr, g in zip(axs[:, 0], letters, gammas):
    ax.set_title(rf"({ltr}) $\gamma = {g}$" + "\n\n", loc = 'left')

#'''
#========================================================================#
# Show and Save Figure                                                   #
#========================================================================#
if save == True:
    name = f"OPDs_grid_{power}.pdf"; name = pltconf.convert_dir(name)
    fig.savefig(name)
#    plt.show()
else:
    plt.show()
