#%% -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
from libs import plotconfig as pltconf
from dataclasses import dataclass
from matplotlib.ticker import FuncFormatter

pltconf.plot_params(True, 10, 0.2, fast = True)

def pi_multiple_formatter(x, pos):
    """
    A custom formatter to display tick labels in multiples of π.
    """
    N = int(np.round(x / np.pi))
    if N == 0:
        return "0"
    elif N == 1:
        return r"$\pi$"
    elif N == -1:
        return r"-$\pi$"
    else:
        return r"${0}\pi$".format(N)


@dataclass
class Var:
    var: str
    name: str


save = False

filenum = 1
#system = "pendulum_EMEH"
#system = "lin_oscillator_gravity"
#system = "pend_oscillator_EH"
system = "multidirect_hybrid_EH"
#system = "multidirect_hybrid_EH_zero_len_pend"
ext = ".png"
simulation = "ftimeseries"

df, df_poinc = pltconf.read_CHAOS_data(system, filenum, simulation)
        
dim = 8
nP = 800
nDiv = 6000
trans = 650
plot_i = nDiv*trans

angles = False
#xvar = Var('x[4]', r"$\bar{\phi}$")
yvar = Var('x[7]', r"$\bar{I}$")
xvar = Var('x[5]', r"$\dot{\bar{\phi}}$")
#yvar = Var('x[6]', r"$\bar{v}$")
#xvar = Var('x[0]', r"$\bar{x}$")
#yvar = Var('x[1]', r"$\dot{\bar{x}}$")
#xvar = Var('x[2]', r"$\bar{z}$")
#yvar = Var('x[3]', r"$\dot{\bar{z}}$")

#=======================================================================#
# Figure Parameters                                                     #
#=======================================================================#
figsize = pltconf.figsize_in_cm(15, 0.6*15)
dpi = pltconf.set_fig_quality(save = save, base_dpi = 200)


fig = plt.figure(1, figsize = figsize, dpi = dpi, layout = "constrained")
fig.set_constrained_layout_pads(hspace=0, wspace=0.1)
grid = fig.add_gridspec(2, 4)

ax1 = fig.add_subplot(grid[0,:-2])
ax2 = fig.add_subplot(grid[1,:-2])
ax3 = fig.add_subplot(grid[:,2:4])

size = 1.5

if angles == True:
    ax1.plot(df['Time'], df[f'{xvar.var}_norm'], rasterized = True, color = "red", linewidth = 1, zorder = 1)
    ax1.scatter(df_poinc['Time'], df_poinc[f'{xvar.var}_norm'], rasterized = True, color = "black", s = size, linewidths = 0, zorder = 2)
else:
    ax1.plot(df['Time'], df[f'{xvar.var}'], rasterized = True, color = "red", linewidth = 1, zorder = 1)
    ax1.scatter(df_poinc['Time'], df_poinc[f'{xvar.var}'], rasterized = True, color = "black", s = size, linewidths = 0, zorder = 2)    

ax1.set_ylabel(xvar.name)
ax1.set_xlabel(r'$\tau$')
ax1.set_xlim(df['Time'].min(), df['Time'].max())
#custom_formatter = FuncFormatter(pi_multiple_formatter)
#ax1.yaxis.set_major_formatter(custom_formatter)


ax2.plot(df['Time'], df[f'{yvar.var}'], rasterized = True, color = "blue", linewidth = 1, zorder = 1)
ax2.scatter(df_poinc['Time'], df_poinc[f'{yvar.var}'], rasterized = True, color = "black", s = size, linewidths = 0, zorder = 2)
ax2.set_ylabel(yvar.name)
ax2.set_xlabel(r'$\tau$')
ax2.set_xlim(df['Time'].min(), df['Time'].max())

if angles == True:
    ax3.scatter(df[f'{xvar.var}_norm'].iloc[plot_i:-1], df[f'{yvar.var}'].iloc[plot_i:-1], rasterized = True, color = "darkorange", s = 0.5, linewidths = 0, zorder = 1)
    ax3.scatter(df_poinc[f'{xvar.var}_norm'], df_poinc[f'{yvar.var}'], rasterized = True, color = "black", s = size, linewidths = 0, zorder = 2)
else:
    ax3.plot(df[f'{xvar.var}'].iloc[plot_i:-1], df[f'{yvar.var}'].iloc[plot_i:-1], rasterized = True, color = "darkorange", linewidth = 1, zorder = 1)
    ax3.scatter(df_poinc[f'{xvar.var}'], df_poinc[f'{yvar.var}'], rasterized = True, color = "black", s = size, linewidths = 0, zorder = 2)
ax3.set_ylabel(yvar.name)
ax3.set_xlabel(xvar.name)
#ax3.set_aspect('equal')

#========================================================================#
# Show and Save Figure                                                   #
#========================================================================#

if save == True:
    pltconf.save_CHAOS_data(fig, system, simulation, ext)
    plt.show()
else:
    plt.show()
