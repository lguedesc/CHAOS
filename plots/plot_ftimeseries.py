#%% -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import os
from libs import plotconfig as pltconf

pltconf.plot_params(True, 10, 0.2)

save = False

system = "bistable_EH"
ext = ".png"
filenum = 0
simulation = "ftimeseries"

df, df_poinc = pltconf.read_CHAOS_data(system, filenum, simulation)
        
nP = 1000
nDiv = 1000
trans = 750
plot_i = nP*trans

#=======================================================================#
# Figure Parameters                                                     #
#=======================================================================#
figsize = pltconf.figsize_in_cm(15, 0.6*15)
dpi = pltconf.set_fig_quality(save = save, base_dpi = 200)


fig = plt.figure(1, figsize = figsize, dpi = dpi, layout = "constrained")
fig.set_constrained_layout_pads(hspace=0, wspace=0.1)
grid = fig.add_gridspec(3, 4)

ax1 = fig.add_subplot(grid[0,:-2])
ax2 = fig.add_subplot(grid[1,:-2])
ax3 = fig.add_subplot(grid[0:2,2:4])
ax4 = fig.add_subplot(grid[2,:-2])

size = 1.5

ax1.plot(df['Time'], df['x[0]'], rasterized = True, color = "red", linewidth = 1, zorder = 1)
ax1.scatter(df_poinc['Time'], df_poinc['x[0]'], rasterized = True, color = "black", s = size, linewidths = 0, zorder = 2)
ax1.set_ylabel(r'$x$')
ax1.set_xlabel(r'$\tau$')
ax1.set_xlim(df['Time'].min(), df['Time'].max())

ax2.plot(df['Time'], df['x[1]'], rasterized = True, color = "blue", linewidth = 1, zorder = 1)
ax2.scatter(df_poinc['Time'], df_poinc['x[1]'], rasterized = True, color = "black", s = size, linewidths = 0, zorder = 2)
ax2.set_ylabel(r'$\dot{x}$')
ax2.set_xlabel(r'$\tau$')
ax2.set_xlim(df['Time'].min(), df['Time'].max())

ax3.plot(df['x[0]'].iloc[plot_i:-1], df['x[1]'].iloc[plot_i:-1], rasterized = True, color = "black", linewidth = 1, zorder = 1)
ax3.scatter(df_poinc['x[0]'], df_poinc['x[1]'], rasterized = True, color = "orange", s = size, linewidths = 0, zorder = 2)
ax3.set_ylabel(r'$\dot{x}$')
ax3.set_xlabel(r'$x$')
#ax3.set_aspect('equal')

ax4.plot(df['Time'], df['LE[0]'], rasterized = True, color = "green", linewidth = 1, zorder = 1, label = "$\lambda_1$")
ax4.plot(df['Time'], df['LE[1]'], rasterized = True, color = "purple", linewidth = 1, zorder = 1, label = "$\lambda_2$")
ax4.hlines(0, df['Time'].min(), df['Time'].max(), color = "black", linewidth = 0.5, zorder = 3)
ax4.set_ylabel(r'$\lambda$')
ax4.set_xlabel(r'$\tau$')
ax4.set_xlim(df['Time'].min(), df['Time'].max())
ax4.legend(loc = "best", prop={'size': 6})

#========================================================================#
# Show and Save Figure                                                   #
#========================================================================#

if save == True:
    pltconf.save_CHAOS_data(fig, system, simulation, ext)
    plt.show()
else:
    plt.show()
