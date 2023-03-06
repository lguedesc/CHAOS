#%% -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from libs import plotconfig as pltconf

pltconf.plot_params(True, 10, 0.2)

save = True

simulation = "poinc"
system = "duffing"
ext = ".pdf"
filenum = 0

df = pltconf.read_CHAOS_data(system, filenum, simulation)
#=======================================================================#
# Figure Parameters                                                     #
#=======================================================================#
figsize = pltconf.figsize_in_cm(15, 0.4*15)
dpi = pltconf.set_fig_quality(save = save, base_dpi = 100)

fig = plt.figure(1, figsize = figsize, dpi = dpi, layout = "constrained")
fig.set_constrained_layout_pads(hspace=0, wspace=0.1)
grid = fig.add_gridspec(2, 4)

ax1 = fig.add_subplot(grid[0,:-2])
ax2 = fig.add_subplot(grid[1,:-2])
ax3 = fig.add_subplot(grid[:,2:4])

size = 0.5

ax1.scatter(df['Time'], df['x[0]'], rasterized = True, color = "red", s = size, linewidths = 0, zorder = 1)
ax1.set_ylabel(r'$x$')
ax1.set_xlabel(r'$\tau$')
ax1.set_xlim(df['Time'].min(), df['Time'].max())

ax2.scatter(df['Time'], df['x[1]'], rasterized = True, color = "blue", s = size, linewidths = 0, zorder = 1)
ax2.set_ylabel(r'$\dot{x}$')
ax2.set_xlabel(r'$\tau$')
ax2.set_xlim(df['Time'].min(), df['Time'].max())

ax3.scatter(df['x[0]'], df['x[1]'], rasterized = True, color = "black", s = size, linewidths = 0, zorder = 1)
ax3.set_ylabel(r'$\dot{x}$')
ax3.set_xlabel(r'$x$')
ax3.set_aspect('equal')
ax3.set_ylim(df['x[1]'].min() - abs(df['x[1]'].min()*0.1), df['x[1]'].max() + abs(df['x[1]'].max()*0.1))
ax3.set_xlim(df['x[0]'].min() - abs(df['x[0]'].min()*0.1), df['x[0]'].max() + abs(df['x[0]'].max()*0.1))

#========================================================================#
# Show and Save Figure                                                   #
#========================================================================#

if save == True:
    pltconf.save_CHAOS_data(fig, system, simulation, ext)
    #plt.show()
else:
    plt.show()

