#%% -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from libs import plotconfig as pltconf

pltconf.plot_params(True, 10, 0.2)

only_steady = False
nP = 5000
nDiv = 1000
trans = 3750
plot_i = nP*trans

save = True

system = "duffing"
simulation = "timeseries"
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

ax1.plot(df['Time'], df['x[0]'], rasterized = True, color = "red", linewidth = 0.5, zorder = 1)
ax1.set_ylabel(r'$x$')
ax1.set_xlabel(r'$\tau$')
ax1.set_xlim(df['Time'].min(), df['Time'].max())

ax2.plot(df['Time'], df['x[1]'], rasterized = True, color = "blue", linewidth = 0.5, zorder = 1)
ax2.set_ylabel(r'$\dot{x}$')
ax2.set_xlabel(r'$\tau$')
ax2.set_xlim(df['Time'].min(), df['Time'].max())

if only_steady == True:
    ax3.plot(df['x[0]'].iloc[plot_i:-1], df['x[1]'].iloc[plot_i:-1], rasterized = True, color = "black", linewidth = 0.5, zorder = 1)
    ax3.set_ylabel(r'$\dot{x}$')
    ax3.set_xlabel(r'$x$')
else:
    ax3.plot(df['x[0]'], df['x[1]'], rasterized = True, color = "black", linewidth = 0.5, zorder = 1)
    ax3.set_ylabel(r'$\dot{x}$')
    ax3.set_xlabel(r'$x$')
#========================================================================#
# Show and Save Figure                                                   #
#========================================================================#
if save == True:
    pltconf.save_CHAOS_data(fig, system, simulation, ext)
    #plt.show()
else:
    plt.show()