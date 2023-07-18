#%% -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from libs import plotconfig as pltconf

pltconf.plot_params(True, 10, 0.2)

only_steady = True
nP = 500
nDiv = 6000
trans = 400
plot_i = nDiv*trans

angles = True
xvar = 'x[4]'
yvar = 'x[5]'

save = False

system = "pend_oscillator_EH"
simulation = "timeseries"
ext = ".pdf"
filenum = 1

df = pltconf.read_CHAOS_data(system, filenum, simulation)
#=======================================================================#
# Figure Parameters                                                     #
#=======================================================================#
figsize = pltconf.figsize_in_cm(15, 0.4*15)
dpi = pltconf.set_fig_quality(save = save, base_dpi = 150)

fig = plt.figure(1, figsize = figsize, dpi = dpi, layout = "constrained")
fig.set_constrained_layout_pads(hspace=0, wspace=0.1)
grid = fig.add_gridspec(2, 4)

ax1 = fig.add_subplot(grid[0,:-2])
ax2 = fig.add_subplot(grid[1,:-2])
ax3 = fig.add_subplot(grid[:,2:4])

print(df['x[4]'])
print(f"plot_i = {plot_i}")

ax1.plot(df['Time'], df['x[4]'], rasterized = True, color = "red", linewidth = 0.5, zorder = 1)
ax1.set_ylabel(r'$x$')
ax1.set_xlabel(r'$\tau$')
ax1.set_xlim(df['Time'].min(), df['Time'].max())

ax2.plot(df['Time'], df['x[5]'], rasterized = True, color = "blue", linewidth = 0.5, zorder = 1)
ax2.set_ylabel(r'$\dot{x}$')
ax2.set_xlabel(r'$\tau$')
ax2.set_xlim(df['Time'].min(), df['Time'].max())

if only_steady == True:
    if angles == True:
        ax3.scatter(df[f'{xvar}_remainder'].iloc[plot_i:-1], df[f'{yvar}'].iloc[plot_i:-1], rasterized = True, color = "black", linewidths = 0, marker='o', s = 0.2, zorder = 1)
    else:
        ax3.plot(df[f'{xvar}'].iloc[plot_i:-1], df[f'{yvar}'].iloc[plot_i:-1], rasterized = True, color = "black", linewidth = 0.5, zorder = 1)
else:
    if angles == True:
        ax3.scatter(df[f'{xvar}_remainder'], df[f'{yvar}'], rasterized = True, color = "black", linewidths = 0, marker='o', s = 0.2, zorder = 1)
    else:
        ax3.plot(df[f'{xvar}'], df[f'{yvar}'], rasterized = True, color = "black", linewidth = 0.5, zorder = 1)
        
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