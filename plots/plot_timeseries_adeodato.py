#%% -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import os
from libs import plotconfig as pltconf

pltconf.plot_params(True, 10, 0.5)

nP = 200
nDiv = 1000
trans = 150
plot_i = nDiv*trans

save = False

system = "adeodato_sma_oscillator"
ext = ".pdf"

readpath = "TimeSeries/out/" + system + "_rk4(10).csv"; readpath = pltconf.convert_dir(readpath)
savepath = "TimeSeries/figs"; savepath = pltconf.convert_dir(savepath)
        
df = pd.read_csv(readpath, delimiter = " ")

#=======================================================================#
# Figure Parameters                                                     #
#=======================================================================#
x_inches2 = 150*(1/25.4)     # [mm]*constant
y_inches2 = x_inches2*(0.4)

if save == True:
    dpi = 300
else:
    dpi = 200

fig = plt.figure(1, figsize = (x_inches2,y_inches2), dpi = dpi, constrained_layout = True)
fig.set_constrained_layout_pads(hspace=0, wspace=0.1)
grid = fig.add_gridspec(2, 4)

ax1 = fig.add_subplot(grid[0,:-2])
ax2 = fig.add_subplot(grid[1,:-2])
ax3 = fig.add_subplot(grid[:,2:4])

ax1.plot(df['Time'], df['x[0]'], rasterized = False, color = "red", linewidth = 0.5, zorder = 1)
ax1.set_ylabel(r'$x$')
ax1.set_xlabel(r'$\tau$')
ax1.set_xlim(9.4, df['Time'].max())
#ax1.set_xlim(df['Time'].min(), df['Time'].max())

ax2.plot(df['Time'], df['dir'], rasterized = True, color = "blue", linewidth = 0.5, zorder = 1)
ax2.plot(df['Time'], df['dir_ant'], rasterized = True, color = "red", linewidth = 0.5, zorder = 1)
ax2.set_ylabel(r'$dir$')
ax2.set_xlabel(r'$\tau$')
ax2.set_xlim(9.4, df['Time'].max())
#ax2.set_xlim(df['Time'].min(), df['Time'].max())

ax3.plot(df['x[0]'].iloc[plot_i:-1], df['x[1]'].iloc[plot_i:-1], rasterized = False, color = "black", linewidth = 0.5, zorder = 1)
ax3.set_ylabel(r'$\dot{x}$')
ax3.set_xlabel(r'$x$')
#ax3.set_aspect('equal')


#========================================================================#
# Show and Save Figure                                                   #
#========================================================================#

if save == True:
    isExist = os.path.exists(savepath)
    if (isExist == False):
        os.makedirs(savepath)
    
    name = "/sample_" + system + ext; name = pltconf.convert_dir(name)
    fig.savefig(savepath + name)
    plt.show()
else:
    plt.show()

