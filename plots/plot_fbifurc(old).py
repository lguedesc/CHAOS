# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import BoundaryNorm, ListedColormap
import pandas as pd
import os
from libs import plotconfig as pltconf

pltconf.plot_params(True, 10, 0.5)

save = False
system = "bistable_EH"
ext = ".pdf"
num = 2

path = "data/FBifurcation/"

if num > 0:
    readpath = f"{path}out/" + system + f"_fbifurc({num}).csv"; readpath = pltconf.convert_dir(readpath)
    readpathpoinc = f"{path}out/" + system + f"_fbifurc_poinc({num}).csv"; readpathpoinc = pltconf.convert_dir(readpathpoinc)
    savepath = f"{path}figs"; savepath = pltconf.convert_dir(savepath)
else:
    readpath = f"{path}out/" + system + "_fbifurc.csv"; readpath = pltconf.convert_dir(readpath)
    readpathpoinc = f"{path}out/" + system + "_fbifurc_poinc.csv"; readpathpoinc = pltconf.convert_dir(readpathpoinc)
    savepath = f"{path}figs"; savepath = pltconf.convert_dir(savepath)
    
df = pd.read_csv(readpath, delimiter = " ")
dfpoinc = pd.read_csv(readpathpoinc, delimiter = " ")
#=======================================================================#
# Figure Parameters                                                     #
#=======================================================================#
x_inches2 = 150*(1/25.4)     # [mm]*constant
y_inches2 = x_inches2*(0.6)

if save == True:
    dpi = 300
else:
    dpi = 200

fig = plt.figure(1, figsize = (x_inches2,y_inches2), dpi = dpi, constrained_layout = True)
fig.set_constrained_layout_pads(hspace=0, wspace=0.1)
grid = fig.add_gridspec(4, 4)

ax1 = fig.add_subplot(grid[0,:-2])
ax2 = fig.add_subplot(grid[1,:-2])
ax3 = fig.add_subplot(grid[2,:-2])
ax4 = fig.add_subplot(grid[3,:-2])
ax5 = fig.add_subplot(grid[0,2:4])
ax6 = fig.add_subplot(grid[1,2:4])
ax7 = fig.add_subplot(grid[2,2:4])
ax8 = fig.add_subplot(grid[3,2:4])

size = 0.5

colormap = pltconf.set_colormap(df['Attractor'])
fillcolor = 'lightgray'
ax1.scatter(dfpoinc['Cpar'], dfpoinc['x[0]'], c = dfpoinc['Attractor'], cmap = colormap, rasterized = True, s = size, linewidths = 0, marker = 'o', zorder = 2)
ax1.plot(df['Cpar'], df['xMAX[0]'], rasterized = True, color = fillcolor, zorder = 1)
ax1.plot(df['Cpar'], df['xMIN[0]'], rasterized = True, color = fillcolor, zorder = 1)
ax1.fill_between(df['Cpar'], df['xMAX[0]'], df['xMIN[0]'], color = fillcolor, zorder = 0)
ax1.set_xlim(dfpoinc['Cpar'].min(), dfpoinc['Cpar'].max())
ax1.set_ylabel(r'$x$')

ax2.scatter(dfpoinc['Cpar'], dfpoinc['x[1]'], c = dfpoinc['Attractor'], cmap = colormap, rasterized = True, s = size, linewidths = 0, marker = 'o', zorder = 2)
ax2.plot(df['Cpar'], df['xMAX[1]'], rasterized = True, color = fillcolor, zorder = 1)
ax2.plot(df['Cpar'], df['xMIN[1]'], rasterized = True, color = fillcolor, zorder = 1)
ax2.fill_between(df['Cpar'], df['xMAX[1]'], df['xMIN[1]'], color = fillcolor, zorder = 0)
ax2.set_xlim(dfpoinc['Cpar'].min(), dfpoinc['Cpar'].max())
ax2.set_ylabel(r'$\dot{x}$')

ax3.scatter(df['Cpar'], df['LE[0]'], c = df['Attractor'], cmap = colormap, rasterized = True, s = size, linewidths = 0, marker = 'o', zorder = 1)
ax3.set_xlim(df['Cpar'].min(), df['Cpar'].max())
ax3.axhline(0, df['Cpar'].min(), df['Cpar'].max(), color = 'black', linewidth = 0.5, zorder = 0)
ax3.set_ylabel(r'$\lambda_1$')

ax4.scatter(df['Cpar'], df['LE[1]'], c = df['Attractor'], cmap = colormap, rasterized = True, s = size, linewidths = 0, marker = 'o', zorder = 1)
ax4.set_ylabel(r'$\lambda_2$')
ax4.set_xlabel(r'$\Omega$')
ax4.set_xlim(df['Cpar'].min(), df['Cpar'].max())

#========================================================================#
# Show and Save Figure                                                   #
#========================================================================#

if save == True:
    isExist = os.path.exists(savepath)
    if (isExist == False):
        os.makedirs(savepath)
    
    name = "/sample_" + system + "_fbifurc" + ext; name = pltconf.convert_dir(name)
    fig.savefig(savepath + name)
    plt.show()
else:
    plt.show()