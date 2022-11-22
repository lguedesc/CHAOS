#%% -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import os
from src.libs import plotconfig as pltconf

pltconf.plot_params(True, 10, 0.5)

save = False

system = "bistable_EH"
ext = ".pdf"

num = 9

readpath = "Bifurcation/out/" + system + f"_bifurc({num}).csv"; readpath = pltconf.convert_dir(readpath)
readpathpoinc = "Bifurcation/out/" + system + f"_bifurc_poinc({num}).csv"; readpathpoinc = pltconf.convert_dir(readpathpoinc)
savepath = "Bifurcation/figs"; savepath = pltconf.convert_dir(savepath)
        
dfpoinc = pd.read_csv(readpathpoinc, delimiter = " ")
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
ax3 = fig.add_subplot(grid[0,2:4])
ax4 = fig.add_subplot(grid[1,2:4])

size = 0.25

ax1.scatter(dfpoinc['Cpar'], dfpoinc['x[0]'], rasterized = True, color = "black", s = size, linewidths = 0, marker = 'o', zorder = 2)
ax1.plot(df['Cpar'], df['xMAX[0]'], rasterized = True, color = 'cyan', lw = 0.5, zorder = 1)
ax1.plot(df['Cpar'], df['xMIN[0]'], rasterized = True, color = 'cyan', lw = 0.5, zorder = 1)
ax1.fill_between(df['Cpar'], df['xMAX[0]'], df['xMIN[0]'], color = "cyan", zorder = 0)
ax1.set_ylabel(r'$x$')
ax1.set_xlabel(r'$\Omega$')
ax1.set_xlim(dfpoinc['Cpar'].min(), dfpoinc['Cpar'].max())

ax2.scatter(dfpoinc['Cpar'], dfpoinc['x[2]'], rasterized = True, color = "black", s = size, linewidths = 0, zorder = 2)
ax2.plot(df['Cpar'], df['xMAX[1]'], rasterized = True, color = 'cyan', lw = 0.5, zorder = 1)
ax2.plot(df['Cpar'], df['xMIN[1]'], rasterized = True, color = 'cyan', lw = 0.5, zorder = 1)
ax2.fill_between(df['Cpar'], df['xMAX[1]'], df['xMIN[1]'], color = "cyan", zorder = 0)
ax2.set_ylabel(r'$\dot{x}$')
ax2.set_xlabel(r'$\Omega$')
ax2.set_xlim(dfpoinc['Cpar'].min(), dfpoinc['Cpar'].max())

ax3.scatter(dfpoinc['Cpar'], dfpoinc['x[2]'], rasterized = True, color = "black", s = size, linewidths = 0, zorder = 2)
ax3.plot(df['Cpar'], df['xMAX[2]'], rasterized = True, color = 'orangered', lw = 0.5, zorder = 1)
ax3.plot(df['Cpar'], df['xMIN[2]'], rasterized = True, color = 'orangered', lw = 0.5, zorder = 1)
ax3.fill_between(df['Cpar'], df['xMAX[2]'], df['xMIN[2]'], color = "lightsalmon", zorder = 0)
ax3.set_ylabel(r'$\nu$')
ax3.set_xlabel(r'$\Omega$')
ax3.set_xlim(dfpoinc['Cpar'].min(), dfpoinc['Cpar'].max())

ax4.plot(df['Cpar'], df['EffAvg'], rasterized = True, color = 'blue', lw = 0.5, zorder = 1, label = "RMS")
ax4.hlines(0, df['Cpar'].min(), df['Cpar'].max(), lw = 0.2, color = 'black')
ax4.set_ylabel(r'customvalue')
ax4.set_xlabel(r'$\Omega$')
ax4.set_xlim(df['Cpar'].min(), df['Cpar'].max())
#ax4.legend()

#========================================================================#
# Show and Save Figure                                                   #
#========================================================================#

if save == True:
    isExist = os.path.exists(savepath)
    if (isExist == False):
        os.makedirs(savepath)
    
    name = "/sample_" + system + "_bifurc" + ext; name = pltconf.convert_dir(name)
    fig.savefig(savepath + name)
    plt.show()
else:
    plt.show()


# %%
