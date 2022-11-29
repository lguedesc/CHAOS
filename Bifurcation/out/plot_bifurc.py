#%% -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import os
from matplotlib.ticker import FormatStrFormatter
from src.libs import plotconfig as pltconf

pltconf.plot_params(True, 10, 0.5)

save = False

system = "pend_oscillator_EH"
ext = ".pdf"

num = 1

if num > 0:
    readpath = "Bifurcation/out/" + system + f"_bifurc({num}).csv"; readpath = pltconf.convert_dir(readpath)
    readpathpoinc = "Bifurcation/out/" + system + f"_bifurc_poinc({num}).csv"; readpathpoinc = pltconf.convert_dir(readpathpoinc)
    savepath = "Bifurcation/figs"; savepath = pltconf.convert_dir(savepath)
else:    
    readpath = "Bifurcation/out/" + system + f"_bifurc.csv"; readpath = pltconf.convert_dir(readpath)
    readpathpoinc = "Bifurcation/out/" + system + f"_bifurc_poinc.csv"; readpathpoinc = pltconf.convert_dir(readpathpoinc)
    savepath = "Bifurcation/figs"; savepath = pltconf.convert_dir(savepath)
    
dfpoinc = pd.read_csv(readpathpoinc, delimiter = " ")
df = pd.read_csv(readpath, delimiter = " ")

#=======================================================================#
# Figure Parameters                                                     #
#=======================================================================#
x_inches2 = 200*(1/25.4)     # [mm]*constant
y_inches2 = x_inches2*(0.4)

if save == True:
    dpi = 300
else:
    dpi = 200

size = 0.25

cols = 4
rows = 2

fig, axs = plt.subplots(rows, cols, figsize = (x_inches2,y_inches2), dpi = dpi, constrained_layout = True)
fig.set_constrained_layout_pads(hspace=0, wspace=0.1)

n = [0, 1, 2, 3, 4, 5, 6, 7]
i = 0
mycolors = ['lightblue', 'lightblue', 'lightgreen', 'lightgreen', 'lightsalmon', 'lightsalmon', 'orange', 'gold']
rmscolors = ['blue', 'blue', 'darkgreen', 'darkgreen', 'red', 'red', 'orangered', 'orange']
for col in range(cols):
    for row in range(rows):
        ax = axs[row, col]
        ax.scatter(dfpoinc['Cpar'], dfpoinc[f'x[{n[i]}]'], rasterized = True, color = "black", s = size, linewidths = 0, marker = '.', zorder = 2)
        ax.plot(df['Cpar'], df[f'xMAX[{n[i]}]'], rasterized = True, color = mycolors[i], lw = 0.5, zorder = 1)
        ax.plot(df['Cpar'], df[f'xMIN[{n[i]}]'], rasterized = True, color = mycolors[i], lw = 0.5, zorder = 1)
        ax.plot(df['Cpar'], df[f'xRMS[{n[i]}]'], rasterized = True, color = rmscolors[i], lw = 0.5, zorder = 3)
        ax.fill_between(df['Cpar'], df[f'xMAX[{n[i]}]'], df[f'xMIN[{n[i]}]'], color = mycolors[i], zorder = 0)
        ax.set_ylabel(f'$x_{n[i]}$')
        ax.set_xlabel(r'$\Omega$')
        ax.set_xlim(dfpoinc['Cpar'].min(), dfpoinc['Cpar'].max())
        ax.set_xticks([0.01, 1, 2])
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.g'))
        i = i + 1

fig2, axs2 = plt.subplots(rows, cols, figsize = (x_inches2,y_inches2), dpi = dpi, constrained_layout = True)
fig.set_constrained_layout_pads(hspace=0, wspace=0.1)

axs2[0, 0].plot(df['Cpar'], df[f'ddX_MIN'], rasterized = True, color = 'red', lw = 0.5, zorder = 1)
axs2[0, 0].plot(df['Cpar'], df[f'ddX_MAX'], rasterized = True, color = 'red', lw = 0.5, zorder = 1)
axs2[0, 0].fill_between(df['Cpar'], df[f'ddX_MAX'], df[f'ddX_MIN'], color = 'red', zorder = 0)


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
