#%% -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import os
from src.libs import plotconfig as pltconf

pltconf.plot_params(True, 10, 0.5)

save = False
system = "duffing"
#system2 = "duffing_cldyn"
ext = ".pdf"

readpath = "LyapunovExp/out/" + system + "_lyap(1).csv"; readpath = pltconf.convert_dir(readpath)
#readpath2 = "LyapunovExp/out/" + system2 + "_lyap(7).csv"; readpath = pltconf.convert_dir(readpath)
savepath = "LyapunovExp/figs"; savepath = pltconf.convert_dir(savepath)



df = pd.read_csv(readpath, delimiter = " ")
#df2 = pd.read_csv(readpath2, delimiter = " ")

#=======================================================================#
# Figure Parameters                                                     #
#=======================================================================#
x_inches2 = 150*(1/25.4)     # [mm]*constant
y_inches2 = x_inches2*(0.4)

if save == True:
    dpi = 2000
else:
    dpi = 200

fig = plt.figure(1, figsize = (x_inches2,y_inches2), dpi = dpi, constrained_layout = True)
fig.set_constrained_layout_pads(hspace=0, wspace=0.1)

ax = fig.add_subplot(1,1,1)

ax.plot(df['Time'], df['LE[0]'], rasterized = True, color = "red", linewidth = 1, zorder = 1, label = "$\lambda_1 \mathrm{(TanMap)}$")
ax.plot(df['Time'], df['LE[1]'], rasterized = True, color = "blue", linewidth = 1, zorder = 1, label = "$\lambda_2 \mathrm{(TanMap)}$")
ax.plot(df['Time'], df['sLE[0]'], rasterized = True, color = "black", linewidth = 1, zorder = 1, label = "$\lambda_1 \mathrm{(TanMap)}$")
ax.plot(df['Time'], df['sLE[1]'], rasterized = True, color = "orange", linewidth = 1, zorder = 1, label = "$\lambda_2 \mathrm{(TanMap)}$")
#ax.plot(df2['Time'], df2['LE[0]'], rasterized = True, color = "black", linestyle = "dashed", linewidth = 0.5, zorder = 2, label = "$\lambda_1 \mathrm{(ClDyn)}$")
#ax.plot(df2['Time'], df2['LE[1]'], rasterized = True, color = "orange", linestyle = "dashed", linewidth = 0.5, zorder = 2, label = "$\lambda_2 \mathrm{(ClDyn)}$")
ax.hlines(0, df['Time'].min(), df['Time'].max(), color = "black", linewidth = 1)
ax.set_ylabel(r'$\lambda$')
ax.set_xlabel(r'$\tau$')
ax.set_xlim(df['Time'].min(), df['Time'].max())
ax.legend(loc = "best")
#========================================================================#
# Show and Save Figure                                                   #
#========================================================================#

if save == True:
    isExist = os.path.exists(savepath)
    if (isExist == False):
        os.makedirs(savepath)
    name = "/sample_" + system + "_lyap" + ext; name = pltconf.convert_dir(name)
    #name = "/sample_" + system + system2 + "_lyap" + ext; name = pltconf.convert_dir(name)
    fig.savefig(savepath + name)
    plt.show()
else:
    plt.show()

