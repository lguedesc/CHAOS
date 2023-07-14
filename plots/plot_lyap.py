#%% -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from libs import plotconfig as pltconf

pltconf.plot_params(True, 10, 0.2)

save = False
system = "duffing"
ext = ".pdf"
filenum = 3
simulation = "lyap"

df = pltconf.read_CHAOS_data(system, filenum, simulation)
#=======================================================================#
# Figure Parameters                                                     #
#=======================================================================#
figsize = pltconf.figsize_in_cm(15, 0.4*15)
dpi = pltconf.set_fig_quality(save = save, base_dpi = 100)

rows = 1; cols = 1
fig, axs = pltconf.makefig_and_axs(figsize, rows, cols, dpi)

axs.plot(df['Time'], df['LE[0]'], rasterized = True, color = "red", linewidth = 1, zorder = 1, label = "$\lambda_1 \mathrm{(TanMap)}$")
axs.plot(df['Time'], df['LE[1]'], rasterized = True, color = "blue", linewidth = 1, zorder = 1, label = "$\lambda_2 \mathrm{(TanMap)}$")
axs.hlines(0, df['Time'].min(), df['Time'].max(), color = "black", linewidth = 1)
axs.set_ylabel(r'$\lambda$')
axs.set_xlabel(r'$\tau$')
axs.set_xlim(df['Time'].min(), df['Time'].max())
axs.legend(loc = "best")
#========================================================================#
# Show and Save Figure                                                   #
#========================================================================#

if save == True:
    pltconf.save_CHAOS_data(fig, system, simulation, ext)
    plt.show()
else:
    plt.show()

