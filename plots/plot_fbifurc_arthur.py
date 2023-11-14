#%% -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
from libs import plotconfig as pltconf
from matplotlib.axes import Axes

pltconf.plot_params(True, 10, 0.2)

def remainder_dataframe(dataframe, divisor):
    # Find the exact values
    exactvalue = dataframe/divisor
    # Round the exact values to the closest integer
    n = round(exactvalue)
    # Find the remainder
    remainder = dataframe - n*divisor
    return remainder

save = False

#system = "lin_oscillator_2DoF"
#system = "lin_2DoF_EH"
#system = "duffing_2DoF_EH"
#system = "bistable_EH"
#system = "duffing"
#system = "pend_oscillator_EH"
#system = "multidirect_hybrid_EH"
#system = "pendulum_EMEH"
system = "pendulum_EMEH_dimensional"
ext = ".pdf"

filenum = 2
simulation = "fbifurc"
dim = 3
angles = True
angles_indexes = { 0 }


df = pd.read_csv("/Users/luaguedescosta/Downloads/Output_f_Lua-2.csv", sep=" ")
#=======================================================================#
# Figure Parameters                                                     #
#=======================================================================#
cols = 1
rows = dim

figsize = pltconf.figsize_in_cm(15, 0.2*rows*15)
figsize2 = pltconf.figsize_in_cm(15, 0.2*len(angles_indexes)*15)
figsize3 = pltconf.figsize_in_cm(15, 0.2*15)
dpi = pltconf.set_fig_quality(save = save, base_dpi = 100)

gama = 0.1
size = 5

fig, axs = pltconf.makefig_and_axs(figsize, rows, cols, dpi, hspace = 0.1, wspace = 0.1)

n = [0, 1, 2, 3, 4, 5, 6, 7]
i = 0
color = 'lightgray'
rmscolor = 'red'
#names = [r'$x$', r'$\dot{x}$', r'$z$', r'$\dot{z}$', r'$\phi$', r'$\dot{\phi}$', r'$v$', r'$I$']
names = [r'$\phi$', r'$\dot{\phi}$', r'$I$']
xticks = [df['om[-]'].min(), 1, df['om[-]'].max()]

for col in range(cols):
    for row in range(rows):
        if i < dim:
            if cols == 1:
                ax = axs[row]
            else:
                ax = axs[row, col]
            ax.plot(df['om[-]'], df[f'y_max[{n[i]}]'], rasterized = True, color = color, lw = 0.5, zorder = 1)
            ax.plot(df['om[-]'], df[f'y_min[{n[i]}]'], rasterized = True, color = color, lw = 0.5, zorder = 1)
            #ax.plot(df['Cpar'], df[f'xRMS[{n[i]}]'], rasterized = True, color = rmscolor, lw = 0.5, zorder = 3)
            ax.fill_between(df['om[-]'], df[f'y_max[{n[i]}]'], df[f'y_min[{n[i]}]'], color = color, zorder = 0)
            ax.set_ylabel(names[i])
            ax.set_xlabel(r'$\Omega$')
            ax.set_xlim(df['om[-]'].min(), df['om[-]'].max())
            ax.set_xticks(xticks)
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.g'))
            i = i + 1
        else:
            pass

plt.show()
#========================================================================#
# Show and Save Figure                                                   #
#========================================================================#

if save == True:
    pltconf.save_CHAOS_data(fig, system, simulation, ext)
    #plt.show()
#else:
#    fig.show()
#    plt.show()

    #plt.show(fig)
    #plt.show(fig2)
    #plt.show(fig3)


# %%
