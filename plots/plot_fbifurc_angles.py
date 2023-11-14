#%% -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
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

varphi_pz = 0.05
varphi_em = 5.0

#system = "lin_oscillator_2DoF"
#system = "lin_2DoF_EH"
#system = "duffing_2DoF_EH"
#system = "bistable_EH"
#system = "duffing"
#system = "pend_oscillator_EH"
system = "multidirect_hybrid_EH"
#system = "pendulum_EMEH"
#system = "pendulum_EMEH_dimensional"
ext = ".pdf"

filenum = 0
simulation = "fbifurc"
dim = 8
angles = True
angles_indexes = { 4 }

df, dfpoinc = pltconf.read_CHAOS_data(system, filenum, simulation)
#=======================================================================#
# Figure Parameters                                                     #
#=======================================================================#
cols = 1
rows = dim

each_row_size_ratio = 0.15
figsize = pltconf.figsize_in_cm(15, each_row_size_ratio*rows*15)
figsize2 = pltconf.figsize_in_cm(15, each_row_size_ratio*len(angles_indexes)*15)
figsize3 = pltconf.figsize_in_cm(15, each_row_size_ratio*6*15)
dpi = pltconf.set_fig_quality(save = save, base_dpi = 100)

gama = 0.1
size = 5

fig, axs = pltconf.makefig_and_axs(figsize, rows, cols, dpi, hspace = 0.1, wspace = 0.1)

n = [0, 1, 2, 3, 4, 5, 6, 7]
i = 0
color = 'lightgray'
rmscolor = 'red'
names = [r'$x$', r'$\dot{x}$', r'$z$', r'$\dot{z}$', r'$\phi$', r'$\dot{\phi}$', r'$v$', r'$I$']
#names = [r'$\phi$', r'$\dot{\phi}$', r'$I$']
xticks = [df['Cpar'].min(), 1, df['Cpar'].max()]


colormap = pltconf.set_colormap(df['Attractor'])

for col in range(cols):
    for row in range(rows):
        if i < dim:
            if cols == 1:
                ax = axs[row]
            else:
                ax = axs[row, col]
            ax.scatter(dfpoinc['Cpar'], dfpoinc[f'x[{n[i]}]'], c = dfpoinc['Attractor'], cmap = colormap, rasterized = True, s = size, linewidths = 0, marker = '.', zorder = 2)
            ax.plot(df['Cpar'], df[f'xMAX[{n[i]}]'], rasterized = True, color = color, lw = 0.5, zorder = 1)
            ax.plot(df['Cpar'], df[f'xMIN[{n[i]}]'], rasterized = True, color = color, lw = 0.5, zorder = 1)
            #ax.plot(df['Cpar'], df[f'xRMS[{n[i]}]'], rasterized = True, color = rmscolor, lw = 0.5, zorder = 3)
            ax.fill_between(df['Cpar'], df[f'xMAX[{n[i]}]'], df[f'xMIN[{n[i]}]'], color = color, zorder = 0)
            ax.set_ylabel(names[i])
            ax.set_xlabel(r'$\Omega$')
            ax.set_xlim(dfpoinc['Cpar'].min(), dfpoinc['Cpar'].max())
            ax.set_xticks(xticks)
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.g'))
            i = i + 1
        else:
            pass

plt.show(block = False)

fig2, axs2 = pltconf.makefig_and_axs(figsize2, len(angles_indexes), cols, dpi, hspace = 0.1, wspace = 0.1)
# Check if axs2 is a single subplot
if isinstance(axs2, Axes):
    axs2 = [axs2]  # Convert single subplot to a list with one element
    
for ax, i in zip(axs2, angles_indexes):
        ax.scatter(dfpoinc['Cpar'], dfpoinc[f'x[{i}]_norm'], rasterized = True, c = dfpoinc['Attractor'], cmap = colormap, s = size, linewidths = 0, marker = '.', zorder = 2)
        ax.plot(df['Cpar'], df[f'xMAX[{i}]_norm'], rasterized = True, color = color, lw = 0.5, zorder = 1)
        ax.plot(df['Cpar'], df[f'xMIN[{i}]_norm'], rasterized = True, color = color, lw = 0.5, zorder = 1)
        #ax.plot(df['Cpar'], df[f'xRMS[{i}]'], rasterized = True, color = rmscolors[i], lw = 0.5, zorder = 3)
        ax.fill_between(df['Cpar'], df[f'xMAX[{i}]_norm'], df[f'xMIN[{i}]_norm'], color = color, zorder = 0)
        ax.set_ylabel(names[i])
        ax.set_xlabel(r'$\Omega$')
        ax.set_xlim(dfpoinc['Cpar'].min(), dfpoinc['Cpar'].max())
        #ax.set_xticks(xticks)
        #ax.xaxis.set_major_formatter(FormatStrFormatter('%.g'))
        ax.set_yticks([-np.pi, 0, np.pi])
        ax.set_yticklabels([r"$-\pi$", 0, r"$\pi$"])

fig3, axs3 = pltconf.makefig_and_axs(figsize2, 6, cols, dpi, hspace = 0.1, wspace = 0.1)

PZ_clr = "red"
EM_clr = "darkorange"
Total_clr = "purple"

scale = 1e3
df["Pout_PZ"] = df["xRMS[6]"].multiply(df["xRMS[6]"].multiply(varphi_pz*scale))
df["Pout_EM"] = df["xRMS[7]"].multiply(df["xRMS[7]"].multiply(varphi_em*scale))

axs3[1].plot(df["Cpar"], df["Pout_PZ"], color = PZ_clr)
axs3[1].fill_between(df['Cpar'], df["Pout_PZ"], 0, color=PZ_clr, alpha=.5)
axs3[1].set_ylabel(r"$\bar{P}_{\mathrm{out}}^{(pz)}$")
axs3[1].set_xlabel(r"$\Omega$")

axs3[2].plot(df["Cpar"], df["Pout_EM"], color = EM_clr)
axs3[2].fill_between(df['Cpar'], df["Pout_EM"], 0, color=EM_clr, alpha=.5)
axs3[2].set_ylabel(r"$\bar{P}_{\mathrm{out}}^{(em)}$")
axs3[2].set_xlabel(r"$\Omega$")


TotalPout = df["Pout_PZ"] + df["Pout_EM"]
axs3[3].plot(df["Cpar"], TotalPout, color = Total_clr)
axs3[3].fill_between(df['Cpar'], TotalPout, 0, color=Total_clr, alpha=.5)
axs3[3].set_ylabel(r"$\bar{P}_{\mathrm{out}}$")
axs3[3].set_xlabel(r"$\Omega$")

ax.plot()

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
