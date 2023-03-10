#%% -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import FormatStrFormatter
from libs import plotconfig as pltconf

pltconf.plot_params(True, 10, 0.2)

def remainder_dataframe(dataframe, divisor):
    # Find the exact values
    exactvalue = dataframe/divisor
    # Round the exact values to the closest integer
    n = round(exactvalue)
    # Find the remainder
    remainder = dataframe - n*divisor
    return remainder

save = True

#system = "lin_oscillator_2DoF"
#system = "lin_2DoF_EH"
#system = "duffing_2DoF_EH"
system = "bistable_EH"
ext = ".pdf"

filenum = 0
simulation = "bifurc"
dim = 3


df, dfpoinc = pltconf.read_CHAOS_data(system, filenum, simulation)
#=======================================================================#
# Figure Parameters                                                     #
#=======================================================================#
figsize = pltconf.figsize_in_cm(15, 0.6*15)
dpi = pltconf.set_fig_quality(save = save, base_dpi = 100)

gama = 0.1
size = 0.25

cols = 1
rows = 3

fig, axs = pltconf.makefig_and_axs(figsize, rows, cols, dpi, hspace = 0.1, wspace = 0.1)

n = [0, 1, 2, 3, 4, 5]
i = 0
mycolors = ['lightblue', 'lightgreen', 'lightsalmon']
rmscolors = ['blue', 'darkgreen', 'red']
names = [r'$x$', r'$\dot{x}$', r'$v$']
xticks = [0.01, 1, 2, 3]

for col in range(cols):
    for row in range(rows):
        if i < dim:
            if cols == 1:
                ax = axs[row]
            else:
                ax = axs[row, col]
            ax.scatter(dfpoinc['Cpar'], dfpoinc[f'x[{n[i]}]'], rasterized = True, color = "black", s = size, linewidths = 0, marker = '.', zorder = 2)
            ax.plot(df['Cpar'], df[f'xMAX[{n[i]}]'], rasterized = True, color = mycolors[i], lw = 0.5, zorder = 1)
            ax.plot(df['Cpar'], df[f'xMIN[{n[i]}]'], rasterized = True, color = mycolors[i], lw = 0.5, zorder = 1)
            #ax.plot(df['Cpar'], df[f'xRMS[{n[i]}]'], rasterized = True, color = rmscolors[i], lw = 0.5, zorder = 3)
            ax.fill_between(df['Cpar'], df[f'xMAX[{n[i]}]'], df[f'xMIN[{n[i]}]'], color = mycolors[i], zorder = 0)
            ax.set_ylabel(names[i])
            ax.set_xlabel(r'$\Omega$')
            ax.set_xlim(dfpoinc['Cpar'].min(), dfpoinc['Cpar'].max())
            ax.set_xticks(xticks)
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.g'))
            i = i + 1
        else:
            pass
'''
Poutname = 'TotalPout'
row1 = 2; col1 = 2
#axs[row1, col1].plot(df['Cpar'], df[f'xMIN[4]'] + df[f'xMIN[5]'], rasterized = True, color = 'lightsalmon', lw = 0.5, zorder = 1)
#xs[row1, col1].plot(df['Cpar'], df[f'xMAX[4]'] + df[f'xMAX[5]'], rasterized = True, color = 'lightsalmon', lw = 0.5, zorder = 1)
#axs[row1, col1].plot(df['Cpar'], df[f'xRMS[4]'] + df[f'xRMS[5]'], rasterized = True, color = 'red', lw = 0.5, zorder = 3)
#axs[row1, col1].fill_between(df['Cpar'], df[f'xMIN[4]'] + df[f'xMIN[5]'], df[f'xMAX[4]'] + df[f'xMAX[5]'], color = 'lightsalmon', zorder = 0)
axs[row1, col1].plot(df['Cpar'], df[Poutname], rasterized = True, color = 'lightsalmon', lw = 0.5, zorder = 1)
#axs[row1, col1].plot(df['Cpar'], df[f'Pout1'], rasterized = True, color = 'red', lw = 0.5, zorder = 2)
#axs[row1, col1].plot(df['Cpar'], df[f'Pout2'], rasterized = True, color = 'blue', lw = 0.5, zorder = 2)
axs[row1, col1].fill_between(df['Cpar'], 0, df[Poutname], color = 'lightsalmon', zorder = 0)
axs[row1, col1].set_ylabel(r'$P_{avg}$')
axs[row1, col1].set_xlabel(r'$\Omega$')
axs[row1, col1].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs[row1, col1].set_xticks(xticks)
axs[row1, col1].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row1 = 2; col1 = 0
axs[row1, col1].plot(df['Cpar'], df[f'ddX1_MIN'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs[row1, col1].plot(df['Cpar'], df[f'ddX1_MAX'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs[row1, col1].plot(df['Cpar'], df[f'ddX1_RMS'], rasterized = True, color = 'blue', lw = 0.5, zorder = 1)
axs[row1, col1].fill_between(df['Cpar'], df[f'ddX1_MAX'], df[f'ddX1_MIN'], color = 'lightblue', zorder = 0)
axs[row1, col1].set_ylabel(r'$\ddot{x}_1$')
axs[row1, col1].set_xlabel(r'$\Omega$')
axs[row1, col1].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs[row1, col1].set_xticks([0.01, 1, 2, 3, 4, 5])
axs[row1, col1].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row2 = 2; col2 = 1
axs[row2, col2].plot(df['Cpar'], df[f'ddX2_MIN'], rasterized = True, color = 'lightgreen', lw = 0.5, zorder = 1)
axs[row2, col2].plot(df['Cpar'], df[f'ddX2_MAX'], rasterized = True, color = 'lightgreen', lw = 0.5, zorder = 1)
axs[row2, col2].plot(df['Cpar'], df[f'ddX2_RMS'], rasterized = True, color = 'green', lw = 0.5, zorder = 1)
axs[row2, col2].fill_between(df['Cpar'], df[f'ddX2_MAX'], df[f'ddX2_MIN'], color = 'lightgreen', zorder = 0)
axs[row2, col2].set_ylabel(r'$\ddot{x}_2$')
axs[row2, col2].set_xlabel(r'$\Omega$')
axs[row2, col2].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs[row2, col2].set_xticks([0.01, 1, 2, 3, 4, 5])
axs[row2, col2].xaxis.set_major_formatter(FormatStrFormatter('%.g'))


fig2, axs2 = plt.subplots(rows+1, cols, figsize = (x_inches2,y_inches2), dpi = dpi, constrained_layout = True)
fig2.set_constrained_layout_pads(hspace=0, wspace=0.1)

row = 0; col = 0
axs2[row, col].plot(df['Cpar'], df[f'Xcm_MIN'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'Xcm_MAX'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'Xcm_RMS'], rasterized = True, color = 'blue', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'Xcm_MAX'], df[f'Xcm_MIN'], color = 'lightblue', zorder = 0)
axs2[row, col].set_ylabel(r'$x_{\mathrm{cm}}$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row = 1; col = 0
axs2[row, col].plot(df['Cpar'], df[f'dXcm_MIN'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'dXcm_MAX'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'dXcm_RMS'], rasterized = True, color = 'blue', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'dXcm_MAX'], df[f'dXcm_MIN'], color = 'lightblue', zorder = 0)
axs2[row, col].set_ylabel(r'$\dot{x}_{\mathrm{cm}}$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row = 2; col = 0
axs2[row, col].plot(df['Cpar'], df[f'ddXcm_MIN'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'ddXcm_MAX'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'ddXcm_RMS'], rasterized = True, color = 'blue', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'ddXcm_MAX'], df[f'ddXcm_MIN'], color = 'lightblue', zorder = 0)
axs2[row, col].set_ylabel(r'$\ddot{x}_{\mathrm{cm}}$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

df['Xb_MAX'] = gama
df['Xb_MIN'] = -gama
df['dXb_MAX'] = gama*df['Cpar']
df['dXb_MIN'] = -gama*df['Cpar']
df['ddXb_MAX'] = gama*df['Cpar']*df['Cpar']
df['ddXb_MIN'] = -gama*df['Cpar']*df['Cpar']

row = 0; col = 1
axs2[row, col].plot(df['Cpar'], df[f'Xrel_MIN'] , rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'Xrel_MAX'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'Xrel_RMS'], rasterized = True, color = 'blue', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'Xrel_MAX'], df[f'Xrel_MIN'], color = 'lightblue', zorder = 0)
axs2[row, col].set_ylabel(r'$x_{rel}$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2, 3, 4, 5])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row = 1; col = 1
axs2[row, col].plot(df['Cpar'], df[f'dXrel_MIN'] , rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'dXrel_MAX'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'dXrel_RMS'], rasterized = True, color = 'blue', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'dXrel_MAX'], df[f'dXrel_MIN'], color = 'lightblue', zorder = 0)
axs2[row, col].set_ylabel(r'$x_{rel}$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2, 3, 4, 5])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row = 2; col = 1
axs2[row, col].plot(df['Cpar'], df[f'ddXrel_MIN'] , rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'ddXrel_MAX'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'ddXrel_RMS'], rasterized = True, color = 'blue', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'ddXrel_MAX'], df[f'ddXrel_MIN'], color = 'lightblue', zorder = 0)
axs2[row, col].set_ylabel(r'$x_{rel}$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2, 3, 4, 5])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row = 0; col = 2
axs2[row, col].plot(df['Cpar'], df[f'Xb_MIN'] , rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'Xb_MAX'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'OVRL_Xb_RMS'], rasterized = True, color = 'blue', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'Xb_MAX'], df[f'Xb_MIN'], color = 'lightblue', zorder = 0)
axs2[row, col].set_ylabel(r'$x_b$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2, 3, 4, 5])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row = 1; col = 2
axs2[row, col].plot(df['Cpar'], df[f'dXb_MIN'] , rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'dXb_MAX'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'OVRL_dXb_RMS'], rasterized = True, color = 'blue', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'dXb_MAX'], df[f'dXb_MIN'], color = 'lightblue', zorder = 0)
axs2[row, col].set_ylabel(r'$\dot{x}_b$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2, 3, 4, 5])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row = 2; col = 2
axs2[row, col].plot(df['Cpar'], df[f'ddXb_MIN'] , rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'ddXb_MAX'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'OVRL_ddXb_RMS'], rasterized = True, color = 'blue', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'ddXb_MAX'], df[f'ddXb_MIN'], color = 'lightblue', zorder = 0)
axs2[row, col].set_ylabel(r'$\ddot{x}_b$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))
'''
#========================================================================#
# Show and Save Figure                                                   #
#========================================================================#

if save == True:
    pltconf.save_CHAOS_data(fig, system, simulation, ext)
    #plt.show()
else:
    plt.show()


# %%
