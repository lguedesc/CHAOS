#%% -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import os
from matplotlib.ticker import FormatStrFormatter
from src.libs import plotconfig as pltconf

pltconf.plot_params(True, 10, 0.5)

def remainder_dataframe(dataframe, divisor):
    # Find the exact values
    exactvalue = dataframe/divisor
    # Round the exact values to the closest integer
    n = round(exactvalue)
    # Find the remainder
    remainder = dataframe - n*divisor
    return remainder


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

gama = 0.2
mu = 45
size = 0.25

cols = 4
rows = 2

df.loc[df['xMAX[4]'] > np.pi, ['xMAX[4]']] = np.pi
df.loc[df['xMAX[4]'] < -np.pi, ['xMAX[4]']] = np.pi

df.loc[df['xMIN[4]'] > np.pi, ['xMIN[4]']] = -np.pi
df.loc[df['xMIN[4]'] < -np.pi, ['xMIN[4]']] = -np.pi

dfpoinc['x[4]'] = remainder_dataframe(dfpoinc['x[4]'], 2*np.pi)
df['xRMS[4]'] = remainder_dataframe(df['xRMS[4]'], 2*np.pi)


fig, axs = plt.subplots(rows+1, cols, figsize = (x_inches2,y_inches2), dpi = dpi, constrained_layout = True)
fig.set_constrained_layout_pads(hspace=0, wspace=0.1)

n = [0, 1, 2, 3, 4, 5, 6, 7]
i = 0
mycolors = ['lightblue', 'lightblue', 'lightgreen', 'lightgreen', 'lightsalmon', 'lightsalmon', 'orange', 'gold']
rmscolors = ['blue', 'blue', 'darkgreen', 'darkgreen', 'red', 'red', 'orangered', 'orange']
names = [r'$x$', r'$\dot{x}$',r'$z$', r'$\dot{z}$', r'$\phi$', r'$\dot{\phi}$', r'$v$', r'$I$']
for col in range(cols):
    for row in range(rows):
        ax = axs[row, col]
        ax.scatter(dfpoinc['Cpar'], dfpoinc[f'x[{n[i]}]'], rasterized = True, color = "black", s = size, linewidths = 0, marker = '.', zorder = 2)
        ax.plot(df['Cpar'], df[f'xMAX[{n[i]}]'], rasterized = True, color = mycolors[i], lw = 0.5, zorder = 1)
        ax.plot(df['Cpar'], df[f'xMIN[{n[i]}]'], rasterized = True, color = mycolors[i], lw = 0.5, zorder = 1)
        if (i != 4):
            ax.plot(df['Cpar'], df[f'xRMS[{n[i]}]'], rasterized = True, color = rmscolors[i], lw = 0.5, zorder = 3)
        ax.fill_between(df['Cpar'], df[f'xMAX[{n[i]}]'], df[f'xMIN[{n[i]}]'], color = mycolors[i], zorder = 0)
        ax.set_ylabel(names[i])
        ax.set_xlabel(r'$\Omega$')
        ax.set_xlim(dfpoinc['Cpar'].min(), dfpoinc['Cpar'].max())
        ax.set_xticks([0.01, 1, 2])
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.g'))
        i = i + 1

row1 = 2; col1 = 0
axs[row1, col1].plot(df['Cpar'], df[f'ddX_MIN'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs[row1, col1].plot(df['Cpar'], df[f'ddX_MAX'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs[row1, col1].plot(df['Cpar'], df[f'ddX_RMS'], rasterized = True, color = 'blue', lw = 0.5, zorder = 1)
axs[row1, col1].fill_between(df['Cpar'], df[f'ddX_MAX'], df[f'ddX_MIN'], color = 'lightblue', zorder = 0)
axs[row1, col1].set_ylabel(r'$\ddot{x}$')
axs[row1, col1].set_xlabel(r'$\Omega$')
axs[row1, col1].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs[row1, col1].set_xticks([0.01, 1, 2])
axs[row1, col1].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row2 = 2; col2 = 1
axs[row2, col2].plot(df['Cpar'], df[f'ddZ_MIN'], rasterized = True, color = 'lightgreen', lw = 0.5, zorder = 1)
axs[row2, col2].plot(df['Cpar'], df[f'ddZ_MAX'], rasterized = True, color = 'lightgreen', lw = 0.5, zorder = 1)
axs[row2, col2].plot(df['Cpar'], df[f'ddZ_RMS'], rasterized = True, color = 'green', lw = 0.5, zorder = 1)
axs[row2, col2].fill_between(df['Cpar'], df[f'ddZ_MAX'], df[f'ddZ_MIN'], color = 'lightgreen', zorder = 0)
axs[row2, col2].set_ylabel(r'$\ddot{z}$')
axs[row2, col2].set_xlabel(r'$\Omega$')
axs[row2, col2].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs[row2, col2].set_xticks([0.01, 1, 2])
axs[row2, col2].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row3 = 2; col3 = 2
axs[row3, col3].plot(df['Cpar'], df[f'ddPhi_MIN'], rasterized = True, color = 'lightsalmon', lw = 0.5, zorder = 1)
axs[row3, col3].plot(df['Cpar'], df[f'ddPhi_MAX'], rasterized = True, color = 'lightsalmon', lw = 0.5, zorder = 1)
axs[row3, col3].plot(df['Cpar'], df[f'ddPhi_RMS'], rasterized = True, color = 'red', lw = 0.5, zorder = 1)
axs[row3, col3].fill_between(df['Cpar'], df[f'ddPhi_MAX'], df[f'ddPhi_MIN'], color = 'lightsalmon', zorder = 0)
axs[row3, col3].set_ylabel(r'$\ddot{\phi}$')
axs[row3, col3].set_xlabel(r'$\Omega$')
axs[row3, col3].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs[row3, col3].set_xticks([0.01, 1, 2])
axs[row3, col3].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row4 = 2; col4 = 3
axs[row4, col4].plot(df['Cpar'], df[f'PoutPZ_Avg'], rasterized = True, color = 'red', lw = 0.5, zorder = 1, label='PZ')
axs[row4, col4].plot(df['Cpar'], df[f'PoutEM_Avg'], rasterized = True, color = 'brown', lw = 0.5, zorder = 1, label='EM')
axs[row4, col4].fill_between(df['Cpar'], df[f'PoutPZ_Avg'] + df[f'PoutEM_Avg'], 0, color = 'pink', zorder = 0)
axs[row4, col4].set_ylabel(r'$P_{\mathrm{avg}}$')
axs[row4, col4].set_xlabel(r'$\Omega$')
axs[row4, col4].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs[row4, col4].set_xticks([0.01, 1, 2])
axs[row4, col4].xaxis.set_major_formatter(FormatStrFormatter('%.g'))


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

row = 0; col = 1
axs2[row, col].plot(df['Cpar'], df[f'Zcm_MIN'], rasterized = True, color = 'lightgreen', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'Zcm_MAX'], rasterized = True, color = 'lightgreen', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'Zcm_RMS'], rasterized = True, color = 'darkgreen', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'Zcm_MAX'], df[f'Zcm_MIN'], color = 'lightgreen', zorder = 0)
axs2[row, col].set_ylabel(r'$z_{\mathrm{cm}}$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row = 1; col = 0
axs2[row, col].plot(df['Cpar'], df[f'dXcm_MIN'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'dXcm_MAX'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'dXcm_RMS'], rasterized = True, color = 'blue', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'dXcm_MAX'], df[f'dXcm_MIN'], color = 'lightblue', zorder = 0)
axs2[row, col].set_ylabel(r'$x_{\mathrm{cm}}$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row = 1; col = 1
axs2[row, col].plot(df['Cpar'], df[f'dZcm_MIN'], rasterized = True, color = 'lightgreen', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'dZcm_MAX'], rasterized = True, color = 'lightgreen', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'dZcm_RMS'], rasterized = True, color = 'darkgreen', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'dZcm_MAX'], df[f'dZcm_MIN'], color = 'lightgreen', zorder = 0)
axs2[row, col].set_ylabel(r'$z_{\mathrm{cm}}$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row = 2; col = 0
axs2[row, col].plot(df['Cpar'], df[f'ddXcm_MIN'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'ddXcm_MAX'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'ddXcm_RMS'], rasterized = True, color = 'blue', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'ddXcm_MAX'], df[f'ddXcm_MIN'], color = 'lightblue', zorder = 0)
axs2[row, col].set_ylabel(r'$x_{\mathrm{cm}}$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row = 2; col = 1
axs2[row, col].plot(df['Cpar'], df[f'ddZcm_MIN'], rasterized = True, color = 'lightgreen', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'ddZcm_MAX'], rasterized = True, color = 'lightgreen', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'ddZcm_RMS'], rasterized = True, color = 'darkgreen', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'ddZcm_MAX'], df[f'ddZcm_MIN'], color = 'lightgreen', zorder = 0)
axs2[row, col].set_ylabel(r'$z_{\mathrm{cm}}$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

df['Xb_MAX'] = gama*np.sin((np.pi/180)*mu)
df['Xb_MIN'] = -gama*np.sin((np.pi/180)*mu)
df['dXb_MAX'] = gama*df['Cpar']*np.sin((np.pi/180)*mu)
df['dXb_MIN'] = -gama*df['Cpar']*np.sin((np.pi/180)*mu)
df['ddXb_MAX'] = gama*df['Cpar']*df['Cpar']*np.sin((np.pi/180)*mu)
df['ddXb_MIN'] = -gama*df['Cpar']*df['Cpar']*np.sin((np.pi/180)*mu)

row = 0; col = 2
axs2[row, col].plot(df['Cpar'], df[f'Xb_MIN'] , rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'Xb_MAX'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'OVRLL_Xb_RMS'], rasterized = True, color = 'blue', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'Xb_MAX'], df[f'Xb_MIN'], color = 'lightblue', zorder = 0)
axs2[row, col].set_ylabel(r'$x_b$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row = 1; col = 2
axs2[row, col].plot(df['Cpar'], df[f'dXb_MIN'] , rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'dXb_MAX'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'OVRLL_dXb_RMS'], rasterized = True, color = 'blue', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'dXb_MAX'], df[f'dXb_MIN'], color = 'lightblue', zorder = 0)
axs2[row, col].set_ylabel(r'$\dot{x}_b$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row = 2; col = 2
axs2[row, col].plot(df['Cpar'], df[f'ddXb_MIN'] , rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'ddXb_MAX'], rasterized = True, color = 'lightblue', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'OVRLL_ddXb_RMS'], rasterized = True, color = 'blue', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'ddXb_MAX'], df[f'ddXb_MIN'], color = 'lightblue', zorder = 0)
axs2[row, col].set_ylabel(r'$\ddot{x}_b$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

df['Zb_MAX'] = gama*np.cos((np.pi/180)*mu)
df['Zb_MIN'] = -gama*np.cos((np.pi/180)*mu)
df['dZb_MAX'] = gama*df['Cpar']*np.cos((np.pi/180)*mu)
df['dZb_MIN'] = -gama*df['Cpar']*np.cos((np.pi/180)*mu)
df['ddZb_MAX'] = gama*df['Cpar']*df['Cpar']*np.cos((np.pi/180)*mu)
df['ddZb_MIN'] = -gama*df['Cpar']*df['Cpar']*np.cos((np.pi/180)*mu)

row = 0; col = 3
axs2[row, col].plot(df['Cpar'], df[f'Zb_MIN'] , rasterized = True, color = 'lightgreen', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'Zb_MAX'], rasterized = True, color = 'lightgreen', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'OVRLL_Zb_RMS'], rasterized = True, color = 'darkgreen', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'Zb_MAX'], df[f'Zb_MIN'], color = 'lightgreen', zorder = 0)
axs2[row, col].set_ylabel(r'$z_b$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row = 1; col = 3
axs2[row, col].plot(df['Cpar'], df[f'dZb_MIN'] , rasterized = True, color = 'lightgreen', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'dZb_MAX'], rasterized = True, color = 'lightgreen', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'OVRLL_dZb_RMS'], rasterized = True, color = 'darkgreen', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'dZb_MAX'], df[f'dZb_MIN'], color = 'lightgreen', zorder = 0)
axs2[row, col].set_ylabel(r'$\dot{z}_b$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row = 2; col = 3
axs2[row, col].plot(df['Cpar'], df[f'ddZb_MIN'] , rasterized = True, color = 'lightgreen', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'ddZb_MAX'], rasterized = True, color = 'lightgreen', lw = 0.5, zorder = 1)
axs2[row, col].plot(df['Cpar'], df[f'OVRLL_ddZb_RMS'], rasterized = True, color = 'darkgreen', lw = 0.5, zorder = 1)
axs2[row, col].fill_between(df['Cpar'], df[f'ddZb_MAX'], df[f'ddZb_MIN'], color = 'lightgreen', zorder = 0)
axs2[row, col].set_ylabel(r'$\ddot{x}_b$')
axs2[row, col].set_xlabel(r'$\Omega$')
axs2[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs2[row, col].set_xticks([0.01, 1, 2])
axs2[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))


rows = 2; cols = 3
fig3, axs3 = plt.subplots(rows, cols, figsize = (x_inches2,y_inches2), dpi = dpi, constrained_layout = True)
fig3.set_constrained_layout_pads(hspace=0, wspace=0.1)

row = 0; col = 0
axs3[row, col].plot(df['Cpar'], df[f'pos_spin'] , rasterized = True, color = 'darkgreen', lw = 0.5, zorder = 1)
axs3[row, col].set_ylabel(r'$\mathrm{spin}(+)$')
axs3[row, col].set_xlabel(r'$\Omega$')
axs3[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs3[row, col].set_xticks([0.01, 1, 2])
axs3[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row = 0; col = 1
axs3[row, col].plot(df['Cpar'], df[f'neg_spin'] , rasterized = True, color = 'darkred', lw = 0.5, zorder = 1)
axs3[row, col].set_ylabel(r'$\mathrm{spin}(-)$')
axs3[row, col].set_xlabel(r'$\Omega$')
axs3[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs3[row, col].set_xticks([0.01, 1, 2])
axs3[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

row = 0; col = 2
axs3[row, col].plot(df['Cpar'], df[f'tflip'] , rasterized = True, color = 'blue', lw = 0.5, zorder = 1)
axs3[row, col].set_ylabel(r'$t_{\mathrm{flip}}$')
axs3[row, col].set_xlabel(r'$\Omega$')
axs3[row, col].set_xlim(df['Cpar'].min(), df['Cpar'].max())
axs3[row, col].set_xticks([0.01, 1, 2])
axs3[row, col].xaxis.set_major_formatter(FormatStrFormatter('%.g'))

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
