# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt    
from matplotlib.ticker import FormatStrFormatter
from libs import plotconfig as pltconf

pltconf.plot_params(True, 9, 0.2)

# =========================================================================== #
# Read Data
# =========================================================================== #
maxper = 6
save = False
system = "pend_oscillator_EH"
ext = ".pdf"

#OMEGA = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 
#         1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00 ]

OMEGA = [1.75]
GAMMA = [0.1]

#path = "/Users/luaguedescosta/Desktop/CHAOS/data/DynDiagram/out/"
path = "data/FDynDiagram/out/"
savepath = f"{path}figs/"
savepath = pltconf.convert_dir(savepath)

for Omega in OMEGA:
    for gamma in GAMMA:
        savename = f"{system}_fdyndiag_O={Omega:.2f}_g={gamma:.1f}"
        file     = f"{system}_fdyndiag_O={Omega:.2f}_g={gamma:.1f}.csv"

        readpath = f"{path}{file}"; readpath = pltconf.convert_dir(readpath)
        
        columnames = ['CparX', 'CparY', 'Attractor', 'PoutPZ_Avg', 'PoutEM_Avg',  'xRMS[0]', 'xRMS[1]', 'xRMS[2]', 'xRMS[3]', 'xRMS[4]', 'xRMS[5]',  'xRMS[6]', 'xRMS[7]', 'xMAX[4]_remainder', 'xMIN[4]_remainder', 'OverallxMAX[4]_remainder', 'OverallxMIN[4]_remainder']
        raw_data = pd.read_csv(readpath, delimiter = " ", usecols = columnames)
        raw_data['TotalPout'] = raw_data['PoutPZ_Avg'] + raw_data['PoutEM_Avg']
        plotcols = ['Attractor', 'PoutPZ_Avg', 'PoutEM_Avg', 'TotalPout',  'xRMS[0]', 'xRMS[1]', 'xRMS[2]', 'xRMS[3]', 'xRMS[4]', 'xRMS[5]',  'xRMS[6]', 'xRMS[7]', 'xMAX[4]_remainder', 'xMIN[4]_remainder', 'OverallxMAX[4]_remainder', 'OverallxMIN[4]_remainder']
        # =========================================================================== #
        # Define figure parameters and Create Figure                         
        # =========================================================================== #
        rows = 4; cols = 4
        figsize = pltconf.figsize_in_cm(20, 0.15*rows*20)
        dpi = pltconf.set_fig_quality(save = save, base_dpi = 200)
        fig, axs = pltconf.makefig_and_axs(figsize, rows, cols, dpi, hspace = 0.1, wspace = 0.1)
        # =========================================================================== #
        # Process and plot data  
        # =========================================================================== #
        for name, ax in zip(plotcols, axs.flat):
            x, y, z = pltconf.process_data(raw_data, 'CparY', 'CparX', name)
            if name == 'Attractor':
                pltconf.plot_attractor_map(fig, ax, x, y, z, maxper, mode = 'lyap')        
            else:
                if (name == 'PoutPZ_Avg') or (name == 'PoutEM_Avg') or (name == 'TotalPout'):
                    z = z*1e3
                plot, cbar = pltconf.plot_rainbow_map(fig, ax, x, y, z)        

            # Customize titles and labels                        
            ax.set_title(f"{name}", fontsize = 8, loc= 'left')
            ax.set_xticks([x.min(), x.max()])
            ax.set_yticks([y.min(), y.max()])
            ax.set_xlabel(r"$\Omega_s$", labelpad = -8)
            ax.set_ylabel(r"$\Omega_\phi$", labelpad = -12)
            ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))    
            if (name == 'PoutPZ_Avg') or (name == 'PoutEM_Avg') or (name == 'TotalPout'):
                cbar.ax.set_title(r"$\times 10^{-3}$", fontsize = 6)
        # =========================================================================== #
        #                                Save Figure                                  #
        # =========================================================================== #
        pltconf.savefigure(savepath, savename, '.png', fig, save = save)
        plt.show()
        plt.close(fig)
        if save == True:
            print(f"saved: O={Omega:.2f}, g={gamma:.1f}")