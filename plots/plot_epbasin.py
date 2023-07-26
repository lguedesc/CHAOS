# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt    
import matplotlib.colors as mpl_col
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
from libs import plotconfig as pltconf

pltconf.plot_params(True, 10, 0.5)

# =========================================================================== #
#                   Create functions to read and handle data                  #
# =========================================================================== #
def read_data(filename):
    df = pd.read_csv(filename, delimiter = " ", dtype="float")

    return df
    
def process_data(data, y_axis, x_axis, z_axis):
    x = data[x_axis].drop_duplicates()
    y = data[y_axis].drop_duplicates()
    z = data.pivot(index = y_axis, columns = x_axis, values = z_axis).values
    
    return x, y, z

# =========================================================================== #
#               Function to Configure Attractors Colorbar 
# =========================================================================== #
def configure_colormap_motion():
    c_list = ['#404040', 'orange', 'red', 'yellow', 'cyan', 'green']
    colormap = ListedColormap([c_list[0],c_list[1],c_list[2],c_list[3], c_list[4], c_list[5]])
    cmap_min = 1
    cmap_max = 6
    return colormap, cmap_min, cmap_max, c_list
# =========================================================================== #
#                           Function to change axis labels 
# =========================================================================== #
def customize_labels(ax, title1):
    ax.set_title(title1, loc = 'left')
    #ax.set_title(title2, loc = 'center')
    ax.set_xlabel(r'$x_1$', labelpad = -10)
    ax.set_xticks([x1.min(), (x1.max() - x1.min())/2, x1.max()])
    ax.set_yticks([y1.min(), (y1.max() - y1.min())/2, y1.max()])
    ax.set_ylabel(r'$x_2$', labelpad = -10)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
# =========================================================================== #
#                               Function to Plot 
# =========================================================================== #
def plot_maps(ax, x, y, z, colormap):
    ax.set_aspect('equal')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size = '2.5%', pad = 0.05)    
    lsize = 7.5
    plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = raster, cmap = colormap, vmin = cmap1_min - 0.5, vmax=cmap1_max + 0.5)
    ticks = [1, 2, 3, 4, 5, 6]
    #labels = ['EP1','EP2','EP3','EP4','EP5']
    cbar = plt.colorbar(plot, ax=ax, cax = cax, orientation = 'vertical', ticks = ticks, drawedges = True)
    #cax.set_yticklabels(labels)
    cbar.ax.tick_params(size = 0, labelsize = lsize)
    cbar.outline.set_edgecolor('white')
    cbar.outline.set_linewidth(1)
    cbar.dividers.set_color('white')
    cbar.dividers.set_linewidth(1)
    cbar.ax.yaxis.set_ticks_position('right')
# =========================================================================== #
#                                    Read Data                                #
# =========================================================================== #
save = False
system = "duffing"
ext = ".pdf"

filenum = 1

if filenum == 0:
    readpath = "data/EPBasin/out/" + system + f"_epbasin.csv"; readpath = pltconf.convert_dir(readpath)
else:
    readpath = "data/EPBasin/out/" + system + f"_epbasin({filenum}).csv"; readpath = pltconf.convert_dir(readpath)

savepath = "data/EPBasin/figs"; savepath = pltconf.convert_dir(savepath)
raw_data = pd.read_csv(readpath, delimiter = " ")

x1, y1, z1 = process_data(raw_data, 'CparY', 'CparX', 'Attractor')
# =========================================================================== #
#                Create custom colormaps and define colormap parameters       #
# =========================================================================== #
colormap1, cmap1_min, cmap1_max, c_list = configure_colormap_motion()
# =========================================================================== #
#                           Define figure parameters                          #
# =========================================================================== #
cm = 1/2.54
x_inches = 6.5     # [mm]*constant
y_inches = x_inches*(1)
raster = True
if save == True:
    dpi = 300
else:
    dpi = 300
# =========================================================================== #
#                               Create figure 
# =========================================================================== #
fig = plt.figure(1, figsize = (x_inches*cm,y_inches*cm), dpi = dpi)#, constrained_layout = True)
#fig.tight_layout()
#fig.set_constrained_layout_pads(hspace = 0, wspace = 0.1, w_pad = 0, h_pad = 0)

fig.subplots_adjust(top=0.9,
                    bottom=0.094,
                    left=0.098,
                    right=0.927,
                    hspace=0.1,
                    wspace=0.26)


lin = 1; col = 1
ax1 = fig.add_subplot(lin,col,1)
# =========================================================================== #
#                               Plot Data  
# =========================================================================== #
plot_maps(ax1, x1, y1, z1, colormap1)
# =========================================================================== #
#                          Customize titles and labels                        #
# =========================================================================== #
customize_labels(ax1, r'Attractors (Equilibrium Points)')
# =========================================================================== #
#                                Save Figure                                  #
# =========================================================================== #
if save == True:
    isExist = os.path.exists(savepath)
    if (isExist == False):
        os.makedirs(savepath)
    
    name = "/sample_" + system + "_EPBasin" + ext; name = pltconf.convert_dir(name)
    fig.savefig(savepath + name)
    plt.show()
else:
    plt.show()

