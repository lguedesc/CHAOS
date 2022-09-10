# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt    
import matplotlib.colors as mpl_col
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
from src.libs import plotconfig as pltconf

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
#                    Functions to Configure Lyap Colorbar                     #
# =========================================================================== #
def configure_colorbar_lyap(z):
    N = 128
    levels1 = np.linspace(z.min(), 0, N)
    levels2 = np.linspace(0, z.max(), N)
    levels = np.hstack((levels1,levels2[1:]))
    
    # Sample the right number of colours 
    # from the right bits (between 0 &amp; 1) of the colormaps we want.
    cmap1 = cmap1 = plt.get_cmap('binary')
    cols1 = cmap1(np.linspace(0, 1, N))
    
    cmap2 = mpl_col.LinearSegmentedColormap.from_list('Origin', ['darkviolet','blue','cyan','limegreen','yellow','darkorange','red','darkred'])
    cols2 = cmap2(np.linspace(0, 1, N))
    
    # Combine them and build a new colormap:
    allcols = np.vstack( (cols1,cols2) )
    cmap = mpl_col.LinearSegmentedColormap.from_list('custom_map', allcols)
        
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)    
    
    return cmap, norm
# =========================================================================== #
#               Function to Configure Attractors Colorbar 
# =========================================================================== #
def configure_colormap_motion():
    c_list = ['#404040', '#FFEE00', '#00DC00', '#FF8000', '#9900FF', '#007BFF', '#FF0000','#700000']
    colormap = ListedColormap([c_list[0],c_list[1], c_list[2], c_list[3], c_list[4], c_list[5], c_list[6], c_list[7]])
    cmap_min = 1
    cmap_max = 8
    return colormap, cmap_min, cmap_max, c_list
# =========================================================================== #
#                           Function to change axis labels 
# =========================================================================== #
def customize_labels(ax, title1, custom = False):
    ax.set_title(title1, loc = 'left')
    #ax.set_title(title2, loc = 'center')
    if custom == '(a)':
        ax.set_xticks([0.01, 1, 2])
        ax.axes.xaxis.set_ticklabels([])
        ax.set_yticks([y2.min(), y2.max()])
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.set_ylabel(r'$\gamma$', labelpad = -10)
    elif custom == '(b)':
        ax.set_xticks([0.01, 1, 2])
        ax.axes.xaxis.set_ticklabels([])
        ax.set_yticks([y2.min(), y2.max()])
        ax.axes.yaxis.set_ticklabels([])
        #ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        #ax.set_ylabel(r'$\gamma$', labelpad = -10)
    elif custom == '(c)':
        ax.set_xlabel(r'$\Omega$')
        ax.set_xticks([0.01, 1, 2])
        ax.set_yticks([y2.min(), y2.max()])
        ax.set_ylabel(r'$\gamma$', labelpad = -10)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    elif custom == '(d)':
        ax.set_xlabel(r'$\Omega$')
        ax.set_xticks([0.01, 1, 2])
        ax.set_yticks([y2.min(), y2.max()])
        ax.axes.yaxis.set_ticklabels([])
        ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
# =========================================================================== #
#                           Function to write in plot 
# =========================================================================== #
def mark_point_in_plot(ax, Omega, gamma, number, color):
    ax.scatter([Omega],[gamma],s = 20,linewidth = 0.5, facecolors = 'none', edgecolors = color)
    ax.text(Omega + 0.03, gamma - 0.09, number, color = color, fontsize=9)
# =========================================================================== #
#                               Function to Plot 
# =========================================================================== #
def plot_maps(ax, x, y, z, colormap, norm, custom = False):
    ax.set_aspect('equal')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size = '2.5%', pad = 0.05)    
    lsize = 7.5
    if custom == 'attractors':
        plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = raster, cmap = colormap, vmin = cmap1_min - 0.5, vmax=cmap1_max + 0.5)
        ticks = [1,2,3,4,5,6,7,8]
        labels = ['P1','P2','P3','P4','P5','MP','Ch','HCh']
        cbar = plt.colorbar(plot, ax=ax, cax = cax, orientation = 'vertical', ticks = ticks, drawedges = True)
        cax.set_yticklabels(labels)
        cbar.ax.tick_params(size = 0, labelsize = lsize)
        cbar.outline.set_edgecolor('white')
        cbar.outline.set_linewidth(1)
        cbar.dividers.set_color('white')
        cbar.dividers.set_linewidth(1)
    elif custom == 'lyapunov':
        plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = raster, cmap = colormap, norm = norm)
        cbar = plt.colorbar(plot, ax = ax, cax = cax, orientation='vertical', format = '${%.2f}$', ticks = [z.min(), 0, z.max()])
        cbar.ax.tick_params(labelsize = lsize)
    else:
        delta_z = z.max() - z.min()
        #cmap = mpl_col.LinearSegmentedColormap.from_list('Origin', ['darkviolet','blue','cyan','limegreen','yellow','darkorange','red','darkred'])
        cmap = 'binary_r'
        plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = raster, cmap = cmap)#, vmax = 0.032)
        cbar = plt.colorbar(plot, ax=ax, cax = cax, orientation = 'vertical', format = '${%.2f}$', ticks = [z.min(),  0.25*delta_z, 0.5*delta_z, 0.75*delta_z, z.max()])
        cbar.ax.tick_params(labelsize = lsize)
    cbar.ax.yaxis.set_ticks_position('right')
# =========================================================================== #
#                                    Read Data                                #
# =========================================================================== #
save = False
system = "vanderpol"
ext = ".pdf"

readpath = "DynDiagram/out/" + system + "_dyndiag(4).csv"; readpath = pltconf.convert_dir(readpath)
savepath = "DynDiagram/figs"; savepath = pltconf.convert_dir(savepath)

raw_data = pd.read_csv(readpath, delimiter = " ")

x1, y1, z1 = process_data(raw_data, 'CparY', 'CparX', 'Attractor')
x2, y2, z2 = process_data(raw_data, 'CparY', 'CparX', 'LE[0]')
#x3, y3, z3 = process_data(raw_data, 'CparY', 'CparX', 'LE[1]')
#x4, y4, z4 = process_data(raw_data, 'CparY', 'CparX', 'LE[2]')
#x5, y5, z5 = process_data(raw_data, 'CparY', 'CparX', 'xRMS[2]')
#x6, y6, z6 = process_data(raw_data, 'CparY', 'CparX', 'OverallxRMS[2]')
#x7, y7, z7 = process_data(raw_data, 'CparY', 'CparX', 'xmax[0]')
#x8, y8, z8 = process_data(raw_data, 'CparY', 'CparX', 'xmin[0]')

# =========================================================================== #
#                Create custom colormaps and define colormap parameters       #
# =========================================================================== #
colormap1, cmap1_min, cmap1_max, c_list = configure_colormap_motion()
colormap2, norm2 = configure_colorbar_lyap(z2)
#colormap3, norm3 = configure_colorbar_lyap(z3)
# =========================================================================== #
#                           Define figure parameters                          #
# =========================================================================== #
cm = 1/2.54
x_inches = 15     # [mm]*constant
y_inches = x_inches*(1.2)
raster = True
if save == True:
    dpi = 300
else:
    dpi = 96
# =========================================================================== #
#                               Create figure 
# =========================================================================== #
fig = plt.figure(1, figsize = (x_inches*cm,y_inches*cm), dpi = dpi)#, constrained_layout = True)
#fig.tight_layout()
#fig.set_constrained_layout_pads(hspace = 0, wspace = 0.1, w_pad = 0, h_pad = 0)


fig.subplots_adjust(top=1,
                    bottom=0.08,
                    left=0.063,
                    right=0.927,
                    hspace=0.1,
                    wspace=0.26)


lin = 4; col = 2
ax1 = fig.add_subplot(lin,col,1)
ax2 = fig.add_subplot(lin,col,2)
ax4 = fig.add_subplot(lin,col,3)
ax3 = fig.add_subplot(lin,col,4)
ax5 = fig.add_subplot(lin,col,5)
ax6 = fig.add_subplot(lin,col,6)
ax7 = fig.add_subplot(lin,col,7)
ax8 = fig.add_subplot(lin,col,8)
# =========================================================================== #
#                               Plot Data  
# =========================================================================== #
plot_maps(ax1, x1, y1, z1, colormap1, norm2, custom = 'attractors')
#plot_maps(ax2, x2, y2, z2, colormap2, norm2, custom = 'lyapunov')
#plot_maps(ax3, x3, y3, z3, colormap1, norm2)
#plot_maps(ax4, x4, y4, z4, colormap1, norm2)
#plot_maps(ax5, x5, y5, z5, colormap1, norm2)
##plot_maps(ax6, x6, y6, z6, colormap1, norm2)
#plot_maps(ax7, x7, y7, z7, colormap1, norm2)
#plot_maps(ax8, x8, y8, z8, colormap1, norm2)
# =========================================================================== #
#                          Customize titles and labels                        #
# =========================================================================== #
customize_labels(ax1, r'(a) Dynamical Attractors', custom = '(a)')
customize_labels(ax2, r'(b) Largest Lyapunov Exponent ($\lambda_1$)', custom = '(b)')
customize_labels(ax3, r'(d) 2nd Lyapunov Exponent ($\lambda_2$)',  custom = '(b)')
customize_labels(ax4, r'(c) 3rd Lyapunov Exponent ($\lambda_2$)', custom = '(a)')
customize_labels(ax5, r'(e) xRMS[2]', custom = '(a)')
customize_labels(ax6, r'(f) Overall xRMS[2]', custom = '(b)')
customize_labels(ax7, r'(g) xmax[0]', custom = '(c)')
customize_labels(ax8, r'(h) xmin[0]', custom = '(d)')
# =========================================================================== #
#                                Save Figure                                  #
# =========================================================================== #
if save == True:
    isExist = os.path.exists(savepath)
    if (isExist == False):
        os.makedirs(savepath)
    
    name = "/sample_" + system + "_dyndiag(serial)" + ext; name = pltconf.convert_dir(name)
    fig.savefig(savepath + name)
    plt.show()
else:
    plt.show()

