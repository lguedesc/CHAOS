import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from matplotlib.colors import BoundaryNorm, ListedColormap
import matplotlib.colors as mpl_col
from mpl_toolkits.axes_grid1 import make_axes_locatable

# =========================================================================== #
#                    Functions to Configure Plot Parameters                   #
# =========================================================================== #
def plot_params(tex, fontsize, lineweight):
    mpl.rcParams['agg.path.chunksize'] = 10000
    plt.rcParams.update({
        'font.family': 'serif',  # use serif/main font for text elements
        'text.usetex': tex,     # use inline math for ticks
        'axes.titlesize': 'medium',
        'font.size': fontsize,
        'axes.linewidth' : lineweight,
        'xtick.major.width': lineweight,
        'ytick.major.width': lineweight,
        'xtick.major.size': 1.75,
        'ytick.major.size': 1.75,
        'xtick.top': True,
        'ytick.right': True,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
    })

def convert_dir(dir):
    if os.name == 'nt':      # For Windows
        if (dir.find('/') > 0):
            dir.replace('/', '\\')
    else:                   # For Linux and MacOs
        if (dir.find('\\') > 0):
            dir.replace('\\', '/')
    return dir

def set_colormap(df):
    if (df.max() >= 8):
        colormap = ListedColormap(['#404040', '#FFEE00', '#00DC00', '#FF8000', '#9900FF', '#007BFF', '#FF0000', '#700000'])
    else:
        colormap = ListedColormap(['#404040', '#FFEE00', '#00DC00', '#FF8000', '#9900FF', '#007BFF', '#FF0000'])
    return colormap

# =========================================================================== #
#                    Functions to Configure Colorbars                         #
# =========================================================================== #
def configure_colorbar_lyap(z):
    N = 128
    levels1 = np.linspace(z.min(), 0, N)
    levels2 = np.linspace(0, z.max(), N)
    levels = np.hstack((levels1,levels2[1:]))
    
    # Sample the right number of colours 
    # from the right bits (between 0 &amp; 1) of the colormaps we want.
    cmap1 = plt.get_cmap('binary')
    cols1 = cmap1(np.linspace(0, 1, N))
    
    cmap2 = mpl_col.LinearSegmentedColormap.from_list('Origin', ['darkviolet','blue','cyan','limegreen','yellow','darkorange','red','darkred'])
    cols2 = cmap2(np.linspace(0, 1, N))
    
    # Combine them and build a new colormap:
    allcols = np.vstack( (cols1,cols2) )
    cmap = mpl_col.LinearSegmentedColormap.from_list('custom_map', allcols)
        
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)    
    
    return cmap, norm

def configure_colormap_motion(z):
    #c_list = ['#404040', '#FFEE00', '#00DC00', '#FF8000', '#9900FF', '#007BFF', '#FF0000','#700000']
    #colormap = ListedColormap([c_list[0],c_list[1], c_list[2], c_list[3], c_list[4], c_list[5], c_list[6], c_list[7]])
    c_list = ['#404040', '#FFEE00', '#00DC00', '#FF8000', '#9900FF', 'red']
    colormap = ListedColormap([c_list[0],c_list[1], c_list[2], c_list[3], c_list[4], c_list[5]])
    cmap_min = z.min()
    cmap_max = z.max()
    return colormap, cmap_min, cmap_max

def configure_rainbow_colormap():
    cmap = mpl_col.LinearSegmentedColormap.from_list('Origin', ['darkviolet','blue','cyan','limegreen','yellow','darkorange','red','darkred'])
    return cmap

# =========================================================================== #
#                    Functions to Plot Diagrams                               #
# =========================================================================== #
def plot_lyap_map(ax, x, y, z):
    raster = True
    ax.set_aspect('equal')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size = '2.5%', pad = 0.05)    
    lsize = 7.5
    
    if (z.min() < 0 and z.max() > 0):    
        colormap, norm = configure_colorbar_lyap(z)
        plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = raster, cmap = colormap, norm = norm)
        cbar = plt.colorbar(plot, ax = ax, cax = cax, orientation='vertical', format = '${%.2f}$', ticks = [z.min(), 0, z.max()])
    elif (z.min() > 0 and z.max() > 0):
        colormap = configure_rainbow_colormap()
        plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = raster, cmap = colormap, vmin = z.min(), vmax = z.max())
        cbar = plt.colorbar(plot, ax = ax, cax = cax, orientation='vertical', format = '${%.2f}$', ticks = [z.min(), (z.max() - z.min())/2, z.max()])
    else:
        plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = raster, cmap = 'binary', vmin = z.min(), vmax = z.max())
        cbar = plt.colorbar(plot, ax = ax, cax = cax, orientation='vertical', format = '${%.2f}$', ticks = [z.min(), (z.max() - z.min())/2, z.max()])
    cbar.ax.tick_params(labelsize = lsize)
    cbar.ax.yaxis.set_ticks_position('right')

def plot_attractor_map(ax, x, y, z, maxper):
    colormap, cmapmin, cmapmax = configure_colormap_motion(z)
    ax.set_aspect('equal')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size = '2.5%', pad = 0.05)    
    lsize = 7.5
    plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = True, cmap = colormap, vmin = cmapmin - 0.5, vmax=cmapmax + 0.5)
    # Copy array without duplicates and sort it 
    ticks = np.sort(z.ravel())
    ticks = np.unique(ticks)
    for i in ticks:
        labels = []
        for i in ticks:
            if (i < maxper):
                name = f"{i}T"
                labels.append(name)
            elif (i == maxper):
                name = f"MP"
                labels.append(name)

    cbar = plt.colorbar(plot, ax=ax, cax = cax, orientation = 'vertical', ticks = ticks, drawedges = True)
    cax.set_yticklabels(labels)
    cbar.ax.tick_params(size = 0, labelsize = lsize)
    cbar.outline.set_edgecolor('white')
    cbar.outline.set_linewidth(1)
    cbar.dividers.set_color('white')
    cbar.dividers.set_linewidth(1)
    cbar.ax.yaxis.set_ticks_position('right')
    
    
def plot_rainbow_map(ax, x, y, z):
    raster = True
    ax.set_aspect('equal')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size = '2.5%', pad = 0.05)    
    lsize = 7.5
    delta_z = abs(z.max()) - abs(z.min())
    cmap = configure_rainbow_colormap()
    plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = raster, cmap = cmap, vmin = z.min(), vmax = z.max())
    cbar = plt.colorbar(plot, ax = ax, cax = cax, orientation='vertical', format = '${%.2f}$', ticks = [z.min(), delta_z/2, z.max()])
    cbar.ax.tick_params(labelsize = lsize)
    cbar.ax.yaxis.set_ticks_position('right')

def plot_neg_pos_map(ax, x, y, z):
    raster = True
    ax.set_aspect('equal')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size = '2.5%', pad = 0.05)    
    lsize = 7.5
    delta_z = abs(z.max()) - abs(z.min())
    #cmap = mpl.cm.coolwarm
    cmap = "RdGy_r"
    plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = raster, cmap = cmap, norm = mpl_col.CenteredNorm())
    #cbar = plt.colorbar(plot, ax = ax, cax = cax, orientation='vertical', format = '${%.2f}$', ticks = [z.min(), delta_z/2, z.max()])
    cbar = plt.colorbar(plot, ax = ax, cax = cax, orientation='vertical', format = '${%.2f}$')
    cbar.ax.tick_params(labelsize = lsize)
    cbar.ax.yaxis.set_ticks_position('right')
# =========================================================================== #
#                    Functions to Handle Data                                 #
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
#                               Misc Functions                                #
# =========================================================================== #
def mark_point_in_plot(ax, Omega, gamma, number, color):
    ax.scatter([Omega],[gamma],s = 20,linewidth = 0.5, facecolors = 'none', edgecolors = color)
    ax.text(Omega + 0.03, gamma - 0.09, number, color = color, fontsize=9)

