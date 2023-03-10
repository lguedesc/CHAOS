import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from matplotlib.colors import BoundaryNorm, ListedColormap
import matplotlib.colors as mpl_col
from mpl_toolkits.axes_grid1 import make_axes_locatable

# =========================================================================== #
#                    Functions to Handle Data                                 #
# =========================================================================== #
def convert_dir(dir):
    if os.name == 'nt':      # For Windows
        if (dir.find('/') > 0):
            dir.replace('/', '\\')
    else:                   # For Linux and MacOs
        if (dir.find('\\') > 0):
            dir.replace('\\', '/')
    return dir

def read_data(filename):
    df = pd.read_csv(filename, delimiter = " ", dtype="float")

    return df

def define_readpath(foldername, simulationname, filenum, system):
    if (filenum > 0):
        readpath = f"data/{foldername}/out/{system}_{simulationname}({filenum}).csv"; readpath = convert_dir(readpath)    
    else:
        readpath = f"data/{foldername}/out/{system}_{simulationname}.csv"; readpath = convert_dir(readpath)
    return readpath

def read_CHAOS_data(system, filenum, simulation):
    if simulation == 'timeseries':
        readpath = define_readpath('TimeSeries', simulation, filenum, system)
        readpath = convert_dir(readpath)
        df = read_data(readpath)
        return df
    elif simulation == 'poinc':
        readpath = define_readpath('PoincareMap', simulation, filenum, system)
        readpath = convert_dir(readpath)
        df = read_data(readpath)
        return df
    elif simulation == 'lyap':
        readpath = define_readpath('LyapunovExp', simulation, filenum, system)
        readpath = convert_dir(readpath)
        df = read_data(readpath)
        return df
    elif simulation == 'ftimeseries':
        readpath = define_readpath('FTimeSeries', simulation, filenum, system)
        readpath_poinc = define_readpath('FTimeSeries', 'poinc', filenum, system)
        readpath = convert_dir(readpath)
        readpath_poinc = convert_dir(readpath_poinc)
        df = read_data(readpath)
        df_poinc = read_data(readpath_poinc)
        return df, df_poinc
    elif simulation == 'bifurc':
        readpath = define_readpath('Bifurcation', simulation, filenum, system)
        readpath_poinc = define_readpath('Bifurcation', 'bifurc_poinc', filenum, system)
        readpath = convert_dir(readpath)
        readpath_poinc = convert_dir(readpath_poinc)
        df = read_data(readpath)
        df_poinc = read_data(readpath_poinc)
        return df, df_poinc
    elif simulation == 'fbifurc':
        readpath = define_readpath('FBifurcation', simulation, filenum, system)
        readpath_poinc = define_readpath('FBifurcation', 'fbifurc_poinc', filenum, system)
        readpath = convert_dir(readpath)
        readpath_poinc = convert_dir(readpath_poinc)
        df = read_data(readpath)
        df_poinc = read_data(readpath_poinc)
        return df, df_poinc
    elif simulation == 'dyndiag':
        readpath = define_readpath('DynDiagram', simulation, filenum, system)
        readpath = convert_dir(readpath)
        df = read_data(readpath)
        return df
    elif simulation == 'fdyndiag':
        readpath = define_readpath('FDynDiagram', simulation, filenum, system)
        readpath = convert_dir(readpath)
        df = read_data(readpath)
        return df
    elif simulation == 'epbasin':
        readpath = define_readpath('EPBasin', simulation, filenum, system)
        readpath = convert_dir(readpath)
        df = read_data(readpath)
        return df
    elif simulation == 'fforcedbasin':
        readpath = define_readpath('FForcBasin', simulation, filenum, system)
        readpath = convert_dir(readpath)
        df = read_data(readpath)
        return df
        
def process_data(data, y_axis, x_axis, z_axis):
    x = data[x_axis].drop_duplicates()
    y = data[y_axis].drop_duplicates()
    z = data.pivot(index = y_axis, columns = x_axis, values = z_axis).values
    
    return x, y, z

def define_savepath(foldername, simulationname, system, ext):
    savepath = f"data/{foldername}/figs/"
    savepath = convert_dir(savepath)
    fullpath = f"data/{foldername}/figs/{system}_{simulationname}{ext}"
    fullpath = convert_dir(fullpath)    
    
    return savepath, fullpath

def save_CHAOS_data(fig, system, simulation, ext):
    savepath = ''
    fullpath = ''
    if simulation == 'timeseries':
        savepath, fullpath = define_savepath('TimeSeries', simulation, system, ext)
    elif simulation == 'poinc':
        savepath, fullpath = define_savepath('PoincareMap', simulation, system, ext)
    elif simulation == 'lyap':
        savepath, fullpath = define_savepath('LyapunovExp', simulation, system, ext)
    elif simulation == 'ftimeseries':
        savepath, fullpath = define_savepath('FTimeSeries', simulation, system, ext)
    elif simulation == 'bifurc':
        savepath, fullpath = define_savepath('Bifurcation', simulation, system, ext)
    elif simulation == 'fbifurc':
        savepath, fullpath = define_savepath('FBifurcation', simulation, system, ext)
    elif simulation == 'dyndiag':
        savepath, fullpath = define_savepath('DynDiagram', simulation, system, ext)
    elif simulation == 'fdyndiag':
        savepath, fullpath = define_savepath('FDynDiagram', simulation, system, ext)
    elif simulation == 'epbasin':
        savepath, fullpath = define_savepath('EPBasin', simulation, system, ext)
    elif simulation == 'fforcedbasin':
        savepath, fullpath = define_savepath('FForcBasin', simulation, system, ext)
    # Check if directory exists
    isExist = os.path.exists(savepath)
    if (isExist == False):
        os.makedirs(savepath)
    # Save Figure
    fig.savefig(fullpath)

def handle_comparison_data(raw_data_negative, raw_data_positive, column_nameX, column_nameY, column_nameZ):
    comp_data_name = 'comparison'
    comp_data = raw_data_negative.loc[:,(column_nameY, column_nameX, column_nameZ)]
    comp_data.rename(columns={column_nameZ:'data_for_negative'}, inplace = True)
    comp_data = pd.concat([comp_data, raw_data_positive[column_nameZ]], axis = 1)
    comp_data[comp_data_name] = ((comp_data[column_nameZ] - comp_data['data_for_negative'])/comp_data['data_for_negative'])*100
    comp_data[comp_data_name] = comp_data[comp_data_name].fillna(0)
    
    return comp_data, comp_data_name     
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
        'axes.xmargin': 0.035,
        'axes.ymargin': 0.035,
    })

def makefig_and_axs(figsize, rows, cols, dpi, hspace = 0, wspace = 0, hpad = 0, wpad = 0):
    fig, axs = plt.subplots(rows, cols, figsize = (figsize[0], figsize[1]), dpi = dpi, layout="constrained")
    fig.get_layout_engine().set(w_pad=wpad, h_pad=hpad, hspace=hspace, wspace=wspace)
    return fig, axs

def makefig_and_grid(figsize, rows, cols, dpi, hspace = 0, wspace = 0, hpad = 0, wpad = 0):
    fig = plt.figure(figsize = (figsize[0], figsize[1]), dpi = dpi, layout="constrained")
    fig.get_layout_engine().set(w_pad=wpad, h_pad=hpad, hspace=hspace, wspace=wspace)
    grid = fig.add_gridspec(rows, cols)
    return fig, grid
    
def figsize_in_cm(xcm, ycm):
    const = 1/2.54
    figsize_in = [xcm*const, ycm*const] # 1[in] = 1[cm]*constant    
    
    return figsize_in
    
def set_colormap(df):
    if (df.max() >= 8):
        colormap = ListedColormap(['#404040', '#FFEE00', '#00DC00', '#FF8000', '#9900FF', '#007BFF', '#FF0000', '#700000'])
    else:
        colormap = ListedColormap(['#404040', '#FFEE00', '#00DC00', '#FF8000', '#9900FF', '#007BFF', '#FF0000'])
    return colormap

def remove_some_colors_from_colormap(min_val, max_val, N, cmap):
    orig_cmap = plt.get_cmap(cmap)
    colors = orig_cmap(np.linspace(min_val, max_val, N))
    cmap = mpl_col.LinearSegmentedColormap.from_list('mycmap', colors)

    return cmap

def set_fig_quality(save = False, base_dpi = 100):
    if save == True:
        dpi = 6*base_dpi
    else:
        dpi = base_dpi
    return dpi
# =========================================================================== #
#                    Functions to Configure Colorbars                         #
# =========================================================================== #
def configure_colorbar_comparison(zmin, zmax):
    N = 128
    levels1 = np.linspace(zmin, 0, N)
    levels2 = np.linspace(0, zmax, N)
    levels = np.hstack((levels1,levels2[1:]))
    
    # Sample the right number of colours 
    # from the right bits (between 0 &amp; 1) of the colormaps we want.    
    cmap1 = plt.get_cmap('Greys_r')
    cols1 = cmap1(np.linspace(0, 0.3, N))
    
    #cmap1 = plt.get_cmap('Reds_r')
    #cols1 = cmap1(np.linspace(0, 0.3, N))
    
    #cmap2 = mpl_col.LinearSegmentedColormap.from_list('Origin', ['darkviolet','blue','cyan','limegreen','yellow','darkorange','red','darkred'])
    cmap2 = plt.get_cmap('Reds')
    cols2 = cmap2(np.linspace(0.4, 0.7, N))
    
    #cmap2 = plt.get_cmap('Greys')
    #cols2 = cmap2(np.linspace(0.7, 1.0, N))
    
    # Combine them and build a new colormap:
    allcols = np.vstack( (cols1,cols2) )
    cmap = mpl_col.LinearSegmentedColormap.from_list('custom_map', allcols)
        
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)    
    
    return cmap, norm

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

def configure_new_colormap_motion(z):
    colors = ['#404040', 'black', '#FFEE00', 'goldenrod', '#00DC00', 'forestgreen', '#FF8000', 'chocolate', 'violet', '#9900FF',
              '#007BFF', 'blue', 'red', '#700000', 'magenta', 'deeppink', 'white'] 
    # Copy array without duplicates and sort it 
    ticks = np.sort(z.ravel())
    ticks = np.unique(ticks)
    # Attribute designated colors to each value of tick
    newclrs = []
    for tks in range(len(ticks)):
        newclrs.append(colors[tks])
    # Make colormap
    colormap = ListedColormap([*newclrs])
    cmap_min = z.min()
    cmap_max = z.max()
    return colormap, cmap_min, cmap_max, ticks

def configure_colormap_motion(z):
    colors = ['#404040', '#FFEE00', '#00DC00', '#FF8000', '#9900FF', '#007BFF', '#FF0000','#700000'] 
    # Copy array without duplicates and sort it 
    ticks = np.sort(z.ravel())
    ticks = np.unique(ticks)
    # Attribute designated colors to each value of tick
    newclrs = []
    for tks in range(len(ticks)):
        newclrs.append(colors[tks])
    # Make colormap
    colormap = ListedColormap([*newclrs])
    cmap_min = z.min()
    cmap_max = z.max()
    return colormap, cmap_min, cmap_max, ticks

def configure_rainbow_colormap():
    cmap = mpl_col.LinearSegmentedColormap.from_list('Origin', ['darkviolet','blue','cyan','limegreen','yellow','darkorange','red','darkred'])
    return cmap

# =========================================================================== #
#                    Functions to Plot Diagrams                               #
# =========================================================================== #
def plot_lyap_map(fig, ax, x, y, z):
    raster = True
    #ax.set_aspect('equal') 
    lsize = 6
    aspect = 20
    pad = 0.025
    if (z.min() < 0 and z.max() > 0):    
        colormap, norm = configure_colorbar_lyap(z)
        plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = raster, cmap = colormap, norm = norm)
        cbar = fig.colorbar(plot, ax = ax, orientation='vertical', aspect = aspect, pad = pad, format = '${%.2f}$', ticks = [z.min(), 0, z.max()])
    elif (z.min() > 0 and z.max() > 0):
        colormap = configure_rainbow_colormap()
        plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = raster, cmap = colormap, vmin = z.min(), vmax = z.max())
        cbar = fig.colorbar(plot, ax = ax, orientation='vertical', aspect = aspect, pad = pad, format = '${%.2f}$', ticks = [z.min(), (z.max() - z.min())/2, z.max()])
    else:
        plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = raster, cmap = 'binary', vmin = z.min(), vmax = z.max())
        cbar = fig.colorbar(plot, ax = ax, orientation='vertical', aspect = aspect, pad = pad, format = '${%.2f}$', ticks = [z.min(), (z.max() - z.min())/2, z.max()])
    cbar.ax.tick_params(labelsize = lsize)
    cbar.ax.yaxis.set_ticks_position('right')
    cbar.ax.minorticks_off()
    
    return plot, cbar
    
def plot_attractor_map(fig, ax, x, y, z, maxper, mode = 'lyap'):
    colormap, cmapmin, cmapmax, ticks = configure_colormap_motion(z)   
    #ax.set_aspect('equal')
    #ax.set_box_aspect(1)
    lsize = 6
    plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = True, cmap = colormap, 
                         vmin = cmapmin - 0.5, vmax=cmapmax + 0.5)    
    '''
    if mode == 'lyap':
        ticks = np.linspace(1, maxper+2, maxper+2, dtype=int)     
        labels = []
        for i in ticks:
            if (i < maxper):
                name = f"{i}T"
                labels.append(name)
            elif (i == maxper):
                name = f"MP"
                labels.append(name)
            elif (i == maxper + 1):
                name = f"CH"
                labels.append(name)
            elif (i == maxper + 2):
                name = f"HC"
                labels.append(name)    
            print(labels)
    else:
        ticks = np.linspace(1, maxper, maxper, dtype=int)
        labels = []
        for i in ticks:
            if (i < maxper):
                name = f"{i}T"
                labels.append(name)
            elif (i == maxper):
                name = f"MP"
                labels.append(name)
    '''
    cbar = fig.colorbar(plot, ax=ax, location = 'right', orientation='vertical', aspect = 20, pad = 0.01,
                        ticks = ticks, drawedges = True)
    #cbar.set_ticklabels(labels)
    cbar.ax.tick_params(size = 0, labelsize = lsize)
    cbar.outline.set_edgecolor('white')
    cbar.outline.set_linewidth(1)
    cbar.dividers.set_color('white')
    cbar.dividers.set_linewidth(1)
    cbar.ax.yaxis.set_ticks_position('right')

def plot_new_attractor_map(fig, ax, x, y, z, maxper, mode = 'lyap'):
    colormap, cmapmin, cmapmax, ticks = configure_new_colormap_motion(z)   
    #ax.set_aspect('equal')
    #ax.set_box_aspect(1)
    lsize = 6
    plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = True, cmap = colormap, 
                         vmin = cmapmin - 0.5, vmax=cmapmax + 0.5)    
    '''
    if mode == 'lyap':
        #ticks = np.linspace(1, maxper+2, maxper+2, dtype=int)     
        ticks = np.linspace(1, 2*(maxper+2)+1, 2*(maxper+2)+1, dtype=int)
        letters = ['L', 'W', 'L', 'W', 'L', 'W', 'L', 'W', 'L', 'W', 'L', 'W', 'L', 'W', 'L', 'W', '']
        print(ticks)
        
        # 12 34 56 78 910 1112 1314 1516 17 
        # 1T 2T 3T 4T 5T  MP   CH   HC   ESCP
        labels = []
        for i, ltr in zip(ticks, letters):
            if (i < 2*(maxper) - 1): #10  
                name = f"{i}T-{ltr}"
                labels.append(name)
            elif (i == 2*maxper or i == 2*maxper - 1): # 12 or 11
                name = f"MP-{ltr}"
                labels.append(name)
            elif (i == 2*maxper + 1 or i == 2*maxper + 2):  #13 14
                name = f"CH-{ltr}"
                labels.append(name)
            elif (i == 2*maxper + 3 or i == 2*maxper + 4):    #14
                name = f"HC-{ltr}"
                labels.append(name)
            else:                        #17
                name = f"ESCP"
                labels.append(name)    
                
            print(labels)

    else:
        ticks = np.linspace(1, maxper, maxper, dtype=int)
        labels = []
        for i in ticks:
            if (i < maxper):
                name = f"{i}T"
                labels.append(name)
            elif (i == maxper):
                name = f"MP"
                labels.append(name)
'''

    cbar = fig.colorbar(plot, ax=ax, location = 'right', orientation='vertical', aspect = 20, pad = 0.01,
                        ticks = ticks, drawedges = True)
    #cbar.set_ticklabels(labels)
    cbar.ax.tick_params(size = 0, labelsize = lsize)
    cbar.outline.set_edgecolor('white')
    cbar.outline.set_linewidth(1)
    cbar.dividers.set_color('white')
    cbar.dividers.set_linewidth(1)
    cbar.ax.yaxis.set_ticks_position('right')

def plot_rainbow_map_old(ax, x, y, z):
    raster = True
    #ax.set_aspect('equal')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size = '2.5%', pad = 0.05)    
    #lsize = 7.5
    delta_z = abs(z.max()) - abs(z.min())
    cmap = configure_rainbow_colormap()
    plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = raster, cmap = cmap, vmin = z.min(), vmax = z.max())
    cbar = plt.colorbar(plot, ax = ax, cax = cax, orientation='vertical', format = '${%.2f}$', ticks = [z.min(), delta_z/2, z.max()])
    #cbar.ax.tick_params(labelsize = lsize)
    #cbar.ax.tick_params()
    cbar.ax.yaxis.set_ticks_position('right')
    
def plot_rainbow_map(fig, ax, x, y, z, colormap = 'myrainbow', clr_min = 0, clr_max = 1):
    raster = True
    #ax.set_aspect('equal')
    pad = 0.025
    aspect = 20
    lsize = 6
    delta_z = abs(z.max()) - abs(z.min())
    if colormap == 'myrainbow':
        cmap = configure_rainbow_colormap()
    else:
        N = len(z)
        cmap = remove_some_colors_from_colormap(clr_min, clr_max, N, colormap)
    plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = raster, cmap = cmap, vmin = z.min(), vmax = z.max())
    cbar = fig.colorbar(plot, ax = ax, location = 'right', orientation='vertical', aspect = aspect, pad = pad,
                        format = '${%.2f}$', ticks = [z.min(), delta_z/2, z.max()])
    cbar.ax.tick_params(labelsize = lsize)
    #cbar.ax.tick_params()
    cbar.ax.yaxis.set_ticks_position('right')
    
    return plot, cbar

def plot_neg_pos_map(ax, x, y, z):
    raster = True
    #ax.set_aspect('equal')
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

def plot_comparison_map(fig, ax, x, y, z, zmin, zmax, orientation = 'vertical'):
    raster = True
    #ax.set_aspect('equal') 
    lsize = 6
    aspect = 20
    pad = 0.025
    format = '${%.0f}$'
    if zmax == z.max():
        maxvalue = z.max()
        biggersymb = ''
        extend = 'neither'
    else:
        maxvalue = zmax   
        biggersymb = '>'
        extend = 'max'
    if zmin == z.min():
        minvalue = z.min()
        lessersymb = ''
        extend = 'neither'
    else:
        minvalue = zmin   
        lessersymb = '<'
        extend = 'min'
    
    if (zmin != z.min() and zmax != z.max()):    
        extend = 'both'
        
    if (z.min() < 0 and z.max() > 0):
        colormap, norm = configure_colorbar_comparison(zmin, zmax)
                
        plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = raster, cmap = colormap, norm = norm)
        cbar = fig.colorbar(plot, ax = ax, orientation=orientation, aspect = aspect, pad = pad, format = format, 
                            extend = extend, ticks = [minvalue, 0, maxvalue])
        cbar.ax.set_yticklabels([rf'${lessersymb}$' + f'${minvalue:.0f}\%$', f'$0\%$', rf'${biggersymb}$' + f'${maxvalue:.0f}\%$'])
    elif (z.min() > 0 and z.max() > 0):
        colormap = configure_rainbow_colormap()
        plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = raster, cmap = colormap, vmin = minvalue, vmax = maxvalue)
        cbar = fig.colorbar(plot, ax = ax, orientation=orientation, aspect = aspect, pad = pad, format = format, 
                            extend = extend, ticks = [minvalue, (maxvalue - minvalue)/2, maxvalue])
        cbar.ax.set_xticklabels([rf'${lessersymb}$' + f'${minvalue:.1f}\%$', f'{(maxvalue - minvalue)/2:.2f}', biggersymb + f'${maxvalue:.2f}\%$'])
    else:
        plot = ax.pcolormesh(x, y, z, shading = 'nearest', rasterized = raster, cmap = 'binary', vmin = minvalue, vmax = maxvalue)
        cbar = fig.colorbar(plot, ax = ax, orientation=orientation, aspect = aspect, pad = pad, format = format, 
                            extend = extend, ticks = [minvalue, (maxvalue - minvalue)/2, maxvalue])
        cbar.ax.set_xticklabels([rf'${lessersymb}$' + f'${minvalue:.1f}\%$', f'{(maxvalue - minvalue)/2:.2f}', biggersymb + f'${maxvalue:.2f}\%$'])
    
    cbar.ax.tick_params(labelsize = lsize)
    cbar.ax.minorticks_off()
    
    if orientation == 'horizontal':
        cbar.ax.xaxis.set_ticks_position('top')
        cbar.ax.xaxis.set_label_position('top')
        tl = cbar.ax.get_xticklabels()
        tl[0].set_horizontalalignment('left')
        tl[-1].set_horizontalalignment('right')
    else:
        cbar.ax.yaxis.set_ticks_position('right')
     
    return plot, cbar

# =========================================================================== #
#                               Misc Functions                                #
# =========================================================================== #
def mark_point_in_plot(ax, Omega, gamma, number, color):
    ax.scatter([Omega],[gamma],s = 20,linewidth = 0.5, facecolors = 'none', edgecolors = color)
    ax.text(Omega + 0.03, gamma - 0.09, number, color = color, fontsize=9)

