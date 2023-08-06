# -*- coding: utf-8 -*-
#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sympy as sym
from prettytable import PrettyTable
plt.rcParams.update({
    'font.family': 'serif',     # use serif/main font for text elements
    'text.usetex': True,        # use inline math for ticks
    'axes.titlesize': 'medium',
    'font.size': 10
})
plt.close('all')
sym.init_printing()

def potential(X, Z):
    Pot = ((g*(mp + ms) + kx*X + (kpz + kz)*Z)*(X**2 + ((g*(mp + ms))/(kpz + kz) + Z)**2))/(2.*np.sqrt(X**2 + (g*(mp + ms) + (kpz + kz)*Z)**2/(kpz + kz)**2))
    return Pot 


# Main Parameters
g = 9.81
mp = 1
ms = 0.5
kx = 1
kpz = 0.5
kz = 0.5

save = False
# ========================================================================== #
# Define Contour Plot Equations and Meshgrid
# ========================================================================== #
npts = 1000
u1 = np.linspace(-2.0,2.0,npts)
u2 = np.linspace(-2.0,2.0,npts)
X, Z = np.meshgrid(u1, u2) 
Pot = potential(X, Z)
# ========================================================================== #
# Set Parameters and Create Figure 
# ========================================================================== #
x_inches = 85*(1/25.4)     # [mm]*constant
y_inches = x_inches*(0.8)
raster = True
if save == True:
    dpi = 300
else:
    dpi = 300
fig = plt.figure(1, figsize = (x_inches,y_inches), dpi = dpi, constrained_layout = True)
ax = fig.add_subplot(1,1,1)
#levels = np.linspace(Pot.min(), 0.0, 21)
plt.rcParams['contour.negative_linestyle'] = 'solid'
# Plot Potential Energy
cs = ax.contourf(X, Z, Pot, cmap=plt.get_cmap('viridis'), zorder=-20)
cb = fig.colorbar(cs, ax=ax)#, shrink=0.6)
cb.ax.set_title(r'$\bar{U}$')
ax.set_xlabel(r'$\bar{x}$')
ax.set_ylabel(r'$\bar{z}$')
# Adjusts to save image with good quality
ax.set_rasterization_zorder(-10)
if raster == True:
    for c in cs.collections:
        c.set_edgecolor("face")
#========================================================================#
# Show and Save Figure                                                   #
#========================================================================#
if save == True:
   plt.savefig(f'fig.pdf' )
plt.show()