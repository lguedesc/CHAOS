import matplotlib.pyplot as plt
from Plots import plotconfig as pltconf
import numpy as np
import pandas as pd

def estimate_start_and_end_of_rotation(df, time_name, angle_name):
    # Check the minimum and maximum values of angle in radians
    rad_i = df[angle_name].min()
    rad_f = df[angle_name].max()
    # Find the multiples of pi
    pi_i = int(rad_i/np.pi) 
    pi_f = int(rad_f/np.pi)
    pi_multiples = np.linspace(pi_i, pi_f, pi_f - pi_i + 1)
    # Determine tha values to search
    values_to_search = pi_multiples*np.pi
    # Calculate the closest values and corresponding times
    closest_values = []
    closest_times = []

    for value in values_to_search:
        closest_index = (df[angle_name] - value).abs().idxmin()
        closest_values.append(df.at[closest_index, angle_name])
        closest_times.append(df.at[closest_index, time_name])

    # Create a new DataFrame
    result_df = pd.DataFrame({'Value_to_Search': values_to_search, 'Closest_angle': closest_values, 'Closest_time': closest_times})

    return result_df['Closest_time']
    
system = "multidirect_hybrid_EH"
filenum = 1

save = False
show = True

# ============================================================================================== #
# Load Data
# ============================================================================================== #
headers = ['Time', 'x[0]', 'x[1]', 'x[2]', 'x[3]', 'x[4]', 'x[5]', 'x[6]', 'x[7]']
df = pltconf.read_data("data/FTimeSeries/out/multidirect_hybrid_EH_ftimeseries(1).csv", cols = headers)

estimated_pi_multiples = estimate_start_and_end_of_rotation(df, "Time", "x[4]")
# ============================================================================================== #
# Create figure and axes
# ============================================================================================== #
nrows = 8; ncols = 1
figwidth = 15
pltconf.plot_params(True, 10, 0.2)
fsize = pltconf.figsize_in_cm(figwidth, 0.15*nrows*figwidth)
dpi = pltconf.set_fig_quality(save = save, base_dpi=120)
fig, axs = pltconf.makefig_and_axs(figsize = fsize, rows = nrows, cols = ncols, dpi = dpi, hpad = 0.01)
# ============================================================================================== #
# Customize axes labels
# ============================================================================================== #
labels = [r"$\bar{x}$", r"$\dot{\bar{x}}$", r"$\bar{z}$", r"$\dot{\bar{z}}$", r"$\bar{\phi}$", r"$\dot{\bar{\phi}}$", r"$\bar{v}$", r"$\bar{I}$"]
for ax, lbl in zip(axs, labels):
    ax.set_ylabel(lbl)
axs[nrows-1].set_xlabel(r"$\tau$")    
# ============================================================================================== #
# Plot data
# ============================================================================================== #
colors = ['black', 'black', 'blue', 'blue', 'purple', 'purple', 'gold', 'orange']
df_headers = ['x[0]', 'x[1]', 'x[2]', 'x[3]', 'x[4]', 'x[5]', 'x[6]', 'x[7]']
for ax, hdr, clr in zip(axs, df_headers, colors):
   ax.plot(df["Time"], df[hdr], color = clr, zorder = 0)
   for i in range(len(estimated_pi_multiples)):
       ax.axvline(x=estimated_pi_multiples[i], lw = 0.5, color = 'black', ls = '--', zorder = 1)
   ax.set_xlim([3650, df["Time"].iloc[-1]])

axs[4].set_ylim(-20120, df["x[4]"].max())

plt.show()