import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

# Load the data from the CSV file
df = pd.read_csv("n_sweep_parellel_update_final.csv")

# Get unique values of n
n_values = df['n'].unique()

# Define the RGBY color scheme
colors = ['red', 'green', 'blue', 'yellow']
tree_types = ['Path', 'Binary', '64ary', 'Star']

# Create a figure with a custom grid for the subplots (3 rows, 3 columns)
fig = plt.figure(figsize=(16, 12)) # Adjusted figure size for better layout
gs = fig.add_gridspec(3, 3) # Changed to 3 rows, 3 columns
fs = 30
plt.rc('font', family='serif')
plt.rc('font', size=fs)          # controls default text sizes
plt.tight_layout()
plt.rc('axes', titlesize=fs * 0.8)     # fontsize of the axes title
plt.rc('axes', labelsize= fs * 1.1)    # fontsize of the x and y labels
plt.rc('xtick', labelsize= fs)    # fontsize of the tick labels
plt.rc('ytick', labelsize=fs * 0.85)    # fontsize of the tick labels (1.5x smaller than 120)
plt.rc('figure', titlesize=0.5 * fs)  # fontsize of the figure title
# Create the axes for the top row (3 plots)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])
top_axes = [ax1, ax2, ax3]

# Create the axes for the middle row (3 plots)
ax4 = fig.add_subplot(gs[1, 0])
ax5 = fig.add_subplot(gs[1, 1])
ax6 = fig.add_subplot(gs[1, 2])
middle_axes = [ax4, ax5, ax6]

# Create the axis for the bottom row (1 plot)
ax7 = fig.add_subplot(gs[2, 0])
bottom_axes = [ax7]

bar_width = 1

# Loop through the first 3 n values for the top row
for i, n_val in enumerate(n_values[:3]):
    ax = top_axes[i]
    df_n = df[df['n'] == n_val].iloc[0]
    time_taken = [df_n[tt] for tt in tree_types]
    x_pos = np.arange(len(tree_types))
    ax.bar(x_pos, time_taken, color=colors, edgecolor='none', width=bar_width)
    ax.set_title(f"n = {n_val}")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(tree_types, rotation=45, ha='right')
    ax.set_xlabel("")
    ax.set_ylabel("Time Taken (s)")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', linestyle='--', alpha=0.7)

# Loop through the next 3 n values for the middle row
for i, n_val in enumerate(n_values[3:6]):
    ax = middle_axes[i]
    df_n = df[df['n'] == n_val].iloc[0]
    time_taken = [df_n[tt] for tt in tree_types]
    x_pos = np.arange(len(tree_types))
    ax.bar(x_pos, time_taken, color=colors, edgecolor='none', width=bar_width)
    ax.set_title(f"n = {n_val}")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(tree_types, rotation=45, ha='right')
    ax.set_xlabel("")
    ax.set_ylabel("Time Taken (s)")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', linestyle='--', alpha=0.7)

# Plot the last n value in the bottom row
ax = bottom_axes[0]
df_n = df[df['n'] == n_values[6]].iloc[0]
time_taken = [df_n[tt] for tt in tree_types]
x_pos = np.arange(len(tree_types))
ax.bar(x_pos, time_taken, color=colors, edgecolor='none', width=bar_width)
ax.set_title(f"n = {n_values[6]}")
ax.set_xticks(x_pos)
ax.set_xticklabels(tree_types, rotation=45, ha='right')
ax.set_xlabel("")
ax.set_ylabel("Time Taken (s)")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='y', linestyle='--', alpha=0.7)


# Create a new axis in the empty space for the legend
ax_legend = fig.add_subplot(gs[2, 1:]) # Using remaining space in the bottom row
ax_legend.axis('off')

# Add the legend to the new axis
handles = [plt.Rectangle((0,0),1,1, color=color) for color in colors]
ax_legend.legend(handles, tree_types, loc='center', title="Tree Type")

# Adjust layout to prevent overlapping titles and labels
plt.tight_layout()

# Save the figure to the PDF file
pdf = PdfPages('n_sweep_plots_legend_fixed.pdf')
pdf.savefig(fig)
pdf.close()
plt.close(fig)

# Save a separate PNG image for direct display
fig.savefig('n_sweep_plots_legend_fixed.png')