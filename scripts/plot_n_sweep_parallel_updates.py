import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

# Force TrueType fonts (Type 42) instead of Type 3
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Load the data from the CSV file
df = pd.read_csv("../results/n_sweep_parallel_update.csv")

# Use a professional-looking style
plt.style.use('seaborn-v0_8-notebook')

# Use a serif font for the entire plot
plt.rc('font', family='serif')

# Create a single figure and axis for the line plot
fig, ax = plt.subplots(figsize=(20, 9))

# Define the tree types and markers by reading the column names from the CSV
tree_types = df.columns.tolist()[1:]
markers = ['o', 's', 'D', '^'] # Circle, Square, Diamond, Triangle

# Plot each tree type as a line with markers
for i, tree_type in enumerate(tree_types):
    ax.plot(df['n'], df[tree_type], label=tree_type, marker=markers[i], markersize=40, linewidth=15)

# Set the plot title and labels
#ax.set_title("Time Taken vs. Number of Updates (n)", fontsize=18)
ax.set_xlabel("Input Size (n)", fontsize=50)
ax.set_ylabel("Time(s)", fontsize=50)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Improve readability
ax.grid(True, linestyle='--', alpha=0.7)
ax.set_xscale('log') # Use a logarithmic scale for the x-axis for better visualization
ax.set_yscale('log') # Set the y-axis to a logarithmic scale as well
ax.tick_params(axis='both', which='major', labelsize=40)

# Add a legend on top, centered on the figure
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, fontsize=44, loc='upper center', bbox_to_anchor=(0.5, 1.04), ncol=4, frameon=False)

# Adjust layout and save the figure
plt.subplots_adjust(left=0.12, right=0.99, top=0.88, bottom=0.20)
fig.savefig('../results/n_sweep.pdf')
print("Plot saved to results/n_sweep.pdf")
