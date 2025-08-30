import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the data from the CSV file
df = pd.read_csv("n_sweep_parellel_update_final.csv")

# Use a professional-looking style
plt.style.use('seaborn-v0_8-notebook')

# Use a serif font for the entire plot
plt.rc('font', family='serif')

# Create a single figure and axis for the line plot
fig, ax = plt.subplots(figsize=(10, 8))

# Define the tree types and markers
tree_types = ['Path', 'Binary', '64ary', 'Star']
markers = ['o', 's', 'D', '^'] # Circle, Square, Diamond, Triangle

# Plot each tree type as a line with markers
for i, tree_type in enumerate(tree_types):
    ax.plot(df['n'], df[tree_type], label=tree_type, marker=markers[i], markersize=12)

# Set the plot title and labels
#ax.set_title("Time Taken vs. Number of Updates (n)", fontsize=18)
ax.set_xlabel("Tree size (n)", fontsize=16)
ax.set_ylabel("Time Taken (s)", fontsize=16)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Improve readability
ax.grid(True, linestyle='--', alpha=0.7)
ax.set_xscale('log') # Use a logarithmic scale for the x-axis for better visualization
ax.set_yscale('log') # Set the y-axis to a logarithmic scale as well
ax.tick_params(axis='both', which='major', labelsize=20)

# Add a legend
ax.legend(title="Tree Type", fontsize=12, title_fontsize=14)

# Adjust layout and save the figure
plt.tight_layout()
fig.savefig('../results/n_sweep_log_log_plot.pdf')
print("Plot saved to results/n_sweep_log_log_plot.pdf")