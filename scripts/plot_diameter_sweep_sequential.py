import pandas as pd
import matplotlib.pyplot as plt

# File names
file1 = '../results/diameter_sweep_update_10000000_final.csv'
file2 = '../results/diameter_sweep_path_query_1000000_final.csv'
file3 = '../results/diameter_sweep_conn_query_1000000_final.csv'

# Read the CSV files
df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)
df3 = pd.read_csv(file3)
    
lcTree = '#000000'
splayTT = '#b6f486'
ettSplay = '#006400'
topologyTree = '#4169E1'
ufoT = '#400e63'
ettTreap = '#8769b6'
ettSkip = '#88d4c3'
rcTree = '#FF609A'
# Define custom color lists for each plot
# The number of colors must match the number of lines (columns to be plotted)
colors1 = [lcTree, splayTT, ettSplay, topologyTree, ufoT, ettTreap, ettSkip, rcTree]
colors2 = [lcTree, ufoT, splayTT, topologyTree, rcTree]
colors3 = colors1 

#Set font for the figure
plt.rc('font', family='serif')
# Create a figure with 3 subplots in a row
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(20, 7))


# Plot for the first CSV
df1.plot(x='Alpha', ax=axes[0], marker='o', color =colors1, markersize=12, legend=False)
axes[0].set_title('Update Queries $n=10^7$')
axes[0].set_yscale('log')
axes[0].grid(True)
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)

# Plot for the second CSV
df2.plot(x='Alpha', ax=axes[1], marker='o', color = colors2, markersize=12,legend=False)
axes[1].set_title('Path Queries $n=10^6$')
axes[1].set_yscale('log')
axes[1].grid(True)
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)

# Plot for the third CSV
df3.plot(x='Alpha', ax=axes[2], marker='o', color=colors3, markersize=12,legend=True)
axes[2].set_title('Connectivity Queries $n=10^6$')
axes[2].set_yscale('log')
axes[2].grid(True)
axes[2].spines['top'].set_visible(False)
axes[2].spines['right'].set_visible(False)

for a in axes:
    a.set_xlabel(r'$\alpha$', fontsize=20)
    a.set_ylabel('Time Taken (s)', fontsize=20)
    a.tick_params(axis='both', which='major', labelsize=20)

# Adjust layout and save the figure
#handles, labels = axes[0].get_legend_handles_labels()
#fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0), ncol=len(labels), frameon=False)

plt.tight_layout()
plt.savefig('../results/diameter_sweep_sequential.pdf')
print("Plot saved to results/diameter_sweep_sequential.pdf")