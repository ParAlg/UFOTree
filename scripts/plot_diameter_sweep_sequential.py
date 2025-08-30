import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler

# File names
file1 = '../results/diameter_sweep_update_10000000_final.csv'
file2 = '../results/diameter_sweep_path_query_1000000_final.csv'
file3 = '../results/diameter_sweep_conn_query_1000000_final.csv'

# Read the CSV files
df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)
df3 = pd.read_csv(file3)
    
lcTree = '#000000'
ufoT = '#b6f486'
splayTT = '#006400'
ettSplay = '#4169E1'
ettTreap = '#400e63'
ettSkip = '#8769b6'
topologyT = '#88d4c3'
rcTree = '#FF609A'
# Define custom color lists for each plot
# The number of colors must match the number of lines (columns to be plotted)
colors1 = [lcTree, ufoT, splayTT, ettSplay, ettTreap, ettSkip, topologyT, rcTree]
colors2 = [lcTree, ufoT, splayTT, topologyT , rcTree]
colors3 = colors1 

# Define custom line styles and markers for the plots
linestyles1 = ['--', '-', '-', '-', '-', '-.', '--', (0, (3, 5, 1, 5, 1, 5))]
markers1 = ['o', '^', 's', 'D', 'x', '*', '^', 'p']

linestyles2 = ['--', '-', '-', '--', (0, (3, 5, 1, 5, 1, 5))]
markers2 = ['o', '^', 's', '^', 'p']

linestyles3 = linestyles1 
markers3 = markers1 

#Set font for the figure
plt.rc('font', family='serif')
plt.rc('legend', fontsize=25)
# Create a figure with 3 subplots in a row and reduced height
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(40,9))
plt.subplots_adjust(wspace=0.5)

# Plot for the first CSV
for i, col in enumerate(df1.columns[1:]):
    axes[0].plot(df1['Alpha'], df1[col], marker=markers1[i], linestyle=linestyles1[i], color=colors1[i], markersize=12, label=col)

axes[0].set_title('Updates', fontsize=20)
axes[0].set_yscale('log')
axes[0].grid(True)
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)

# Plot for the second CSV
for i, col in enumerate(df2.columns[1:]):
    axes[1].plot(df2['Alpha'], df2[col], marker=markers2[i], linestyle=linestyles2[i], color=colors2[i], markersize=12, label=col)
axes[1].set_title('Path Queries',fontsize=20)
axes[1].set_yscale('log')
axes[1].grid(True)
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)

# Plot for the third CSV
for i, col in enumerate(df3.columns[1:]):
    axes[2].plot(df3['Alpha'], df3[col], marker=markers3[i], linestyle=linestyles3[i], color=colors3[i], markersize=12, label=col)

axes[2].set_title('Connectivity Queries', fontsize=20)
axes[2].set_yscale('log')
axes[2].grid(True)
axes[2].spines['top'].set_visible(False)
axes[2].spines['right'].set_visible(False)

for a in axes:
    a.set_xlabel(r'$\alpha$', fontsize=25)
    a.set_ylabel('Time Taken (s)', fontsize=25)
    a.tick_params(axis='both', which='major', labelsize=25)

# Adjust layout and save the figure
handles, labels = axes[2].get_legend_handles_labels()
#fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0), ncol=len(labels), frameon=False)
fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1), ncol=4, frameon=False, columnspacing=1.5)

#plt.tight_layout()
plt.savefig('../results/diameter_sweep_sequential.pdf')
print("Plot saved to results/diameter_sweep_sequential.pdf")