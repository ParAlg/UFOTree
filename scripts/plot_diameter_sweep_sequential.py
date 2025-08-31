import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# File names
file_paths = {
    'Updates': '../results/diameter_sweep_update_10000000_final.csv',
    'Path Queries': '../results/diameter_sweep_path_query_1000000_final.csv',
    'Connectivity Queries': '../results/diameter_sweep_conn_query_1000000_final.csv'
}

# Read data from files
dataframes = {key: pd.read_csv(path) for key, path in file_paths.items()}

# Define custom colors, styles, and markers
plot_styles = {
    'all': {
        'colors': ['#000000', '#b6f486', '#006400', '#4169E1', '#400e63', '#8769b6', '#88d4c3', '#FF609A'],
        'linestyles': ['--', '-', '-', '-', '-', '-.', '--', (0, (3, 5, 1, 5, 1, 5))],
        'markers': ['o', '^', 's', 'D', 'x', '*', '^', 'p']
    },
    'path_queries': {
        'colors': ['#000000', '#b6f486', '#006400', '#88d4c3', '#FF609A'],
        'linestyles': ['--', '-', '-', '--', (0, (3, 5, 1, 5, 1, 5))],
        'markers': ['o', '^', 's', '^', 'p']
    }
}

# Set global plot parameters
plt.rc('font', family='serif')
plt.rc('legend', fontsize=25)

# Create figure and subplots
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(25, 7))
plt.subplots_adjust(wspace=0.2, left=0.06, right=0.98, top=0.75, bottom=0.15)

# Map subplot titles to their data and styles
plot_configs = {
    '(a) Updates': {'df': dataframes['Updates'], 'styles': plot_styles['all']},
    '(b) Connectivity Queries': {'df': dataframes['Connectivity Queries'], 'styles': plot_styles['all']},
    '(c) Path Queries': {'df': dataframes['Path Queries'], 'styles': plot_styles['path_queries']}
}

# Loop through subplots and plot data
for i, (title, config) in enumerate(plot_configs.items()):
    df = config['df']
    styles = config['styles']
    ax = axes[i]
    
    for j, col in enumerate(df.columns[1:]):
        ax.plot(df['Alpha'], df[col], 
                marker=styles['markers'][j],
                linestyle=styles['linestyles'][j],
                color=styles['colors'][j],
                markersize=20,
                linewidth=9,
                label=col)

    ax.set_title(title, fontsize=30)
    ax.set_yscale('log')
    ax.grid(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel(r'$\alpha$', fontsize=30)
    ax.set_ylabel('Time (s)', fontsize=35)
    ax.tick_params(axis='both', which='major', labelsize=25)

# Adjust layout and save the figure
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.03), ncol=4, frameon=False, columnspacing=1.5)

plt.savefig('../results/diameter_sweep_sequential.pdf')
print("Plot saved to results/diameter_sweep_sequential.pdf")