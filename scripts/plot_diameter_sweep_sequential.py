import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sys # Import the sys module to handle command-line arguments

# --- Argument Parsing ---
# The script expects 4 arguments: 3 input CSV files and 1 output PDF path
if len(sys.argv) != 5:
    print("Usage: python3 plot_diameter_sweep.py <updates_csv> <path_queries_csv> <connectivity_queries_csv> <output_pdf_path>")
    sys.exit(1)

# Assign command-line arguments to descriptive variables
updates_path = sys.argv[1]
path_queries_path = sys.argv[2]
conn_queries_path = sys.argv[3]
output_pdf_path = sys.argv[4]

# File names are now taken from command-line arguments
file_paths = {
    'Updates': updates_path,
    'Path Queries': path_queries_path,
    'Connectivity Queries': conn_queries_path
}

# Read data from files
try:
    dataframes = {key: pd.read_csv(path) for key, path in file_paths.items()}
except FileNotFoundError as e:
    print(f"Error: One of the input files was not found: {e}")
    sys.exit(1)
except pd.errors.EmptyDataError as e:
    print(f"Error: One of the input files is empty: {e}")
    sys.exit(1)


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

# Use the command-line argument for the output path
plt.savefig(output_pdf_path)
print(f"Plot saved to {output_pdf_path}")