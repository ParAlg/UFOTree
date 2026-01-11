import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import textwrap

# Force TrueType fonts (Type 42) instead of Type 3
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Set default font size and font family for all text elements
# Removed specific 'Times New Roman' as it's not found, will use Matplotlib's default serif font
plt.rcParams.update({
    'font.size': 30,
    'font.family': 'serif',  # Set font family to generic serif
    # 'font.serif': ['Times New Roman'] # Removed this line
}) 

# Check if the file path is provided
if len(sys.argv) < 3:
    print("Usage: python3 plot_results.py <input_csv_path> <output_pdf_path>")
    sys.exit(1)

# Get the file path from the first two command-line arguments
input_path = sys.argv[1]
output_path = sys.argv[2]

# Load the CSV file and take only the first 9 rows
data = pd.read_csv(input_path).head(9)

# Drop the extra column if it's the last one and contains only NaN in the rest of the rows
data = data.loc[:, ~data.columns.str.contains('Unnamed')]  # Drop columns with 'Unnamed'

# Extract data
categories = data.iloc[:, 0]  # First column for x-axis labels
# Convert all other columns to numeric, coercing errors to NaN
values = data.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')

# Number of categories and number of groups
n_categories = len(categories)
n_groups = len(values.columns)

# Bar width and positions - Halved from 0.2 to 0.1
bar_width = 0.1
index = np.arange(n_categories)

# Wrap long category labels to force two lines (adjust width as needed based on label content)
wrapped_labels = [textwrap.fill(label, width=10) for label in categories]

# Create the plot - Increased figsize width from 12 to 24
fig, ax = plt.subplots(figsize=(24, 12)) 

# Set y-axis to log scale
# Note: Log scales cannot display zero or negative values.
# If your data contains zeros, they will not be plotted or may cause errors.
ax.set_yscale('log')

# Loop through each column and create a bar for it
for i in range(n_groups):
    plt.bar(index + i * bar_width, values.iloc[:, i], bar_width, label=values.columns[i])

def format_title(filename):
    filename = os.path.basename(filename)
    base_name = filename.replace('.csv', '')
    parts = base_name.split('_')
    formatted_string = f"{parts[0].capitalize()} {parts[1].capitalize()}, n={parts[2]}"
    return formatted_string

# Add labels, title, and legend
title = format_title(input_path)
plt.xlabel('Input Case') 
if title.startswith("Peak Space"):
    plt.ylabel('Space (Bytes)') 
else: # Assuming 'Time (Seconds)'
    plt.ylabel('Time (Seconds)') 
    max_y_value = values.values.max()

    if not np.isnan(max_y_value):
        y_tick_interval = 0.5
    else:
        print("Warning: Max y-value is NaN. Y-axis ticks will be auto-determined.")

plt.title(title)
plt.xticks(index + bar_width * (n_groups - 1) / 2, wrapped_labels, rotation=0) 
plt.legend(fontsize=24) # Set legend font size explicitly to 24

# Automatically adjust subplot params for a tight layout
fig.tight_layout()

# Save the plot as a PDF file, ensuring content is not cut off
plt.savefig(output_path, format='pdf', bbox_inches='tight')
