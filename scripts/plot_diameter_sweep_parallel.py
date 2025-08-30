import pandas as pd
import matplotlib.pyplot as plt

# Load the data from the CSV file
df = pd.read_csv('../results/diameter_sweep_parellel_update_10000000_1000000_final.csv')

# Plot the data
plt.figure(figsize=(10, 6))

colors = [
        '#8769b6',
        '#b6f486',
        '#88d4c3',
        '#FF609A'
        ]
# Plot each line with a marker
plt.plot(df['Test Case'], df['ETT (Skip List)'], label='ETT (Skip List)', marker='o', color=colors[0], markersize=12)
plt.plot(df['Test Case'], df['UFO Tree'], label='UFO Tree', marker='o', color=colors[1], markersize=12)
plt.plot(df['Test Case'], df['Topology Tree'], label='Topology Tree', marker='o', color=colors[2], markersize=12)
plt.plot(df['Test Case'], df['Rake-Compress Tree'], label='Rake-Compress Tree', marker='o', color=colors[3], markersize=12)

# Add a horizontal line at y=1 to show the 1-second cutoff
#plt.axhline(y=1, color='r', linestyle='--', label='1-second cutoff')

# Set y-axis limits

# Add labels and title
plt.xscale("linear")
plt.yscale('log')

ax = plt.gca()

# Disable the top and right spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=20)
ax.set_xlabel(r'$\alpha$', fontsize=16)
ax.set_ylabel('Time Taken (s)', fontsize=16)

plt.legend()
plt.grid(True)

plt.tight_layout()
# Save the plot to a PDF file
plt.savefig('../results/diameter_sweep_parallel_updates.pdf')

print("Plot saved to results/diameter_sweep_parallel_updates.pdf")
