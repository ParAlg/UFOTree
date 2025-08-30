import pandas as pd
import matplotlib.pyplot as plt

# Load the data from the CSV file
df = pd.read_csv('diameter_sweep_parellel_update_10000000_1000000_final.csv')

# Plot the data
plt.figure(figsize=(10, 6))

colors = [
        '#8769b6',
        '#b6f486',
        '#88d4c3',
        '#FF609A'
        ]
# Plot each line with a marker
plt.plot(df['Test Case'], df['ETT (Skip List)'], label='ETT (Skip List)', marker='o', color=colors[0])
plt.plot(df['Test Case'], df['UFO Tree'], label='UFO Tree', marker='o', color=colors[1])
plt.plot(df['Test Case'], df['Topology Tree'], label='Topology Tree', marker='o', color=colors[2])
plt.plot(df['Test Case'], df['Rake-Compress Tree'], label='Rake-Compress Tree', marker='o', color=colors[3])

# Add a horizontal line at y=1 to show the 1-second cutoff
plt.axhline(y=1, color='r', linestyle='--', label='1-second cutoff')

# Set y-axis limits
plt.ylim(0, 25)

# Add labels and title
plt.xlabel('Test Case')
plt.ylabel('Time Taken')
plt.title('Line Plot with Data Points and 1-second Cutoff')
plt.legend()
plt.grid(True)

# Save the plot to a PDF file
plt.savefig('diameter_sweep_parallel_update_line_graph_with_points_and_cutoff.pdf')

print("Plot saved to diameter_sweep_parallel_update_line_graph_with_points_and_cutoff.pdf")
