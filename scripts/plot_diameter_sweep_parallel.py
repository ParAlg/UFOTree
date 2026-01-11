import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

# Force TrueType fonts (Type 42) instead of Type 3
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Load the data
df = pd.read_csv('../results/diameter_sweep_parellel_update_10000000_1000000_final.csv')

# FIX: If Pandas treated the first column as an index (The Shift)
if df.index.name is None and df.iloc[:, -1].isna().all():
    # 1. Pull the index (0.00, 0.25...) back into the dataframe as a column
    df = df.reset_index()
    
    # 2. The shift caused the last column to be NaN (garbage). Drop it.
    df = df.iloc[:, :-1]
    
    # 3. Re-assign the correct column headers
    df.columns = ['Test Case', 'ETT (Skip List)', 'UFO Tree', 'Topology Tree', 'Rake-Compress Tree']

print(df.head()) # Verify the fix

# ... (Previous data loading code) ...

# Plot the data
plt.figure(figsize=(10, 7)) # Increased height slightly to accommodate the top legend

colors = ['#8769b6', '#b6f486', '#88d4c3', '#FF609A']

# Plot each line
plt.plot(df['Test Case'], df['ETT (Skip List)'], label='ETT (Skip List)', marker='o', color=colors[0], markersize=12)
plt.plot(df['Test Case'], df['UFO Tree'], label='UFO Tree', marker='o', color=colors[1], markersize=12)
plt.plot(df['Test Case'], df['Topology Tree'], label='Topology Tree', marker='o', color=colors[2], markersize=12)
plt.plot(df['Test Case'], df['Rake-Compress Tree'], label='Rake-Compress Tree', marker='o', color=colors[3], markersize=12)

# Axis setup
plt.xscale("linear")
plt.yscale('log')

ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=20)
ax.set_xlabel(r'$\alpha$', fontsize=20)   # Increased to match ticks
ax.set_ylabel('Time Taken (s)', fontsize=20)

plt.grid(True)

# --- NEW LEGEND SETUP ---
# bbox_to_anchor=(x, y): (0.5, 1.02) places it in the middle x, just above the top y
# loc='lower center': The anchor point refers to the bottom-center of the legend box
# ncol=2: Arranges the 4 items into 2 rows of 2 (cleaner than one long line)
plt.legend(loc='lower center', bbox_to_anchor=(0.5, 1.02), 
           ncol=2, fontsize=16, frameon=False)

# Adjust layout to make sure the external legend isn't cut off
plt.tight_layout()

# Save
plt.savefig('../results/diameter_sweep_parallel_updates.pdf', bbox_inches='tight')
print("Plot saved.")
