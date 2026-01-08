import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
from matplotlib.ticker import FuncFormatter, MaxNLocator

def read_and_process_csv(filepath: str) -> tuple:
    try:
        with open(filepath, mode='r', newline='') as file:
            reader = csv.reader(file)
            header = next(reader)
            bar_groups = header[1:]
            data_dict = {group: [] for group in bar_groups}
            category_labels = []
            for row in reader:
                test_case_name = row[0]
                category_labels.append(test_case_name)
                for i, value in enumerate(row[1:]):
                    data_dict[bar_groups[i]].append(float(value))
            return data_dict, category_labels, bar_groups
    except Exception as e:
        print(f"Error: {e}")
        return None, None, None

def format_scientific(x, pos):
    if x == 0:
        return "0"
    exponent = int(np.floor(np.log10(abs(x))))
    coeff = x / 10**exponent
    return f"${coeff:g} \\cdot 10^{{{exponent}}}$"

def plot_multi_bar_chart(data: dict, categories: list, output_filepath: str):
    if not data or not categories:
        return
    
    bar_group_names = list(data.keys())
    num_bar_groups = len(bar_group_names)
    colors = ['#000000', '#b6f486', '#006400', '#4169E1', '#400e63', '#8769b6', '#88d4c3', '#FF609A']
    bar_width = 0.11
    
    plt.rc('font', family='serif', size=120)
    plt.rc('axes', titlesize=144, labelsize=120)
    plt.rc('xtick', labelsize=100)
    plt.rc('ytick', labelsize=80)
    plt.rc('legend', fontsize=90)

    fig, ax = plt.subplots(1, 1, figsize=(84, 16))

    x_positions = np.arange(len(categories))
    all_values = []
    for i, group_name in enumerate(bar_group_names):
        offset = bar_width * i
        values = data[group_name]
        all_values.extend(values)
        ax.bar(x_positions + offset, values, bar_width, label=group_name, color=colors[i % len(colors)])

    ax.set_ylabel("Memory(B)") 
    ax.set_yscale('linear') 
    
    # Bold Spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(5)
    ax.spines['left'].set_linewidth(5)
    ax.tick_params(axis='both', which='major', width=4, length=15, pad=20)

    if all_values:
        max_val = max(all_values)
        # Ensure headroom for legend
        ax.set_ylim(0, max_val * 1.25)
        
        # --- SET TO EXACTLY 4 INTERVALS ---
        # nbins=4 creates 4 gaps between 5 labels
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 5, 10]))
        ax.yaxis.set_major_formatter(FuncFormatter(format_scientific))

    ax.set_xticks(x_positions + bar_width * (num_bar_groups - 1) / 2)
    ax.set_xticklabels(categories, ha="center", rotation=0)

    ax.grid(axis='y', linestyle='--', color='#cccccc', linewidth=2, alpha=0.8)
    ax.set_axisbelow(True)

    plt.tight_layout(pad=1.0)
    plt.subplots_adjust(top=0.70, bottom=0.15)

    fig.legend(bar_group_names, loc='upper center', bbox_to_anchor=(0.5, 1.0), 
               ncol=4, frameon=False, columnspacing=1.5)

    buffer = 0.2
    ax.set_xlim(x_positions[0] - buffer, x_positions[-1] + num_bar_groups * bar_width + buffer)

    plt.savefig(output_filepath)
    print(f"Plot saved to {output_filepath}")
    plt.close()

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script.py <input.csv> <output_file>")
        sys.exit(1)
    d, l, g = read_and_process_csv(sys.argv[1])
    if d:
        plot_multi_bar_chart(d, l, sys.argv[2])
