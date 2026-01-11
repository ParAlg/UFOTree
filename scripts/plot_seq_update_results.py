# Import the necessary libraries
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import csv
import sys

# Force TrueType fonts (Type 42) instead of Type 3
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# A complete, self-contained script to create a multi-bar plot.

def read_and_process_csv(filepath: str) -> tuple:
    """
    Reads a CSV file and processes the data for plotting a multi-bar chart.

    Args:
        filepath (str): The path to the CSV file.

    Returns:
        tuple: A tuple containing a dictionary of bar data and a list of
               category labels. Returns (None, None) if the file is not found
               or an error occurs.
    """
    try:
        with open(filepath, mode='r', newline='') as file:
            reader = csv.reader(file)
            header = next(reader)  # Read the first row (data structures)
            bar_groups = header[1:] # Store data structures as bar groups

            # Dictionary to hold the data for each bar group (data structure)
            data_dict = {group: [] for group in bar_groups}
            
            # List to hold the x-axis labels (test case names)
            category_labels = []

            # Process each row after the header
            for row in reader:
                test_case_name = row[0]
                category_labels.append(test_case_name)
                
                # Append the values for each data structure
                for i, value in enumerate(row[1:]):
                    data_dict[bar_groups[i]].append(float(value))
                    
            return data_dict, category_labels, bar_groups

    except FileNotFoundError:
        print(f"Error: The file at {filepath} was not found.")
        return None, None, None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None, None, None


def plot_multi_bar_chart(data: dict, categories: list, title: str, xlabel: str, ylabel: str, output_filepath: str):
    """
    Creates and displays a multi-bar plot from the given data, split into two rows.

    Args:
        data (dict): A dictionary where keys are the bar group labels (e.g., 'Array')
                     and values are a list of data points for each category.
        categories (list): A list of category labels for the x-axis.
        title (str): The title of the plot.
        xlabel (str): The label for the x-axis.
        ylabel (str): The label for the y-axis.
        output_filepath (str): The path to save the generated plot image.
    """
    # Check if the input data is valid
    if not data or not categories:
        print("Error: The input data is empty.")
        return
    
    # Get the names of the bar groups
    bar_group_names = list(data.keys())
    num_bar_groups = len(bar_group_names)
    
    # Define a custom color palette based on your previous scheme, with reordered colors
    # to increase contrast between successive bars.
    colors = [
        '#000000',
        '#b6f486',
        '#006400',
        '#4169E1',
        '#400e63',
        '#8769b6',
        '#88d4c3',
        '#FF609A'
    ]

    # Split the categories and data into two halves
    half_categories = len(categories) // 2
    
    categories_part1 = categories[:half_categories]
    categories_part2 = categories[half_categories:]
    
    data_part1 = {group: values[:half_categories] for group, values in data.items()}
    data_part2 = {group: values[half_categories:] for group, values in data.items()}
    
    # Set the width of each individual bar
    bar_width = 0.11
    
    # Use a serif font for the entire plot
    plt.rc('font', family='serif')
    
    # Adjust font sizes
    plt.rc('font', size=120)          # controls default text sizes
    plt.rc('axes', titlesize=144)     # fontsize of the axes title
    plt.rc('axes', labelsize=120)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=100)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=80)    # fontsize of the tick labels (1.5x smaller than 120)
    plt.rc('legend', fontsize=90)    # legend fontsize (1.5x bigger than previous value)
    plt.rc('figure', titlesize=144)  # fontsize of the figure title

    # Create the figure and two subplots in two rows
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(84, 24))

    # Plot the first half of the groups on the top subplot (ax1)
    x_positions1 = np.arange(len(categories_part1))
    for i, group_name in enumerate(bar_group_names):
        offset = bar_width * i
        ax1.bar(x_positions1 + offset, data_part1[group_name], bar_width, label=group_name, color=colors[i % len(colors)])

    # Plot the second half of the groups on the bottom subplot (ax2)
    x_positions2 = np.arange(len(categories_part2))
    for i, group_name in enumerate(bar_group_names):
        offset = bar_width * i
        ax2.bar(x_positions2 + offset, data_part2[group_name], bar_width, label=group_name, color=colors[i % len(colors)])

    # Configure the first subplot (top row)
    ax1.set_ylabel("Time(s)")
    ax1.set_yscale('log')
    ax1.set_xticks(x_positions1 + bar_width * (num_bar_groups - 1) / 2)
    ax1.set_xticklabels(categories_part1, ha="center", rotation=0)
    ax1.grid(axis='y', linestyle='--', alpha=0.7)
    # Add extra space between the x-axis and the labels
    ax1.tick_params(axis='x', which='major', pad=20)
    # Add back the border for the left and bottom spines
    ax1.spines['bottom'].set_visible(True)
    ax1.spines['bottom'].set_linewidth(2)
    ax1.spines['left'].set_visible(True)
    ax1.spines['left'].set_linewidth(2)
    # Ensure top and right spines remain invisible
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Configure the second subplot (bottom row)
    ax2.set_ylabel("Time(s)")
    ax2.set_yscale('log')
    ax2.set_xticks(x_positions2 + bar_width * (num_bar_groups - 1) / 2)
    ax2.set_xticklabels(categories_part2, ha="center", rotation=0)
    ax2.grid(axis='y', linestyle='--', alpha=0.7)
    # Add extra space between the x-axis and the labels
    ax2.tick_params(axis='x', which='major', pad=20)
    # Add back the border for the left and bottom spines
    ax2.spines['bottom'].set_visible(True)
    ax2.spines['bottom'].set_linewidth(2)
    ax2.spines['left'].set_visible(True)
    ax2.spines['left'].set_linewidth(2)
    # Ensure top and right spines remain invisible
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    # Use tight_layout to handle subplot spacing
    plt.tight_layout(pad=1.0, h_pad=0.5)
    
    # Manually adjust the top and bottom of the plot to control space
    plt.subplots_adjust(top=0.8, bottom=0.075)

    # Place the legend at the top of the figure and arrange entries horizontally
    fig.legend(bar_group_names, loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=4, frameon=False, columnspacing=1.5)

    # Reduce space around the edges by setting explicit x-axis limits
    # This centers the plot on the bar groups without excess space
    buffer = 0.2
    ax1.set_xlim(x_positions1[0] - buffer, x_positions1[-1] + num_bar_groups * bar_width + buffer)
    ax2.set_xlim(x_positions2[0] - buffer, x_positions2[-1] + num_bar_groups * bar_width + buffer)

    # Save the plot to a file
    plt.savefig(output_filepath)
    print(f"Plot saved to {output_filepath}")

    # Close the plot to free up memory
    plt.close()

if __name__ == '__main__':
    # Check if the correct number of command-line arguments were provided
    if len(sys.argv) != 4:
        print("Usage: python script_name.py <input_csv_file1> <input_csv_file2> <output_pdf_file>")
        sys.exit(1)

    # Get the file paths from the command-line arguments
    input_file_path1 = sys.argv[1]
    input_file_path2 = sys.argv[2]
    output_file_path = sys.argv[3]
    
    # Read the data from the two CSV files
    data1, labels1, groups1 = read_and_process_csv(input_file_path1)
    data2, labels2, groups2 = read_and_process_csv(input_file_path2)

    # Exit if either file could not be read or if they have different column headers
    if not data1 or not data2:
        print("Error: Could not read one or both input files. Exiting.")
        sys.exit(1)
    
    if groups1 != groups2:
        print("Error: The two CSV files must have the same data structure columns.")
        sys.exit(1)

    # Combine the data from both files
    combined_data = {group: data1[group] + data2[group] for group in groups1}
    combined_labels = labels1 + labels2

    # Plot the combined data
    plot_multi_bar_chart(
        data=combined_data,
        categories=combined_labels,
        title='Performance Comparison',
        xlabel='Test Cases',
        ylabel='Performance Metric',
        output_filepath=output_file_path
    )
