import matplotlib as mpl
import re
from collections import defaultdict

def run():
    # Put data into dictionary
    data = open("nums.txt", encoding='utf-8')
    benchmarks = []
    test_type, n = '', ''
    numbers = defaultdict(list) 
    data_structure = ""
    for line in data:
        match = re.search(r"\[ RUNNING (.+)? BENCHMARK WITH n=([0-9]+) \]", line) 
        if match:
            test_type = match.group(1)
            n = match.group(2)
            numbers = defaultdict(list)
            benchmarks.append((test_type, n, numbers))
            continue
        match = re.search(r"Init Space: (.+) Bytes|Total Time: (.+)|Peak Space: (.+) Bytes", line) 
        if match:
            if match.group(1):
                numbers[data_structure].append(match.group(1))
            elif match.group(2):
                numbers[data_structure].append(match.group(2))
            else:
                numbers[data_structure].append(match.group(3))
        else:
            match = re.search(r"([A-Z|a-z|0-9]+)(?s:.)*", line)
            data_structure = match.group(1) if match else None
    
    plot_benchmark_data(benchmarks)

def plot_benchmark_data(benchmarks):
    return 0

run()
