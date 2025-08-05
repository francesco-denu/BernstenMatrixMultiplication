import numpy as np

def calculate_stats(filename):
    with open(filename, 'r') as file:
        content = file.read()

    sets = content.strip().split('\n\n')
    results = []

    for data_set in sets:
        values = list(map(float, data_set.split()))
        avg = np.mean(values)
        var = np.var(values)
        results.append((avg, var))

    return results

def save_stats(output_filename, all_results):
    with open(output_filename, 'w') as file:
        for file_index, results in enumerate(all_results, 1):
            file.write(f"Results for file {file_index}:\n")
            for i, (avg, var) in enumerate(results, 1):
                file.write(f"  Set {i}:\n")
                file.write(f"    Average: {avg}\n")
                file.write(f"    Variance: {var}\n")
            file.write("\n")

def main():
    filenames = ['sequential_1', 'bernsten_8_cores_1', 'cannon_9_cores_1']
    output_filename = '_final_BENCHMARKS.txt'
    all_results = []

    for filename in filenames:
        results = calculate_stats(filename)
        all_results.append(results)

    save_stats(output_filename, all_results)

if __name__ == '__main__':
    main()


"""import numpy as np
import matplotlib.pyplot as plt

def calculate_stats(filename):
    with open(filename, 'r') as file:
        content = file.read()

    sets = content.strip().split('\n\n')
    results = []

    for data_set in sets:
        values = list(map(float, data_set.split()))
        avg = np.mean(values)
        var = np.var(values)
        results.append((avg, var))

    return results

def plot_stats(all_results):
    colors = ['r', 'g', 'b']  # Different colors for each file

    plt.figure(figsize=(10, 6))

    for file_index, results in enumerate(all_results, 1):
        color = colors[file_index % len(colors)]
        for set_index, (avg, var) in enumerate(results):
            x_position = set_index + 1  # X position for the set
            plt.errorbar(x_position, avg, yerr=np.sqrt(var), fmt='o', color=color, 
                         ecolor=color, elinewidth=2, capsize=5, label=f'File {file_index}' if set_index == 0 else "")

    plt.xlabel('Set Index')
    plt.ylabel('Values')
    plt.title('Average Values with Variance by Set')
    plt.legend()
    plt.show()

def main():
    filenames = ['_BERNSTEN_1_CORE.txt', '_BERNSTEN_8_CORES.txt', '_CANNON_9_CORES.txt']
    all_results = []

    for filename in filenames:
        results = calculate_stats(filename)
        all_results.append(results)

    plot_stats(all_results)

if __name__ == '__main__':
    main()




"""