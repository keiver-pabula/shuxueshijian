import matplotlib.pyplot as plt
import numpy as np
import os
import glob

# Read spline file
def read_spline(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    intervals = []
    coefficients = []

    # Special case for Point 2 files
    if "P2_s23" in file_path:
        i = 0
        while i < len(lines):
            # Extract interval
            interval = list(map(float, lines[i].strip().split(",")))
            intervals.append(interval)
            i += 1

            # Extract coefficients
            coeff = list(map(float, lines[i].strip().split(" ")))
            coefficients.append(coeff)
            i += 1
    else:
        for line in lines:
            if line.startswith("Interval:"):
                # Extract interval details
                interval = list(map(float, line.split("[")[1].split("]")[0].split(",")))
                intervals.append(interval)
            elif "Polynomial Coefficients:" in line:
                # B-Spline Format
                coeff = list(map(float, line.split(":")[1].strip().split(", ")))
                coefficients.append(coeff)
            elif line.startswith("y ="):
                # PP-Spline Format
                parts = line[4:].replace("x", "").split("+")
                slope = float(parts[0].strip())
                intercept = float(parts[1].strip())
                coefficients.append([intercept, slope])  # Order: [b, a]

    if not intervals or not coefficients:
        raise ValueError("Unsupported spline file format.")

    return intervals, coefficients

# Plot spline
def plot_spline(intervals, coefficients, output_path, title):
    plt.figure()

    for (start, end), coeff in zip(intervals, coefficients):
        x = np.linspace(start, end, 100)
        y = np.polyval(list(reversed(coeff)), x)  # Reverse coeff for numpy's format
        plt.plot(x, y)

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(title)
    plt.grid(True)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path)
    plt.close()

# Main function to process all .txt files
def process_all_splines(input_dir, output_dir):
    txt_files = glob.glob(os.path.join(input_dir, "*.txt"))

    if not txt_files:
        print(f"No .txt files found in directory: {input_dir}")
        return

    for file_path in txt_files:
        try:
            file_name = os.path.basename(file_path).replace(".txt", "")
            output_path = os.path.join(output_dir, f"{file_name}.png")

            print(f"Processing file: {file_path}")
            intervals, coefficients = read_spline(file_path)
            plot_spline(intervals, coefficients, output_path, title=f"Spline Visualization: {file_name}")
            print(f"Saved visualization: {output_path}")
        except Exception as e:
            print(f"Error processing file {file_path}: {e}")

if __name__ == "__main__":
    input_directory = "output/Check"  # Directory containing .txt files
    output_directory = "figure/Check"  # Directory to save images

    process_all_splines(input_directory, output_directory)
