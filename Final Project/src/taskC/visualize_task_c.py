import pandas as pd
import matplotlib.pyplot as plt
import os

def load_data(file_path):
    """
    Load data from a .txt file with space-separated values.
    """
    return pd.read_csv(file_path, delim_whitespace=True)

def plot_data(data, output_folder, output_file):
    """
    Plot ppForm and B-Spline interpolation data.
    """
    plt.figure(figsize=(8, 6))
    plt.plot(data['t'], data['ppForm'], label='pp-Form (Cubic Spline)', linestyle='--', color='blue')
    plt.plot(data['t'], data['bSpline'], label='B-Spline', linestyle='-', color='red')

    # Add labels, legend, and title
    plt.title("Comparison of pp-Form and B-Spline Interpolation")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.grid()

    # Save the plot
    os.makedirs(output_folder, exist_ok=True)
    plt.savefig(os.path.join(output_folder, output_file))
    plt.close()

def main():
    input_file = "output/taskC/comparison_pp_b_spline.txt"
    output_folder = "figure/taskC"
    output_file = "comparison_pp_b_spline.png"

    if not os.path.exists(input_file):
        print(f"Error: {input_file} does not exist.")
        return

    data = load_data(input_file)
    plot_data(data, output_folder, output_file)
    print(f"Plot saved to {os.path.join(output_folder, output_file)}")

if __name__ == "__main__":
    main()
