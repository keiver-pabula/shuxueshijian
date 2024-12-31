import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

# Define parametric equations for overlay
def r1(t):
    x = np.sqrt(3) * np.cos(t)
    y = 2.0 / 3.0 * (np.sqrt(np.sqrt(3) * np.abs(np.cos(t))) + np.sqrt(3) * np.sin(t))
    return x, y

def r2(t):
    x = np.sin(t) + t * np.cos(t)
    y = np.cos(t) - t * np.sin(t)
    return x, y

def r3(t):
    x = np.sin(np.cos(t)) * np.cos(np.sin(t))
    y = np.sin(np.cos(t)) * np.sin(np.sin(t))
    z = np.cos(np.cos(t))
    return x, y, z

# Function to read data from the generated files
def read_curve_data(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if "x(t)" in line and "y(t)" in line:  # For 2D or 3D curves
                parts = line.split(",")
                x = float(parts[0].split("=")[1].strip())
                y = float(parts[1].split("=")[1].strip())
                if len(parts) == 3:  # Check if z(t) exists for 3D curves
                    z = float(parts[2].split("=")[1].strip())
                    data.append((x, y, z))
                else:
                    data.append((x, y))
    return np.array(data)

# Function to plot 2D or 3D curves and save them to the specified folder
def plot_curve(data, curve_name, resolution, parameterization, output_dir, is_3d=False):
    os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists
    if is_3d:
        t = np.linspace(0, 2 * np.pi, 1000)
        x, y, z = r3(t)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(data[:, 0], data[:, 1], data[:, 2], label=f"{curve_name}_{parameterization}_N{resolution}", color="blue")
        ax.plot(x, y, z, label=f"{curve_name}(t)", linestyle="--", color="orange")
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(f'{curve_name} ({parameterization}, 3D)')
        plt.legend()
        plt.savefig(os.path.join(output_dir, f"{curve_name}_{parameterization}_N{resolution}.png"))
        plt.close()
    else:
        if curve_name == 'r1':
            t = np.linspace(-np.pi, np.pi, 1000)
            x, y = r1(t)
        elif curve_name == 'r2':
            t = np.linspace(0, 6 * np.pi, 1000)
            x, y = r2(t)
        plt.figure()
        plt.plot(data[:, 0], data[:, 1], label=f"{curve_name}_{parameterization}_N{resolution}", color="blue")
        plt.plot(x, y, label=f"{curve_name}(t)", linestyle="--", color="orange")
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title(f'{curve_name} ({parameterization}, 2D)')
        plt.legend()
        plt.grid()
        plt.savefig(os.path.join(output_dir, f"{curve_name}_{parameterization}_N{resolution}.png"))
        plt.close()

# Main function to process all resolutions and curve files
def main():
    resolutions = [10, 40, 160]  # Different resolutions
    curves = ['r1', 'r2', 'r3']  # Curves to process
    curve_types = ['chord', 'unit']  # Both chordal and unit parameterizations

    output_dir = "figure/TaskE"  # Base output directory
    for curve in curves:
        for resolution in resolutions:
            for curve_type in curve_types:
                file_name = f"output/taskE/{curve}_{curve_type}_N{resolution}.txt"
                if os.path.exists(file_name):
                    try:
                        data = read_curve_data(file_name)
                        is_3d = (curve == 'r3')  # Only Curve 3 is 3D
                        plot_curve(data, curve, resolution, curve_type, output_dir, is_3d)
                    except Exception as e:
                        print(f"Error processing file {file_name}: {e}")
                else:
                    print(f"File not found: {file_name}")

if __name__ == "__main__":
    main()
