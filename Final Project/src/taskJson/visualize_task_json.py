import json
import os
import matplotlib.pyplot as plt
import numpy as np


def load_json_config(file_path):
    """
    Loads a JSON configuration file.
    """
    with open(file_path, 'r') as file:
        return json.load(file)


def generate_and_visualize(config):
    """
    Generates and visualizes curves based on the JSON configuration.
    """
    # Create output directories for text files and figures
    txt_output_dir = os.path.join("Output", "taskJson")
    figure_output_dir = os.path.join("figure", "taskJson")
    os.makedirs(txt_output_dir, exist_ok=True)
    os.makedirs(figure_output_dir, exist_ok=True)

    print(f"Generating outputs in {txt_output_dir} and visualizations in {figure_output_dir}...\n")

    for curve in config["curves"]:
        name = curve["name"]
        parameterization = curve["parameterization"]
        resolution = curve["resolution"]

        # Generate data points based on parameterization
        t = np.linspace(0, 1, resolution)
        if name == "r1":
            x = np.sqrt(3) * np.cos(t * 2 * np.pi)
            y = (2.0 / 3.0) * (np.sqrt(np.sqrt(3) * np.abs(np.cos(t * 2 * np.pi))) + np.sqrt(3) * np.sin(t * 2 * np.pi))
            data = np.column_stack((x, y))
        elif name == "r2":
            x = np.sin(t * 2 * np.pi) + t * np.cos(t * 2 * np.pi)
            y = np.cos(t * 2 * np.pi) - t * np.sin(t * 2 * np.pi)
            data = np.column_stack((x, y))
        elif name == "r3":
            x = np.sin(np.cos(t * 2 * np.pi)) * np.cos(np.sin(t * 2 * np.pi))
            y = np.sin(np.cos(t * 2 * np.pi)) * np.sin(np.sin(t * 2 * np.pi))
            z = np.cos(np.cos(t * 2 * np.pi))
            data = np.column_stack((x, y, z))
        else:
            print(f"Unknown curve name: {name}. Skipping...")
            continue

        # Save data to a .txt file
        txt_file_path = os.path.join(txt_output_dir, f"{name}_{parameterization}_N{resolution}.txt")
        with open(txt_file_path, "w") as txt_file:
            if name == "r3":
                txt_file.write("x(t), y(t), z(t)\n")
                for row in data:
                    txt_file.write(f"{row[0]:.6f}, {row[1]:.6f}, {row[2]:.6f}\n")
            else:
                txt_file.write("x(t), y(t)\n")
                for row in data:
                    txt_file.write(f"{row[0]:.6f}, {row[1]:.6f}\n")
        print(f"Saved {txt_file_path}")

        # Plot the data
        plt.figure()
        if name == "r3":
            # 3D plot for r3
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot(data[:, 0], data[:, 1], data[:, 2], label=f"{name} ({parameterization}, N={resolution})", color="blue")
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.set_zlabel("Z")
            ax.set_title(f"3D Visualization for {name}")
        else:
            # 2D plot for r1 and r2
            plt.plot(data[:, 0], data[:, 1], label=f"{name} ({parameterization}, N={resolution})", color="blue")
            plt.xlabel("X")
            plt.ylabel("Y")
            plt.title(f"2D Visualization for {name}")
            plt.grid()

        plt.legend()
        image_file_path = os.path.join(figure_output_dir, f"{name}_{parameterization}_N{resolution}.png")
        plt.savefig(image_file_path)
        print(f"Saved {image_file_path}")
        plt.close()

    print("All outputs and visualizations have been generated successfully!")


def main():
    # Path to the JSON configuration file
    json_file_path = "config/config.json"

    # Check if the file exists
    if not os.path.exists(json_file_path):
        print(f"Error: Configuration file not found at {json_file_path}.")
        return

    # Load the JSON configuration
    config = load_json_config(json_file_path)

    # Generate visualizations and output .txt files
    generate_and_visualize(config)


if __name__ == "__main__":
    main()
