import pandas as pd
import matplotlib.pyplot as plt
import os

# Create directories for visualization
os.makedirs("figure/taskD", exist_ok=True)

# Load data
data = pd.read_csv("output/taskD/task_d_spline_comparison.txt", delim_whitespace=True)

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(data["t"], data["cubic"], label="Cubic Spline", linestyle="--", color="blue")
plt.plot(data["t"], data["quadratic"], label="Quadratic B-Spline", linestyle="-", color="green")

# Add markers for nodes (subset for better visualization)
plt.scatter(data["t"][::10], data["cubic"][::10], color="red", label="Nodes (Cubic)")
plt.scatter(data["t"][::10], data["quadratic"][::10], color="orange", label="Nodes (Quadratic)")

# Add labels and legend
plt.title("Comparison of Cubic Spline and Quadratic B-Spline (Task D)")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.legend()
plt.grid()

# Save the plot
output_file = "figure/taskD/task_d_spline_visualization.png"
plt.savefig(output_file)
print(f"Visualization saved to '{output_file}'")
plt.show()
