import pandas as pd
import matplotlib.pyplot as plt
import os

# Ensure output and figure directories exist
output_dir = './output/taskA'
figure_dir = './figure/taskA'
os.makedirs(output_dir, exist_ok=True)
os.makedirs(figure_dir, exist_ok=True)

subdivisions = [6, 11, 21, 41, 81]
colors = ['cyan', 'green', 'red', 'blue', 'magenta']

plt.figure(figsize=(8, 6))

# Plot exact function
exact_data = pd.read_csv("comparison_n81.csv")
t = exact_data['t']
exact = exact_data['exact']
plt.plot(t, exact, label='Exact Function', linestyle='--', color='black')

# Save exact data to a .txt file
with open(f"{output_dir}/exact.txt", 'w') as f:
    f.write("t, exact\n")
    for t_val, exact_val in zip(t, exact):
        f.write(f"{t_val:.6f}, {exact_val:.6f}\n")

# Plot splines for different subdivisions
for n, color in zip(subdivisions, colors):
    data = pd.read_csv(f"comparison_n{n}.csv")
    t_values = data['t']
    spline_values = data['spline']
    plt.plot(t_values, spline_values, label=f'N={n}', color=color)
    
    # Save each subdivision data to a .txt file
    with open(f"{output_dir}/spline_n{n}.txt", 'w') as f:
        f.write("t, spline\n")
        for t_val, spline_val in zip(t_values, spline_values):
            f.write(f"{t_val:.6f}, {spline_val:.6f}\n")

# Add labels, title, and legend
plt.title("Problem A")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid()

# Save the plot to a .png file
plt.savefig(f"{figure_dir}/problem_a.png")
plt.close()

print("Outputs generated:")
print(f"Data files saved in {output_dir}/")
print(f"Figure saved as {figure_dir}/problem_a.png")
