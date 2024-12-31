import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime

def read_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    intervals = []
    polynomials = []
    current_intervals = []
    current_polynomials = []
    
    for line in lines:
        line = line.strip()
        # Skip lines that don't contain intervals or polynomial coefficients
        if line.startswith("Fitting") or line.startswith("Knots") or line == "":
            continue
        if line.startswith("Interval:"):
            # Parse the interval
            interval = list(map(float, line.split(":")[1].strip().strip("[]").split(",")))
            current_intervals.append(interval)
        elif line.startswith("Polynomial Coefficients:"):
            # Parse the coefficients
            coefficients = list(map(float, line.split(":")[1].strip().split(",")))
            current_polynomials.append(coefficients)
        elif line == "====":
            if current_intervals and current_polynomials:
                intervals.append(current_intervals)
                polynomials.append(current_polynomials)
                current_intervals = []
                current_polynomials = []
    
    # Add the last set of intervals and polynomials if they exist
    if current_intervals and current_polynomials:
        intervals.append(current_intervals)
        polynomials.append(current_polynomials)
    
    return intervals, polynomials

def plot_polynomials(intervals, polynomials, output_path):    

    # Plot the piecewise polynomials
    for idx, (interval_set, polynomial_set) in enumerate(zip(intervals, polynomials)):
        for interval, coeffs in zip(interval_set, polynomial_set):
            x = np.linspace(interval[0], interval[1], 400)
            y = np.polyval(coeffs[::-1], x)
            plt.plot(x, y)  # Removed the label

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Piecewise Polynomial Plot')
    plt.grid(True)

    
    plt.savefig(output_path)
    plt.close()

def log_message(message):
    log_file = 'logs/log.txt'
    timestamp = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    if not os.path.exists(os.path.dirname(log_file)):
        os.makedirs(os.path.dirname(log_file))
    with open(log_file, 'a') as log:
        log.write(f'{timestamp} : {message}\n')

def main():
    input_file = './output/taskTest/func_results.txt'
    output_dir = './figure/taskTest'
    output_file = os.path.join(output_dir, 'func_plot.png')
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    intervals, polynomials = read_data(input_file)
    plot_polynomials(intervals, polynomials, output_file)
    log_message(f'Successfully saved {output_file}')

if __name__ == '__main__':
    main()
