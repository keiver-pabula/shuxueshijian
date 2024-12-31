import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime

def read_data(file_path):
    """
    Reads intervals and polynomial coefficients from a file.

    Args:
        file_path (str): Path to the input file.

    Returns:
        tuple: A tuple containing a list of intervals and a list of polynomials.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    intervals = []
    polynomials = []
    
    i = 0
    while i < len(lines):
        try:
            # Try to parse the interval
            interval = list(map(float, lines[i].strip().split(',')))
            intervals.append(interval)
            
            # Parse the polynomial coefficients
            polynomial = list(map(float, lines[i + 1].strip().split()))
            polynomials.append(polynomial)
            
            i += 2  # Move to the next pair of lines
        except ValueError:
            # Skip lines that cannot be parsed
            i += 1
            continue
    
    return intervals, polynomials


def original_function(x):
    """
    The original function to visualize alongside the polynomials.
    """
    return 1 / (1 + 25 * x**2)

def plot_polynomials(intervals, polynomials, output_path, plot_type):
    """
    Plots the original function and the piecewise polynomials.

    Args:
        intervals (list): List of intervals.
        polynomials (list): List of polynomial coefficients.
        output_path (str): Path to save the plot.
        plot_type (str): Type of the plot (e.g., 'cubic', 'quadratic').
    """
    plt.figure()
    
    # Plot the original function
    x = np.linspace(-1, 1, 1000)
    y = original_function(x)
    plt.plot(x, y, 'k--', label='Original function')
    
    # Plot the piecewise polynomials
    for interval, coeffs in zip(intervals, polynomials):
        x = np.linspace(interval[0], interval[1], 400)
        y = np.polyval(coeffs[::-1], x)
        plt.plot(x, y, label=f'Interval {interval}')
    
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'Piecewise Polynomial Plot ({plot_type.capitalize()})')
    plt.legend()
    plt.grid(True)
    plt.savefig(output_path)
    plt.close()

def log_message(message):
    """
    Logs a message to a log file with a timestamp.
    """
    log_file = 'logs/log.txt'
    timestamp = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    with open(log_file, 'a') as log:
        log.write(f'{timestamp} : {message}\n')

def main():
    input_dir = './output/TaskB'
    output_dir = './figure/TaskB'
    plot_types = ['cubic', 'quadratic']
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for plot_type in plot_types:
        input_file = os.path.join(input_dir, f'{plot_type}.txt')
        output_file = os.path.join(output_dir, f'{plot_type}.png')
        
        if not os.path.exists(input_file):
            log_message(f'Error: {input_file} does not exist.')
            continue
        
        log_message(f'Reading data from {input_file}')
        intervals, polynomials = read_data(input_file)
        
        log_message(f'Plotting data and saving to {output_file}')
        plot_polynomials(intervals, polynomials, output_file, plot_type)
        
        log_message(f'Successfully saved {output_file}')

if __name__ == '__main__':
    main()
