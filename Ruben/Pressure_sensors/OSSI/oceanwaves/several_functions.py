#functions to use in the laboratorioes of ocean waves
import numpy as np
import matplotlib.pyplot as plt
import pickle

#solve the dispertion relation using Newton-Rhapson
def waveNumber_dispertionRelation_NewtonRhapson(T,d,tolerance):
    """
    Calculate wave number using the Newton-Raphson method for dispersion relation.
    
    Args:
        T (float): Wave period.
        d (float): Water depth.
        tolerance (float): Tolerance for convergence.
        
    Returns:
        float: Calculated wave number.
    """
    g = 9.81  # Acceleration due to gravity
    omega = 2 * np.pi / T  # Angular frequency
    k0 = omega * omega / g  # Initial guess for wave number
    k_guess = k0  # Initialize the initial guess for the wave number

    # Iterating using the Newton-Raphson method to find the wave number
    for i in range(1000):
        fx = k_guess * np.tanh(k_guess * d) - k0  # Function value
        fxp = d * (k_guess - k_guess * (np.tanh(k_guess * d)) ** 2) + np.tanh(k_guess * d)  # Derivative of the function
        k = k_guess - fx / fxp  # Updated guess for the wave number

        # Checking if the relative change in the wave number is within tolerance
        if np.abs((k_guess - k) / k) < tolerance:
            break  # If tolerance is met, exit the loop
        k_guess = k  # Update the guess for the next iteration

    return k  # Return the final calculated wave number

#plot in one line x vs y 
def easyplot(x, y, x_limits=None, y_limits=None, xlabel=None, legend=None, title=None ,variable=None):
    """
    Plots the relationship between x and y with optional x-axis limits and adjusts the axis based on one or two limits

    Args:
        x (list or numpy.ndarray): The x-values.
        y (list or numpy.ndarray): The y-values.
        x_limits/y_limits (tuple or None): Tuple specifying the lower and upper limits of the axis (default: None).
        xlabel/legend/title = Strings to identify the plot
        
    Returns:
        None
    """
    line1, = plt.plot(x, y)
    
    if x_limits:
        plt.xlim(x_limits[0], x_limits[1])
        line1_xdata = line1.get_xdata()
        line1_ydata = line1.get_ydata()

        mask = (line1_xdata >= plt.xlim()[0]) & (line1_xdata <= plt.xlim()[1])
        line1_ymin = line1_ydata[mask].min()
        line1_ymax = line1_ydata[mask].max()

        # Automatically adjust the y-axis limits 
        plt.ylim(line1_ymin, line1_ymax)
        
    if y_limits:
        plt.ylim(y_limits[0], y_limits[1])
        line1_xdata = line1.get_xdata()
        line1_ydata = line1.get_ydata()

        mask = (line1_ydata >= plt.ylim()[0]) & (line1_ydata <= plt.ylim()[1])
        line1_xmin = line1_xdata[mask].min()
        line1_xmax = line1_xdata[mask].max()

        # Automatically adjust the x-axis limits 
        plt.xlim(line1_xmin, line1_xmax)
        
    if xlabel:
        plt.xlabel(xlabel)
    if legend:
        plt.legend([legend])
    if title:
        plt.title(title)
    if variable:
        line1.set_label(variable)
    plt.grid(True)
    plt.legend()
    plt.legend(frameon=False)
    return

# function
def compute_variables(d, T, Hoff):
    alpha = 9.81
    beta = 2 * np.pi / T
    gamma = beta * beta / alpha
    delta = gamma * d
    epsilon = delta * (np.tanh(delta)) ** -0.5
    zeta = (delta + epsilon ** 2 * np.cosh(epsilon) ** -2) / (np.tanh(epsilon) + epsilon * np.cosh(epsilon) ** -2) / d

    eta = 2 * np.pi / zeta
    theta = eta / T
    iota = 0.5 + (zeta * d / np.sinh(2 * zeta * d))
    kappa = iota * theta

    lambda_var = np.sqrt(kappa[0] / kappa)
    mu = Hoff * lambda_var

    nu = 0.8
    xi = nu * d
    mu[mu > xi] = xi[mu > xi]

    return eta, iota, theta, kappa, lambda_var, mu
    
    
# Save variables to a file
def save_variables(filename, **kwargs):
    with open(filename, 'wb') as file:
        pickle.dump(kwargs, file)

# Read variables from a file and assign them different names
def read_variables(filename, **new_names):
    read_vars = {}
    with open(filename, 'rb') as file:
        data = pickle.load(file)

    for old_name, new_name in new_names.items():
        if old_name in data:
            read_vars[new_name] = data[old_name]
        else:
            print(f"Variable '{old_name}' does not exist in the file.")

    return read_vars