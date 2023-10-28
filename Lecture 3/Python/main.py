from gauss_seidel import gauss_seidel
from newton_raphson import newton_raphson

import pickle
import numpy as np


# Load the system data
with open("9 bus system.pkl", "rb") as file:
    system = pickle.load(file)

# Perform the power flow analysis using Gauss-Seidel's method
V_gs, theta_gs, iteration_gs = gauss_seidel(system)

# Perform the power flow analysis using DistFlow
V_nr, theta_nr, P, Q, iteration_nr = newton_raphson(system)

# Compare the voltage magnitudes
dV = np.abs(np.abs(V_nr) - V_gs)
mean_dV = np.mean(dV)
print(f"The mean voltage deviation is: {mean_dV}")

# Compare the number of iterations
print("The number of iterations to convergence:")
print(f"Gauss-Seidel: {iteration_gs}")
print(f"Newton-Raphson: {iteration_nr}")
