from shirmohammadi import shirmohammadi
from distflow import distflow

import pickle
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


# Load the system data
with open("13 bus system.pkl", "rb") as file:
    system = pickle.load(file)

# Visualize the test system
G = nx.Graph()
edges = [(row[1], row[2]) for row in system['branches']]
G.add_edges_from(edges)
plt.figure()
pos = nx.spring_layout(G, seed=42)
nx.draw(G, pos, with_labels=True)
plt.title("Test System")
plt.show()

# Perform the power flow analysis using Shirmohammadi's method
V_sh, iteration_sh = shirmohammadi(system)

# Perform the power flow analysis using DistFlow
V_df, iteration_df = distflow(system)

# Compare the voltage magnitudes
dV = np.abs(np.abs(V_sh) - V_df)
mean_dV = np.mean(dV)
print(f"The mean voltage deviation is: {mean_dV}")

# Compare the number of iterations
print("The number of iterations to convergence:")
print(f"Shirmohammadi - {iteration_sh}")
print(f"DistFlow - {iteration_df}")
