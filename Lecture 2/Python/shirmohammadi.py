import numpy as np
import pickle


def shirmohammadi(system):
    """
    Perform power flow calculations using Shirmohammadi's method.
    """

    # Extract the bus data
    number_of_buses = system['buses'].shape[0]
    p_load = system['buses'][:, 1]
    q_load = system['buses'][:, 2]

    # Extract the branch data
    number_of_branches = system['branches'].shape[0]
    from_bus = system['branches'][:, 1].astype(int)
    to_bus = system['branches'][:, 2].astype(int)
    r = system['branches'][:, 3]
    x = system['branches'][:, 4]

    # Initialize the calculation variables
    S_node = - p_load - 1j * q_load
    Z = r + 1j * x
    V = np.ones(number_of_buses, dtype=complex)
    I_node = np.zeros(number_of_buses, dtype=complex)
    J_branch = np.zeros(number_of_branches, dtype=complex)
    V_old = np.zeros(number_of_buses, dtype=complex)
    tolerance = 1e-6
    iteration = 0

    # Main loop
    while np.any(np.abs(V - V_old) > tolerance):

        # Update the iteration variables
        V_old = np.copy(V)
        iteration += 1

        # Calculate the complex node current injections
        for i in range(number_of_buses):
            I_node[i] = np.conj(S_node[i] / V[i])

        # Backward sweep: calculate the complex branch currents
        for i in range(number_of_branches - 1, -1, -1):
            J_branch[i] = - I_node[to_bus[i]]

            for j in range(number_of_branches):
                if to_bus[i] == from_bus[j]:
                    J_branch[i] += J_branch[j]

        # Forward sweep: calculate the complex node voltages
        for i in range(number_of_branches):
            V[to_bus[i]] = V[from_bus[i]] - Z[i] * J_branch[i]

    return V, iteration


# Sample usage
if __name__ == '__main__':
    # Load the system data
    with open("13 bus system.pkl", "rb") as file:
        system = pickle.load(file)
    # Perform the power flow
    V, iteration = shirmohammadi(system)
    print("Node voltages:", np.abs(V))
    print("Number of iterations:", iteration)
