import numpy as np
import pickle


def DC_power_flow(system):
    """
    Perform the power flow using the linear DC method.
    """

    # Determine the number of buses and the number of branches in the system
    number_of_buses = system['buses'].shape[0]
    number_of_branches = system['branches'].shape[0]

    # Extract the system data
    p_gen = system['buses'][:, 6]
    p_load = system['buses'][:, 4]
    from_bus = system['branches'][:, 0].astype(int)
    to_bus = system['branches'][:, 1].astype(int)
    x = system['branches'][:, 3]

    # Create the bus-susceptance matrix
    B = np.zeros((number_of_buses, number_of_buses))

    for i in range(number_of_branches):
        # Off-diagonal elements
        B[from_bus[i], to_bus[i]] += 1 / x[i]
        B[to_bus[i], from_bus[i]] += 1 / x[i]

        # Diagonal elements
        B[from_bus[i], from_bus[i]] -= 1 / x[i]
        B[to_bus[i], to_bus[i]] -= 1 / x[i]

    # Remove the row and column associated with the slack bus to create a reduced bus-susceptance matrix
    B_r = np.delete(B, 0, axis=0)
    B_r = np.delete(B_r, 0, axis=1)

    # Create the bus injection vector
    p = p_gen - p_load

    # Remove the row associated with the slack bus to create a reduced bus-injection vector
    p_r = p[1:]

    # Determine the voltage phase angles
    theta = np.zeros(number_of_buses)
    theta[1:] = -np.linalg.inv(B_r) @ p_r

    # Calculate the branch active power flows
    p_branch = np.zeros(number_of_branches)

    for i in range(number_of_branches):
        p_branch[i] = (theta[from_bus[i]] - theta[to_bus[i]]) / x[i]

    return theta, p_branch


if __name__ == '__main__':
    # Load the system data
    with open("9 bus system.pkl", "rb") as file:
        system = pickle.load(file)
    # Perform the DC power flow
    theta, p_branch = DC_power_flow(system)
    # Print the results
    print(theta)
    print(p_branch)
