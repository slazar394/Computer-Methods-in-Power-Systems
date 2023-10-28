import numpy as np
import pickle


def create_Yb(system):
    """
    Create the bus-admittance matrix using the non-singular transformation method.
    """

    # Determine the number of independent buses and the number of branches in the system
    number_of_buses = system['buses'].shape[0]
    number_of_branches = system['branches'].shape[0]

    # Extract the branch data
    from_bus = system['branches'][:, 0].astype(int)
    to_bus = system['branches'][:, 1].astype(int)
    r = system['branches'][:, 2]
    x = system['branches'][:, 3]
    b = system['branches'][:, 4]

    # Initialize the bus-admittance matrix
    Yb = np.zeros((number_of_buses, number_of_buses), dtype=complex)

    # Form the bus-admittance matrix
    for i in range(number_of_branches):
        # Off-diagonal elements
        Yb[from_bus[i], to_bus[i]] -= 1 / (r[i] + 1j * x[i])
        Yb[to_bus[i], from_bus[i]] -= 1 / (r[i] + 1j * x[i])

        # Diagonal elements
        Yb[from_bus[i], from_bus[i]] += 1 / (r[i] + 1j * x[i]) + 1j * b[i] / 2
        Yb[to_bus[i], to_bus[i]] += 1 / (r[i] + 1j * x[i]) + 1j * b[i] / 2

    return Yb


if __name__ == "__main__":
    # Load the system data
    with open("9 bus system.pkl", "rb") as file:
        system = pickle.load(file)
    # Create the bus-admittance matrix
    Y_b = create_Yb(system)
    # Print the matrix
    print(Y_b)
