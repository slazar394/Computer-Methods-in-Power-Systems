import numpy as np
import pickle


def singular_transformation(system):
    """
    Create the bus-admittance matrix using the singular transformation method.
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

    # Form the branch-bus incidence matrix A and the branch admittance matrix Y
    A = np.zeros((number_of_buses, number_of_branches))
    Y = np.zeros(number_of_branches, dtype=complex)

    for i in range(number_of_branches):
        # Series elements between two independent buses
        A[from_bus[i], i] = 1
        A[to_bus[i], i] = -1
        Y[i] = 1 / (r[i] + 1j * x[i])

        # Parallel elements between the independent and the reference bus
        if b[i] != 0:
            A = np.hstack((A, np.zeros((number_of_buses, 2))))
            A[from_bus[i], -2] = 1
            A[to_bus[i], -1] = 1
            Y = np.append(Y, [0, 0])
            Y[-2] = 1j * b[i] / 2
            Y[-1] = 1j * b[i] / 2

    # Convert the branch admittance vector to an equivalent diagonal matrix
    Y = np.diag(Y)

    # Form the bus admittance matrix
    Yb = A @ Y @ A.T

    # return the bus admittance matrix
    return Yb


def nonsingular_transformation(system):
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
    # Create the bus-admittance matrix using both methods
    Y_b_1 = singular_transformation(system)
    Y_b_2 = nonsingular_transformation(system)
    # Check the difference between the two matrices
    dY = np.sum(np.abs(Y_b_1 - Y_b_2))
    print(dY)
