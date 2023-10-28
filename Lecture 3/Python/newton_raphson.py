from create_Yb import create_Yb

import numpy as np
import pickle


def newton_raphson(system):
    """
    Perform power flow analysis using Newton-Raphson's power flow method.
    """

    # Extract the system data
    number_of_buses = system['buses'].shape[0]
    bus_type = system['buses'][:, 1].astype(int)
    pq_buses = np.where(bus_type == 3)[0]
    v_specified = system['buses'][:, 2]
    theta = system['buses'][:, 3]
    p_load = system['buses'][:, 4]
    q_load = system['buses'][:, 5]
    p_gen = system['buses'][:, 6]
    q_gen = system['buses'][:, 7]

    # Calculate the specified active and reactive power injections
    p_specified = p_gen - p_load
    q_specified = q_gen - q_load

    # Create the bus-admittance matrix and evaluate its real and imaginary parts
    Y_b = create_Yb(system)
    G = np.real(Y_b)
    B = np.imag(Y_b)

    # Initialize the calculation variables
    V = np.copy(v_specified)
    V_old = np.zeros_like(V)
    iteration = 0
    tolerance = 1e-6

    # Main loop
    while np.any(np.abs(V * np.exp(1j * theta) - V_old) > tolerance):
        # Update the iteration variables
        V_old = V * np.exp(1j * theta)
        iteration += 1

        # Determine active and reactive power injections
        P = np.zeros(number_of_buses)
        Q = np.zeros(number_of_buses)
        for i in range(number_of_buses):
            for j in range(number_of_buses):
                P[i] += V[i] * V[j] * (G[i, j] * np.cos(theta[i] - theta[j]) + B[i, j] * np.sin(theta[i] - theta[j]))
                Q[i] += V[i] * V[j] * (G[i, j] * np.sin(theta[i] - theta[j]) - B[i, j] * np.cos(theta[i] - theta[j]))

        # Determine the active and reactive power injection deviations
        dP = p_specified - P
        dP = dP[1:]
        dQ = q_specified - Q
        dQ = dQ[pq_buses]

        # Create J1 = dP/dtheta
        J1 = np.zeros((number_of_buses, number_of_buses))
        for i in range(number_of_buses):
            for j in range(number_of_buses):
                if i == j:
                    for k in range(number_of_buses):
                        J1[i, i] += V[i] * V[k] * (
                                    -G[i, k] * np.sin(theta[i] - theta[k]) + B[i, k] * np.cos(theta[i] - theta[k]))
                    J1[i, i] -= V[i] ** 2 * B[i, i]
                else:
                    J1[i, j] = V[i] * V[j] * (
                                G[i, j] * np.sin(theta[i] - theta[j]) - B[i, j] * np.cos(theta[i] - theta[j]))
        J1 = J1[1:, 1:]

        # Create J2 = dP/dV
        J2 = np.zeros((number_of_buses, number_of_buses))
        for i in range(number_of_buses):
            for j in range(number_of_buses):
                if i == j:
                    for k in range(number_of_buses):
                        J2[i, i] += V[k] * (G[i, k] * np.cos(theta[i] - theta[k]) + B[i, k] * np.sin(theta[i] - theta[k]))
                    J2[i, i] += G[i, i] * V[i]
                else:
                    J2[i, j] = V[i] * (G[i, j] * np.cos(theta[i] - theta[j]) + B[i, j] * np.sin(theta[i] - theta[j]))
        J2 = J2[1:, pq_buses]

        # Create J3 = dQ/dtheta
        J3 = np.zeros((number_of_buses, number_of_buses))
        for i in range(number_of_buses):
            for j in range(number_of_buses):
                if i == j:
                    for k in range(number_of_buses):
                        J3[i, i] += V[i] * V[k] * (G[i, k] * np.cos(theta[i] - theta[k]) + B[i, k] * np.sin(theta[i] - theta[k]))
                    J3[i, i] -= G[i, i] * V[i] ** 2
                else:
                    J3[i, j] = -V[i] * V[j] * (G[i, j] * np.cos(theta[i] - theta[j]) + B[i, j] * np.sin(theta[i] - theta[j]))
        J3 = J3[pq_buses, 1:]

        # Create J4 = dQ/dtheta
        J4 = np.zeros((number_of_buses, number_of_buses))
        for i in range(number_of_buses):
            for j in range(number_of_buses):
                if i == j:
                    for k in range(number_of_buses):
                        J4[i, i] += V[k] * (G[i, k] * np.sin(theta[i] - theta[k]) - B[i, k] * np.cos(theta[i] - theta[k]))
                    J4[i, i] -= B[i, i] * V[i]
                else:
                    J4[i, j] = V[i] * (G[i, j] * np.sin(theta[i] - theta[j]) - B[i, j] * np.cos(theta[i] - theta[j]))
        J4 = J4[pq_buses, :][:, pq_buses]

        # Create the complete Jacobian matrix
        J = np.block([[J1, J2], [J3, J4]])

        # Solve the linear system of equations
        dX = np.linalg.solve(J, np.concatenate([dP, dQ]))

        # Update the voltage magnitudes and voltage phase angles
        dTh = dX[:number_of_buses - 1]
        dV = dX[number_of_buses - 1:]
        theta[1:] += dTh
        V[pq_buses] += dV

    # Convert theta to degrees
    theta = np.rad2deg(theta)

    return V, theta, P, Q, iteration


# Sample usage
if __name__ == '__main__':
    # Load the system data
    with open("9 bus system.pkl", "rb") as file:
        system = pickle.load(file)
    # Perform the power flow
    v, theta, P, Q, iteration = newton_raphson(system)
    print("Voltage magnitudes:", v)
    print("Voltage phase angles:", theta)

