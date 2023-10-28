from create_Yb import create_Yb

import numpy as np
import cmath
import pickle


def gauss_seidel(system):
    """
    Perform power flow analysis using Gauss-Seidel's power flow method.
    """

    # Extract the system data
    number_of_buses = system['buses'].shape[0]
    bus_type = system['buses'][:, 1].astype(int)
    v_specified = system['buses'][:, 2]
    theta = system['buses'][:, 3]
    p_load = system['buses'][:, 4]
    q_load = system['buses'][:, 5]
    p_gen = system['buses'][:, 6]
    q_gen = system['buses'][:, 7]
    q_gen_min = system['buses'][:, 8]
    q_gen_max = system['buses'][:, 9]

    # Create the bus-admittance matrix
    Y_b = create_Yb(system)

    # Initialize the calculation variables
    p = p_gen - p_load
    q = q_gen - q_load
    v = v_specified * np.exp(1j * np.radians(theta))
    v_old = np.zeros_like(v)
    iteration = 0
    tolerance = 1e-6

    # Precompute reciprocal of diagonal elements of Y_b
    Y_b_diag_inv = 1.0 / np.diag(Y_b)

    # Main loop
    while np.any(np.abs(v - v_old) > tolerance):
        # Update the iteration variables
        v_old = np.copy(v)
        iteration += 1

        # Apply the Gauss-Seidel power flow equations to each node except the slack node
        for i in range(1, number_of_buses):
            # Calculate the sum of node current contributions
            S = np.dot(Y_b[i, :], v)
            S -= Y_b[i, i] * v[i]

            if bus_type[i] == 2:
                # PV node
                q[i] = -np.imag(np.conj(v[i]) * (S + Y_b[i, i] * v[i]))
                q_gen[i] = q[i] + q_load[i]

                # Check the generator reactive power limits
                if q_gen[i] < q_gen_min[i] or q_gen[i] > q_gen_max[i]:
                    q_gen[i] = max(min(q_gen[i], q_gen_max[i]), q_gen_min[i])
                    q[i] = q_gen[i] - q_load[i]
                    v[i] = Y_b_diag_inv[i] * ((p[i] - 1j * q[i]) / np.conj(v[i]) - S)
                else:
                    v[i] = Y_b_diag_inv[i] * ((p[i] - 1j * q[i]) / np.conj(v[i]) - S)
                    v[i] = v_specified[i] * cmath.exp(1j * cmath.phase(v[i]))

            else:
                # PQ node
                v[i] = Y_b_diag_inv[i] * ((p[i] - 1j * q[i]) / np.conj(v[i]) - S)

    # Evaluate the voltage magnitudes and phase angles
    theta = np.angle(v, deg=True)
    v = np.abs(v)

    return v, theta, iteration


# Sample usage
if __name__ == '__main__':
    # Load the system data
    with open("9 bus system.pkl", "rb") as file:
        system = pickle.load(file)
    # Perform the power flow
    v, theta, iteration = gauss_seidel(system)
    print("Voltage magnitudes:", v)
    print("Voltage phase angles:", theta)
