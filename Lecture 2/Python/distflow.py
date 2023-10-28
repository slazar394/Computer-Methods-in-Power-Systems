import numpy as np
import pickle


def distflow(system):
    """
    Perform power flow calculations using the DistFlow method.
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
    V = np.ones(number_of_buses)
    V_old = np.zeros(number_of_buses)
    P = np.zeros(number_of_branches)
    Q = np.zeros(number_of_branches)
    P_rec = np.zeros(number_of_branches)
    Q_rec = np.zeros(number_of_branches)
    tolerance = 1e-6
    iteration = 0

    # Main loop
    while np.any(np.abs(V - V_old) >= tolerance):

        # Update the iteration variables
        V_old = np.copy(V)
        iteration += 1

        # Backward sweep: calculate active and reactive powers at the sending and receiving ends
        for i in range(number_of_branches - 1, -1, -1):
            P_rec[i] = p_load[to_bus[i]]
            Q_rec[i] = q_load[to_bus[i]]

            for j in range(number_of_branches):
                if to_bus[i] == from_bus[j]:
                    P_rec[i] += P[j]
                    Q_rec[i] += Q[j]

            # Power flows at the sending end
            P[i] = P_rec[i] + r[i] * (P_rec[i] ** 2 + Q_rec[i] ** 2) / V[to_bus[i]]
            Q[i] = Q_rec[i] + x[i] * (P_rec[i] ** 2 + Q_rec[i] ** 2) / V[to_bus[i]]

        # Forward sweep: calculate node voltages
        for i in range(number_of_branches):
            V[to_bus[i]] = V[from_bus[i]] - 2 * (P[i] * r[i] + Q[i] * x[i]) + \
                (r[i] ** 2 + x[i] ** 2) * (P[i] ** 2 + Q[i] ** 2) / V[from_bus[i]]

    # Perform voltage correction
    V = np.sqrt(V)

    return V, iteration


# Sample usage
if __name__ == '__main__':
    with open("13 bus system.pkl", "rb") as file:
        system = pickle.load(file)
    V, iteration = distflow(system)
    print("Node voltages:", V)
    print("Number of iterations:", iteration)
