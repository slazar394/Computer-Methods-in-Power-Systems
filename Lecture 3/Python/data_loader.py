import pandas as pd
import pickle

# Load the bus and branch data from the .xlsx file and convert the dataframes to arrays
buses = pd.read_excel("9 bus system.xlsx", sheet_name="Bus data").values
branches = pd.read_excel("9 bus system.xlsx", sheet_name="Branch data").values

# Renumber the bus ids to start from 0
buses[:, 0] -= 1
branches[:, [0, 1]] -= 1

# Store the bus and branch data in a dictionary
system = {
    'buses': buses,
    'branches': branches
}

print(system)

# Save the data using pickle
with open("9 bus system.pkl", "wb") as file:
    pickle.dump(system, file)
