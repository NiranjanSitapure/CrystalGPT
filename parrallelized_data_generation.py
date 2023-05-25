# %%
import numpy as np 
import scipy.io
from IPython import get_ipython

# from matplotlib import pyplot as plt
# from matplotlib.pyplot import figure

#Import local customizable functions
import generalized_crys_kinetics as kinetics
import system_functions as sys_funcs
import pandas as pd
import random
import time

# %%
# %% [markdown]
pairs_G_and_B = 5 #Number of different (G,B) systems
N_RUN = 5000 # Number of runs for each (G,B)

# %%
# Generate 20 different [a,b,c] pairs
growth_params = []
nucleation_params = []
for i in range(pairs_G_and_B):
    growth_params.append(kinetics.random_growth_generator())
    nucleation_params.append(kinetics.random_nucleation_generator())

# Create a dictionary with the data to be exported
data_dict = {"growth_params": growth_params, "nucleation_params": nucleation_params}

# Export data to an Excel file with two sheets
with pd.ExcelWriter("different_G_and_B_params_v2.xlsx") as writer:
    for sheet_name, data in data_dict.items():
        df = pd.DataFrame(data, columns=["a", "b", "c"])
        df.to_excel(writer, sheet_name=sheet_name, index=False)

# %% [markdown]
#  Data Generation Loop (PARALLEL)
# 1. Loops over different pairs of (G,B)
# 2. Loops over 'N_RUN' number of arbitary tempearture profiles

# This uses the concurrent futures packages to use 100+ CPU (or 'N' cores) to generated
# simulation data in parallel for 20 different crystal systems, each with 5000 operating conditions 

# %%
import concurrent.futures
import simulation

if __name__ == '__main__':
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for pair in range(len(growth_params)):
            futures = [executor.submit(simulation.run_simulation, pair, run, growth_params, nucleation_params) for run in range(N_RUN)]
            concurrent.futures.wait(futures)




