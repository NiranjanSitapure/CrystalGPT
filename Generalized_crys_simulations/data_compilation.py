# %% [markdown]
# ### Mega Data Compilation File
# 1. Here we take the data from 20/40 different crystal systems (each with a different G and B pair), and then make separate dataframes to save together in a HDF5 file (combined dataset)
# 2. Each of the dfs have data from 5000 different trials.
# 3. Plot data for 5 random cases.
# 4. Parallelize the code for efficiency

# %%
import os
import pandas as pd
import numpy as np

# %% [markdown]
# 1. Iterating over different folders to combine all the .csv in one df_combined. 
# 2. Each df from the .csv is processed to get the modified df.
# 3. df_combined for each 'job' (pair of G and B) is saved with the corresponding key in a single .H5 file.

import time
import os
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

start_time = time.time()

# Define a function to process a single job
def process_job(job):
    # Create an empty list to store the dataframes for this 'pair' folder
    pair_dfs = []

    # Loop through all the CSV files in the 'pair' folder
    directory = f'pair_{job}'
    for i, filename in enumerate(os.listdir(directory)):
        if filename.endswith('.csv'):
            # Load the CSV file into a pandas DataFrame and append it to the list for this 'pair' folder
            filepath = os.path.join(directory, filename)
            df = pd.read_csv(filepath)

            ######### Make DataFrame Changes ##########
            df = df.iloc[-13:,:]
            df = df.transpose()
            df.columns = ['conc.','temp.',
                        'moment_0','moment_1','moment_2','moment_3',
                        'suspension_density','SMD','time','T_jacket','supersat',
                        'G_rates', 'B_rates']
            df['SMD'] = df['SMD']*10**6
            # df = df.round(1)
            setpoint = df['SMD'].iloc[-1]
            error = setpoint- df['SMD']
            df['setpoint'] = setpoint*np.ones(len(error))
            df['error'] = error

            ###### Append dfs for folder 'pair' ########
            pair_dfs.append(df)

        # Print a statement every 'N' iterations
        if (i + 1) % 1000 == 0:
            print(f'Processing file {i + 1}...')

    # Concatenate all the dataframes for this 'pair' folder into a single dataframe
    df_combined = pd.concat(pair_dfs, ignore_index=True)
    print(f'########### Pair {job} Done! ##########')

    return df_combined

# Create a ThreadPoolExecutor with 10 worker threads
with ThreadPoolExecutor(max_workers=500) as executor:
    # Submit each job to the executor
    futures = [executor.submit(process_job, job) for job in range(20, 25)]

    # Create an empty dictionary to store the DataFrames
    dfs = {}

    # Wait for all the jobs to complete and get their results
    for future in as_completed(futures):
        job = future.result()

        # Store the DataFrame in the dictionary with the key 'df_i'
        dfs[f'df_{job}'] = job
       

# Save all the combined dataframes as a single HDF5 file
filename = 'finetuning_5_pairs.h5'
with pd.HDFStore(filename, mode='w') as store:
    for key, value in dfs.items():
        store[key] = value

end_time = time.time()

print(f"Elapsed time: {end_time - start_time:.2f} seconds")



# %% [markdown]
# #### Plotting 
# Plotting the combined dataset from each 'job'

# %%
# from matplotlib import pyplot as plt
# import seaborn as sns


# # Load the HDF5 file into a dictionary of pandas DataFrames
# filename = 'megadata.h5'
# dfs = {}
# with pd.HDFStore(filename, mode='r') as store:
#     for key in store.keys():
#         dfs[key[1:]] = store[key]


# # Loop over each DataFrame in the dictionary and create a plot for each one
# for key in dfs.keys():

#     plt.figure()
#     sns.lineplot(x=dfs[key]['time']/3600, y=dfs[key]['SMD'], alpha=0.05, color='blue')
#     plt.xlabel('time (h)', fontsize=18)
#     plt.ylabel('Crystal Size ($\mu m$)', fontsize=18)
#     plt.xlim([0, 25])
#     plt.ylim([0.75*dfs[key]['SMD'].min(), 1.25*dfs[key]['SMD'].max()])
#     plt.tick_params(axis='both', labelsize=18)
#     plt.grid(False)

#     plt.figure()
#     sns.lineplot(x=dfs[key]['time']/3600, y=dfs[key]['conc.'], alpha=0.05, color='pink')
#     plt.xlabel('time (h)',fontsize=18)
#     plt.ylabel('Solute Conc. (kg/kg)',fontsize=18)
#     plt.xlim([0, 25])
#     plt.ylim([0.75*dfs[key]['conc.'].min(), 1.25*dfs[key]['conc.'].max()])
#     plt.tick_params(axis='both', labelsize=18)
#     plt.grid(False)

#     plt.figure()
#     sns.lineplot(x=dfs[key]['time']/3600, y=dfs[key]['suspension_density'], alpha=0.05, color='green')
#     plt.xlabel('time (h)',fontsize=18)
#     plt.ylabel('Suspension Density (kg/kg)',fontsize=18)
#     plt.xlim([0, 25])
#     plt.ylim([0.75*dfs[key]['suspension_density'].min(), 1.25*dfs[key]['suspension_density'].max()])
#     plt.tick_params(axis='both', labelsize=18)
#     plt.grid(False)

# %% [markdown]
# ## Misc. Code

# %%
# Create an empty list to store the dataframes
# dfs = []

# # Loop through all the CSV files in the directory
# directory = 'pair_20'

# for filename in os.listdir(directory):
#     if filename.endswith('.csv'):
#         # Load the CSV file into a pandas DataFrame and append it to the list
#         filepath = os.path.join(directory, filename)
#         df = pd.read_csv(filepath)
#         ######### Make DataFrame Changes ##########
#         df = df.iloc[-13:,:]
#         df = df.transpose()
#         df.columns = ['conc.','temp.',
#                     'moment_0','moment_1','moment_2','moment_3',
#                     'suspension_density','SMD','time','T_jacket','supersat',
#                     'G_rates', 'B_rates']
#         df['SMD'] = df['SMD']*10**6
#         df = df.round(1)
#         setpoint = df['SMD'].iloc[-1]
#         error = setpoint- df['SMD']
#         df['setpoint'] = setpoint*np.ones(len(error))
#         df['error'] = error

#         ###### Append dfs for folder 'pair' ########
#         dfs.append(df)

# # Concatenate all the dataframes into a single dataframe
# combined_df = pd.concat(dfs, ignore_index=True)
# print('Dataset Size: ',combined_df.shape)
# combined_df.head()


# %% [markdown]
# How to retrieve a dataframe back from a HDF5 file. 

# %%
# # Load the HDF5 file into a pandas DataFrame
# filename = 'megadata.h5'
# df = pd.read_hdf(filename, key='df_20')

# # Print the first 5 rows of the DataFrame
# df.head()


