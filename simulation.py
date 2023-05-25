
from functools import partial
import system_functions as sys_funcs
import generalized_crys_kinetics as kinetics


import random
import numpy as np
import os

from scipy.integrate import odeint, solve_ivp
from scipy.optimize import minimize_scalar


import signal
# Allows to skip a certain simulation if it takes too much time or is stuck for certain 'N' seconds. 
def signal_handler(signum, frame):
    raise TimeoutError("solve_ivp timed out")

# Model Initialization (same as the Crystallization Files)
R  = 8.314 # Universal gas constant
volume_shape_factor = 1; #Shape factor for cubic crystals
size_factor = 10**-6 # Converting the um into cm
rho_crystal = 1540 #kg/m^3 #Crystal density for Dextrose Monohydrate

#Decide the size range of the CSD based on prior/literature
minsize = 1 # um0
maxsize = 500 # um
intervals = 250 # This is the # intervals to bin the crystals by their size


# Initial Concentration
init_solution = random.uniform(0.75, 0.95)  # Initial weight of the solution is 1-0.125 (seed crystals)
init_concentration = random.uniform(0.6, 0.8)  # kg_dextrose/kg_slurry
init_DX = init_concentration / init_solution

# Initial CSD
init_crystals = 5e8  # (#/m^3)
size_range = [50, 200]  # um
sigma_range = [5, 25]  # um

initial_mean_size = np.random.uniform(size_range[0], size_range[1])
initial_sigma_size = np.random.uniform(sigma_range[0], sigma_range[1])
initial_y = np.linspace(minsize*size_factor, maxsize*size_factor, num=intervals) 
initial_boundaries = np.linspace(minsize, maxsize, intervals+1)*size_factor
initial_F = sys_funcs.normal_dist(initial_y/size_factor/maxsize, initial_mean_size/maxsize, initial_sigma_size/maxsize)



 ################ SOLVING TIME AND DESCRETIZATION ########################## 
time_to_solve = 3600*24 # in seconds # This is the crystallization time
descretize_time = 600 #seconds = 10 min.
sampling_time = 60*60 # 30 mins (when Tj changes)

time_intervals = round(time_to_solve/descretize_time) + 1  # This is the number of steps in the ODE solution 
time_series = np.linspace(0, time_to_solve, int(time_intervals)) 


# Need to make function to convert a sampling_temp value into an array of sampling minutes
def make_temp_series(Temp_value):
    temp_series = np.ones(int(sampling_time/descretize_time))*Temp_value
    temp_series= np.round(temp_series,2)
    return temp_series

# Define a function to generate a random temperature profile given a time array
def temp_generation(time_array):
    # Get the time samples for temperature measurement
    time_samples = time_array[time_array % sampling_time == 0]
    # Initialize the temperature list with a random initial temperature
    temperature = []
    max_deltaT = random.uniform(2, 5)
    deltaT = random.uniform(-max_deltaT, max_deltaT*3/5)
    T_initial = random.uniform(35, 40) + deltaT
    temp_now = T_initial
    temperature.append(temp_now)
    # Generate the temperature profile for the first time interval
    temperature.extend(make_temp_series(temp_now))
    # Generate the temperature profile for subsequent time intervals
    for t in time_samples[1:-1]:
        deltaT = random.uniform(-max_deltaT, max_deltaT*3/5)
        temp_now = temp_now + deltaT
        # Keep the temperature between 5 and 40 degrees Celsius
        if temp_now < 5:
            temp_now = 5
        elif temp_now > 40:
            temp_now = 40
        temperature.extend(make_temp_series(temp_now))
    return np.array(np.round(temperature, 2))

# Generate the initial cooling profile using the temp_generation function
initial_cooling_profile = np.asarray(temp_generation(time_series))
initial_cooling_profile = np.round(initial_cooling_profile, 2)


# Constants
rho_crystal = 1540  # Density for dextrose crystals
volume_shape_factor = 1  # For cubic crystals
heat_transfer_coeff = 200  # W/m^2/K
Cp_slurry = 4182  # J/kg/k
jacket_area = 0.01  # m^2
mass_slurry = 1  # kg (this is our basis for the calculation)

def differential_function(t, X, sample_y, cooling_profile, param_growth, param_nucleation):

    """
    This function takes in time (t), concentrations/size distribution (X), initial boundaries for crystal sizes,
    crystal size distribution (sample_y), and the cooling profile for the reactor (cooling_profile), and returns
    the derivatives of each variable.
    """
    y = sample_y
    Dy = np.diff(initial_boundaries)
    ya = np.hstack((y[1:None], y[-1]))
    yb = np.hstack((y[0], y[0:-2], y[-2]))

    # Unpack variables
    F = X[0:-2]  # CSD at time 't'
    c = X[-2]  # Solute concentration
    T = X[-1] + 273.15  # Slurry temperature

    # Determine the temperature of the jacket based on the cooling profile
    T_jacket = cooling_profile[int(np.floor(t / descretize_time))] + 273

    # Calculate moments and solubility
    M_T = sys_funcs.moments(F, 3, initial_boundaries, sample_y) * volume_shape_factor * rho_crystal
    c_dextrose = c / (1 - M_T)  # converting kg_dext/kg_slurry to kg_dext/kg_solution
    S = c_dextrose / kinetics.solubility(T)  # Solubility and saturation is measured per kg_solution

    # Growth terms
    G = kinetics.growth_rate(S, T, y,param_growth)
    Ga = kinetics.growth_rate(S, T, ya,param_growth)
    Fa = np.hstack((F[1:None], F[-1]))
    Gb = kinetics.growth_rate(S, T, yb,param_growth)
    Fb = np.hstack((F[0], F[0:-2], F[-2]))

    # Nucleation term
    J = kinetics.nucleation_rate(S, M_T, param_nucleation)
    Fb[0] = J / G[0]  # Nucleation boundary condition

    # Growth derivative
    dF = - (Ga * Fa - Gb * Fb) / (ya - yb)
    dF[np.isnan(dF)] = 0

    # Concentration derivative
    dc = -3 * volume_shape_factor * rho_crystal * np.sum(G * F * y ** 2 * Dy) - J * y[0] ** 3 * volume_shape_factor * rho_crystal

    # Temperature derivative
    dT = heat_transfer_coeff * jacket_area * (T_jacket - T) / (mass_slurry * Cp_slurry)

    # Combine the differentials in one array
    dXdt = np.hstack((dF, dc, dT))

    return dXdt


def run_simulation(pair, run, growth_params, nucleation_params):
        
        ################ Model Initialization [Basis of 1 kg_solution] ##########################

        R  = 8.314 # Universal gas constant
        volume_shape_factor = 1; #Shape factor for cubic crystals
        size_factor = 10**-6 # Converting the um into cm
        rho_crystal = 1540 #kg/m^3 #Crystal density for Dextrose Monohydrate

        #Decide the size range of the CSD based on prior/literature
        minsize = 1 # um0
        maxsize = 500 # um
        intervals = 250 # This is the # intervals to bin the crystals by their size
       
        # Initial Concentration
        init_solution = random.uniform(0.75, 0.95)  # Initial weight of the solution is 1-0.125 (seed crystals)
        init_concentration = random.uniform(0.6, 0.85)  # kg_dextrose/kg_slurry
        init_DX = init_concentration / init_solution

        # Initial CSD
        init_crystals = 5e8  # (#/m^3)
        size_range = [50, 150]  # um
        sigma_range = [5, 20]  # um

        initial_mean_size = np.random.uniform(size_range[0], size_range[1])
        initial_sigma_size = np.random.uniform(sigma_range[0], sigma_range[1])
        initial_y = np.linspace(minsize*size_factor, maxsize*size_factor, num=intervals) 
        initial_boundaries = np.linspace(minsize, maxsize, intervals+1)*size_factor
        initial_F = sys_funcs.normal_dist(initial_y/size_factor/maxsize, initial_mean_size/maxsize, initial_sigma_size/maxsize)


        # Optimize number of initial crystals to match the (1-init_solition) (w/w) seed
        def objective(init_crystals):
            temp_F = initial_F/sum(initial_F)*init_crystals
            initial_F[initial_F<0] = 0  # Cannot have negative sized crystals
            temp_weight = sys_funcs.moments(temp_F, 3, initial_boundaries, initial_y)*volume_shape_factor*rho_crystal
            return (temp_weight - (1-init_solution))**2

        # Minimize the function
        result = minimize_scalar(objective, method='brent')

        # Update initial CSD with optimized number of crystals
        opt_x, opt_y = result['x'], result['fun']
        initial_F = initial_F/sum(initial_F)*round(opt_x)
        initial_F[initial_F<0] = 0  # Cannot have negative sized crystals

        # Print out optimization results
        print('Number of initial seed crystals: %.6f' % round(opt_x))
        print('Initial seed mass: %.6f' % (opt_y + (1-init_solution)))
        print('Total Evaluations n: %d' % result['nfev'])

        
        # Initial crystallizer temperature
        cooling_profile = temp_generation(time_series)
        T_cry_init = cooling_profile[0]

        X0 = np.append(initial_F, (init_concentration, T_cry_init))
        t_span = (0.0, time_to_solve)
        signal.signal(signal.SIGALRM, signal_handler)
        signal.alarm(30)  # set a timer for 200 seconds
        try:
            result_solve_ivp = solve_ivp(differential_function, t_span, X0,
                                        t_eval=time_series, method='BDF',
                                        args=(initial_y, cooling_profile,
                                            growth_params[pair], nucleation_params[pair]))
        except TimeoutError:
            print(f"TimeoutError: solve_ivp timed out for Pair {pair+1} - Run {run+1}")
            return
        finally:
            signal.alarm(0)  # cancel the timer


        X_out = result_solve_ivp.y  # Solver output
        conc_profiles = X_out[-2, :]
        temp_profiles = X_out[-1, :]
        final_dists = X_out[0:-2, :]  # Temporal CSD
        final_times = time_series
        
        # Writing data to CSV or .MAT
        M_T_array = sys_funcs.moments_array(final_dists[:, :], 3, initial_boundaries, initial_y) * volume_shape_factor * rho_crystal
        SMD = sys_funcs.moments_array(final_dists[:, :], 4, initial_boundaries, initial_y) / (
                    sys_funcs.moments_array(final_dists[:, :], 3, initial_boundaries, initial_y))
        
        moments_0 = sys_funcs.moments_array(final_dists[:, :], 0, initial_boundaries, initial_y)
        moments_1 = sys_funcs.moments_array(final_dists[:, :], 1, initial_boundaries, initial_y)
        moments_2 = sys_funcs.moments_array(final_dists[:, :], 2, initial_boundaries, initial_y)
        moments_3 = M_T_array / volume_shape_factor / rho_crystal
        suspension_density = M_T_array

        c_dextrose = conc_profiles/(1-M_T_array)
        supersat_array = c_dextrose/kinetics.solubility(temp_profiles+273.15)
        growth_array = kinetics.growth_rate(supersat_array,temp_profiles,SMD,growth_params[pair])
        nucleation_array = kinetics.nucleation_rate(supersat_array,M_T_array, nucleation_params[pair])

        folder_name = f"pair_{pair+20}"
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        file_name = os.path.join(folder_name, f"trial_{run}.csv")
        result_data = np.vstack((X_out,moments_0, moments_1, moments_2, moments_3, suspension_density, SMD, final_times, cooling_profile, supersat_array,growth_array,nucleation_array))
        try:
            np.savetxt(file_name, result_data, delimiter=',')
            print(f"File saved successfully: {file_name}")
        except Exception as e:
            print(f"Error saving file: {file_name}")
            print(e)

        #Printing simulation progress to the command window/console/log file. 
        print('##############################################################')
        print(f"Run {run+1} for Pair {pair+1} is done.")
        print(f"growth_params = {growth_params[pair]}")
        print(f"nucleation_params = {nucleation_params[pair]}")
        print('##############################################################')

    