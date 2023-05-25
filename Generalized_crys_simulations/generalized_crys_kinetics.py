# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 14:45:40 2022

@author: niranjan_sitapure
"""



#Growth and nucleation kinetics can be changed based on your system

#Growth kinetics can also have a crystal size dependence.
import numpy as np
import random

# Here S = C/C_sat; This is the saturation ratio

#%% 
# Growth Rate
per_deviaition = 0.05

# def growth_rate(supersat, temp, sizes, params):
        
#     """
#     This function takes in supersaturation (S), temperature (T), sizes (y), and parameters [a_g, b_g, c_g],
#     and returns a result based on those inputs.
#     """
#     # unpack parameters
#     a_g, b_g, c_g = params

#     if supersat<1:
#         G=0.0*np.ones(len(sizes))
#     else:
#         G = a_g*np.exp(b_g/8.314/temp)*(supersat-1)**c_g*np.ones(len(sizes))
#     return G

def growth_rate(supersat, temp, sizes, params):
    """
    This function takes in supersaturation (S), temperature (T), sizes (y), and parameters [a_g, b_g, c_g],
    and returns a list of growth rates based on those inputs.
    """
    # unpack parameters
    a_g, b_g, c_g = params

    # Check if inputs are scalar or vector
    if np.isscalar(supersat) and np.isscalar(temp):
        if supersat < 1:
            G = np.zeros(len(sizes))
        else:
            G = a_g * np.exp(b_g / 8.314 / temp) * (supersat - 1) ** c_g*np.ones(len(sizes))
    else:
        # Create an array of zeros for supersat < 1
        G = np.zeros(len(sizes))
        # Compute growth rates for supersat >= 1
        mask = (supersat >= 1)
        G[mask] = a_g * np.exp(b_g / 8.314 / temp[mask]) * (supersat[mask] - 1) ** c_g
    
    return np.array(G)




import random

# Generate Random Parameters for Growth Rate
def random_growth_generator(a_mean = 1.14*10**-3, b_mean = -29549, c_mean = 1.05):
    """
    This function takes in the mean values of a, b, and c, and generates a random number (sigma)
    around those mean values.
    """
    # calculate standard deviations for a, b, and c
    a_std = per_deviaition * a_mean
    b_std = per_deviaition * abs(b_mean)
    c_std = per_deviaition * c_mean
    
    # generate random numbers for a, b, and c
    a = random.normalvariate(a_mean, a_std)
    b = random.normalvariate(b_mean, b_std)
    c = random.normalvariate(c_mean, c_std)
    
    # return the random numberss
    return a, b, c


#%% 

# Growth Rate
per_deviaition = 0.1

def nucleation_rate(supersat, sus_density, params):
        
    """
    This function takes in supersaturation (S), suspension density (M_T), and parameters [a_g, b_g, c_g],
    and returns a result based on those inputs.
    """
    # unpack parameters
    a_B, b_B, c_B = params

# Check if inputs are scalar or vector
    if np.isscalar(supersat) and np.isscalar(sus_density):
        if supersat < 1:
            B = 0.0
        else:
            B = a_B*sus_density**b_B*(supersat-1)**c_B
    else:
        # Create an array of zeros for supersat < 1
        B = np.zeros(len(supersat))
        # Compute growth rates for supersat >= 1
        mask = (supersat >= 1)
        B[mask] = a_B*sus_density[mask]**b_B*(supersat[mask]-1)**c_B
    
    return np.array(B)

# Generate Random Parameters for Growth Rate
def random_nucleation_generator(a_mean = 4.50*10**4, b_mean = 0.49, c_mean = 1.41):
    """
    This function takes in the mean values of a, b, and c, and generates a random number (sigma)
    around those mean values.
    """
    # calculate standard deviations for a, b, and c
    a_std = per_deviaition * a_mean
    b_std = per_deviaition * abs(b_mean)
    c_std = per_deviaition * c_mean
    
    # generate random numbers for a, b, and c
    a = random.normalvariate(a_mean, a_std)
    b = random.normalvariate(b_mean, b_std)
    c = random.normalvariate(c_mean, c_std)
    
    # return the random numbers
    return a, b, c

     

#%% 

# Modify solubulity based on the system    
def solubility(T):
    # T can be a np array 
    global rhoc
    T = T- 273 # converting to Celsius
    c_sat = (0.804*T+10)*0.01  #UNITS: c_sat is in kg_dextrose/kg_solution
    # c_sat = 1.0036*10**(-3)*T**3 + 1.4*10**(-2)*T**2 -0.128*T + 3.4613 
    return c_sat


# Linear temperature profile
# def Tprofile(t): # Output UNITS: Kelvin
#     if (t%(30*60)):
#         #change the slope every 30 minutes
#         slope = random.uniform(0.1, 1)
    
#     T = 42-slope*(t/3600)   # Linear profile with temperature decreasing from 42C to 33 C in 24 h 
#     # T   = 42-0.375*(t/3600)+ 2*np.sin(t/3600) #Sinusoidal Temperature Profile
#     T =T+273
#     return T