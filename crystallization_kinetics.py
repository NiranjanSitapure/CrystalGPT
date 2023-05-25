# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 14:45:40 2022

@author: niranjan_sitapure
"""



#Growth and nucleation kinetics can be changed based on your system

#Growth kinetics can also have a crystal size dependence.
import numpy as np
import system_fuctions

rho_crystal = 1540  #Density for Lysozyme crystals
volume_shape_factor = np.pi/6 #For tetragonal crystals

# Here S = C/C_sat; This is the saturation ratio


def growth_rate(S,T,y): # UNITS: in m/s    
    G = (S>1)*1.14*10**-3*np.exp(-29549/8.314/T)*(S-1)**1.05*np.ones(len(y));
    return G

# Heterogenous nucleation rate for dextrose
def nucleation_rate(conc,S,T,F,initial_boundaries,initial_y): # UNITS: #/s per kg_solution
    # T can be a np array
    M_T = system_fuctions.moments(F,3,initial_boundaries,initial_y)*volume_shape_factor*rho_crystal
    B= (S>1)*4.50*10**4*M_T**0.49*(S-1)**1.41
    
    
    return B
     
# Modify solubulity based on the system    
def solubility(T):
    # T can be a np array 
    global rhoc
    T = T- 273 # converting to Celsius
    c_sat = (0.804*T+29.94)*0.01  #UNITS: c_sat is in kg_dextrose/kg_solution
    return c_sat


# Linear temperature profile
def Tprofile(t): # Output UNITS: Kelvin
    # T = 42-0.375*(t/3600)   # Linear profile with temperature decreasing from 42C to 33 C in 24 h 
    T   = 42-0.375*(t/3600)+ 2*np.sin(t/3600) #Sinusoidal Temperature Profile
    T =T+273
    return T