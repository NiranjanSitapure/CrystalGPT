# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 14:29:43 2022

@author: niranjan_sitapure
"""

# %% Basic Functions required for GECRYS

import numpy as np

##########################################################
def kelvin_to_fahrenheit(kelvin):
    c = kelvin - 273
    f = c * 9 / 5 + 32 
    return f

def fahrenheit_to_kelvin(fahrenheit):

    f = fahrenheit
    c = (f - 32) * 5 / 9
    K = c +273
    return K

##########################################################
# The moments function calculates the j^th moment for a given distribution
# Use this if the input distribution is a one dimensional array
def moments(F,j,initial_boundaries,sample_y):
    Dy = np.diff(initial_boundaries); 
    Fmo = np.sum(F*sample_y**j*Dy) 
    return Fmo

#If the input distirbution is an array of different distirbutions then use moments_array
def moments_array(F,j,initial_boundaries,sample_y): 
    times = F.shape[1]
    F = np.transpose(F);
    Dy = np.array([np.diff(initial_boundaries)]*times);   
    sample_y = np.array([sample_y]*times); 
    Fmo = np.sum(F*sample_y**j*Dy,axis = 1);
    return Fmo

############################################################
# Normal distribution function
def normal_dist(x , mean , sd):
    prob_density = 1/(sd * np.sqrt(2 * np.pi)) *np.exp( - (x - mean)**2 / (2 * sd**2) )
    return prob_density
