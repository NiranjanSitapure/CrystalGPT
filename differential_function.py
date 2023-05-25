

# %% Calculating differential for the PBE solver
import numpy as np
import crystallization_kinetics as kinetics
import system_fuctions
rho_crystal = 1540  #Density for Lysozyme crystals
volume_shape_factor = np.pi/6 #For tetragonal crystals


def differential_function(t,X,initial_boundaries,sample_y):
    y = sample_y;
    Dy = np.diff(initial_boundaries);
    ya = np.hstack((y[1:None],y[-1]))
    yb = np.hstack(( y[0],y[0:-2], y[-2]))
    
    
    F = X[0:-1];#CSD at time 't'
    c = X[-1]; #solute concentration 
    T = kinetics.Tprofile(t);
    c_dextrose = c/0.875 #converting kg_dext/kg_slurry to kg_dext/kg_solution
    S = c_dextrose/kinetics.solubility(T); #Solubility and saturation is measured per kg_solution
    
    # Growth Terms
    G = kinetics.growth_rate(S,T,y);
    
    Ga = kinetics.growth_rate(S,T,ya);
    Fa = np.hstack((F[1:None], F[-1]))
    
    Gb = kinetics.growth_rate(S,T,yb);
    Fb = np.hstack((F[0],F[0:-2], F[-2]))
  
    # Nucleation Term
    J = kinetics.nucleation_rate(c,S,T,F,initial_boundaries,sample_y); #This comes from the nucleation rate equation
    Fb[0] = J/G[0] #Nucleation boundary condition 
    
    
    # Growth Derivative
    dF = - ( Ga*Fa - Gb*Fb )/(ya-yb) ;
    dF[np.isnan(dF)] = 0;
    
    #Concentration Derivative
    dc = -3*volume_shape_factor*rho_crystal*np.sum(G*F*y**2*Dy)-J*y[0]**3*volume_shape_factor*rho_crystal; #UNITS: kg_dextros/kg_solution
    # dc = -3*volume_shape_factor*rho_crystal*np.sum(G*F*y**2*Dy)
    #Combine the differentials in one array
    dXdt = np.hstack((dF,dc))
    
    return dXdt