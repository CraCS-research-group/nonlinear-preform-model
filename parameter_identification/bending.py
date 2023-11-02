# -*- coding: utf-8 -*-
"""
This script calculates parameters for the multi-linear asymmetric constitutive 
law, which represents the non-constant bending stiffness of fabrics based on 
experimental moment-curvature data. It uses numerical optimisation techniques 
to fit the model to the experimental data and provides the compressive 
stiffness, softening stiffness, and softening strain parameters.

This was developed for the paper: 
[1] Broberg PH, Lindgaard E, Thompson AJ, Belnoue JP-H, Hallett SR, 
    and Bak BLV (2003). "An accurate forming model for capturing the nonlinear 
    material behaviour of multilayered binder-stabilised fabrics and predicting
    fibre wrinkling" (Manuscript submitted)
    
For further explanation of the functions presented here, please take a look at
the reference.
    
Please feel free to use and adapt the codes but remember to give proper 
attribution.  

Input parameters:
-----------
t : float
    Thickness of the ply [m].

E_t : float
    Tensile modulus of the fabric [Pa].
    
test_name : str
    name of the csv file containing the experimental moment-curvature relation.

Returns:
--------
Compressive stiffness, softening stiffness, and softening strain printed.

Usage:
------
1. Ensure 'moment_curv.csv' file containing experimental moment-curvature data 
   is present in the 'characterised_properties' directory.
2. Adjust the parameters 't' and 'E_t' according to the fabric properties.
3. Run the script to obtain the compressive stiffness (E_c), softening 
   stiffness (E_s), and softening strain (esoft) parameters.
4. The script also generates a plot comparing the experimental data with the 
   model prediction and prints the obtained parameter values.

@authors: Peter H Broberg, Esben Lindgaard, Adam J Thompson, 
          Jonathan P-H Belnoue, Stephen R Hallett, Brian LV Bak
          Aalborg University/University of Bristol

Date: November 2023
"""

import scipy.optimize as op
import numpy as np
import matplotlib.pyplot as plt
import os
from module import *

####### Find parameters given a measure moment-curvature relationship ########
# Parameters to input for computation of the model
t   = 0.001   # Thickness of the ply
E_t = 1.407e9 # Tensile modulus Pa 
test_name = 'moment_curv.csv' # Name of the file containing the characterised curve

# Fitting parameters - Can be changed to obtain a better fit
start_fit    = 90   # Entry to start fitting the moment-curvature relation
scale_stiff  = 1e6  # Scale on the stiffness
scale_strain = 1e-5 # Scale on the strain
x0 = [1.0, 1.0]     # Starting guess

# Load experimental curve
cwd = os.getcwd()
characterised_properties_path = os.path.join(cwd, 'characterised_properties')
test_path = os.path.join(characterised_properties_path, test_name)

data = np.genfromtxt(test_path, delimiter=',')
data = np.transpose(data)
curvature_exp = data[0] 
moment_exp = data[1]

# Find the compressive modulus
B =  moment_exp[1] / curvature_exp[1] # Bending stiffness for the linear part 
                                      # of moment-curvature relation
E_c = find_compressive_mod(E_t, t, B)

# Optimize parameters
para = op.minimize(ls_min,x0, args=(moment_exp[start_fit:],                   \
                                    curvature_exp[start_fit:], E_t, E_c, t,   \
                                    scale_stiff, scale_strain))

E_s = para.x[0]*scale_stiff
esoft = para.x[1]*scale_strain
ksoft = calc_ksoft(E_c, E_t, t, esoft)

##### Plot the moment-curvature relationship and print values #################
# Plot parameters
n_plot = 100
max_k = 23

ksoft = calc_ksoft(E_c, E_t, t, esoft)
M_plot = np.zeros(n_plot)
k = np.linspace(ksoft, max_k, n_plot)

for i in range(len(k)):
    M_plot[i] = M_nl(k[i],E_t,E_c,E_s,esoft,t)

M_plot = np.insert(M_plot, 0, 0) # To plot the linear part of moment-curvature
k = np.insert(k, 0, 0) # To plot the linear part of moment-curvature

plt.figure()
plt.plot(curvature_exp, moment_exp, color = 'tab:orange',                     \
         label='Characterised behaviour')
plt.plot(k, M_plot, '--k', label = 'Model behaviour')
plt.xlabel(r'Curvature, $\kappa$ [1/m]')
plt.ylabel(r'Moment, $M$ [N m/m]')
plt.legend()

# Print values 
print('Parameters obtained from the data fitting: ')
print('Compressive stiffness: ' + str(E_c))
print('Softening stiffness: ' + str(E_s))
print('Softening strain: ' + str(esoft))



