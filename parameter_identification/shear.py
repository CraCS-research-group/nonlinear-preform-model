# -*- coding: utf-8 -*-
"""
This script computes and saves material laws for Mode I and Mode II cohesive 
interfaces to be used in ABAQUS/Explicit built-in laws based on given input 
data. The Mode I cohesive interface is modeled using a bilinear 
traction-separation law, while the Mode II interface uses experimental 
transverse shear stress-strain data and is adjusted for use as a cohesive law. 
The resulting tabular damage parameters are saved as a CSV file.

This was developed for the paper: 
[1] Broberg PH, Lindgaard E, Thompson AJ, Belnoue JP-H, Hallett SR, 
    and Bak BLV (2003). "An accurate forming model for capturing the nonlinear 
    material behaviour of multilayered binder-stabilised fabrics and predicting
    fibre wrinkling" (Manuscript submitted)
    
For further explanation of the functions presented here, please take a look at
the reference.
    
Please feel free to use and adapt the codes but remember to give proper 
attribution.  

@authors: Peter H Broberg, Esben Lindgaard, Adam J Thompson,  
          Jonathan Belnoue, Stephen R Hallett, Brian LV Bak
          Aalborg University/University of Bristol

Date: November 2023
"""

import numpy as np 
import matplotlib.pyplot as plt
import os
from module import *

h    = 0.001 # Thickness of layer

# Mode-I parameters 
tau_m = 60000
Gc    = 200
KI    = 5.0e9

# Mode-II parameters 
test_name = 'shear_stress_strain.csv'
slip_end = 3e-3
lin_start = 3e-4
n_lin = 100
C_shear = 1.0
delta_t = 1e-5

########################## Cohesive interface ################################
#### Mode-I ####
## Original parameters
d_0I  = tau_m / KI 
d_f   = 2*Gc / tau_m
n     = 500

## Plot the desired traction-separation law ###
d = np.linspace(0, d_f + d_0I, n)
t = np.zeros(n)
for i in range(n):
    t[i] = f_bilinear(d[i], KI, d_0I, d_f)

plt.figure()
plt.plot(d, t*1e-6)
plt.xlabel('Normal separation, $\Delta_n$  [m]')
plt.ylabel('Normal traction, $t_n$ [MPa]')
plt.legend(['Mode-I cohesive law'])

## Calculate damage as a function of separations
de = np.linspace(0, 1, n) # This is added for more points where changes are large
d_e = de*de * (d_f - d_0I)

D = np.zeros(n)

for i in range(n):
    t_bar = KI * (d_e[i] + d_0I)
    t     = f_bilinear(d_e[i] + d_0I, KI, d_0I, d_f)
    D[i]  = 1 - (t / t_bar)

##  Save the material law
mode = np.zeros(len(D))
input_array_0 = np.transpose(np.array([D, d_e, mode, np.zeros(len(D))]))
#np.savetxt('dam0.csv', input_array, delimiter=',')

##### Mode-II ####
## Original parameters
cwd = os.getcwd()
characterised_properties_path = os.path.join(cwd, 'characterised_properties')
test_path = os.path.join(characterised_properties_path, test_name)

my_data = np.genfromtxt(test_path, delimiter=',')
shear = np.transpose(my_data)

slip = make_stress_disp_curve(shear[0], h)
traction = shear[1] * 1e6


slip_adjusted, traction_adjusted = adjust_curve(slip, traction, C_shear, 
                                                delta_t, lin_start, n_lin, slip_end)

plt.figure()
plt.plot(slip_adjusted, traction_adjusted*1e-6)
plt.xlabel('Tansverse separation, $\Delta_s$ [m]')
plt.ylabel('Shear traction, $t_s$ [MPa]')
plt.legend(['Mode-II cohesive law'])

### Calculate damage as a function of separations
KII = initial_stiffness = traction_adjusted[1] / slip_adjusted[1]
d_0II = slip_adjusted[1]

D = np.zeros(len(traction_adjusted)-1)
d_e = slip_adjusted[1:] - slip_adjusted[1]

for i in range(len(traction_adjusted)-1):
    t_bar =  initial_stiffness * (d_e[i] + d_0II)
    t     = traction_adjusted[i+1]
    D[i]  = 1 - (t / t_bar)

mode = np.ones(len(D))
input_array_1 = np.transpose(np.array([D, d_e, mode, np.zeros(len(D))]))
input_array = np.concatenate((input_array_0,input_array_1))

model_parameters_path = os.path.join(cwd, 'model_parameters')
save_name = 'damage_table.csv'
save_path = os.path.join(model_parameters_path, save_name)

np.savetxt(save_path, input_array, delimiter=',')

# Print values 
print('Parameters obtained from the data fitting: ')
print('Mode-I stiffness, KI: ' + str(KI))
print('Mode-II stiffness, KII: ' + str(KII))
print('Mode-I separation at damge onset, d_0I: ' + str(d_0I))
print('Mode-II separation at damge onset, d_0II: ' + str(d_0II))


