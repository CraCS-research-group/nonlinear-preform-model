# -*- coding: utf-8 -*-
"""
This module contains the functions used for parameter identification for the 
nonlinear bending model and the nonlinear transverse shear model developed 
for the paper: 
[1] Broberg PH, Lindgaard E, Thompson AJ, Belnoue JP-H, Hallett SR, 
    and Bak BLV (2003). "An accurate forming model for capturing the nonlinear 
    material behaviour of multilayered binder-stabilised fabrics and predicting
    fibre wrinkling" (Manuscript submitted)
    
For further explanation of the functions presented here, please take a look at
the reference.

Feel free to use and adapt the codes but remember to give proper attribution.  

@authors: Peter H Broberg, Esben Lindgaard, Adam J Thompson, 
          Jonathan P-H Belnoue, Stephen R Hallett, Brian LV Bak
          Aalborg University/University of Bristol

Date: November 2023
"""

import scipy.optimize as op
import numpy as np
from scipy.integrate import quad

def calc_z0(E_c, E_t, t):
    """
    Calculate the neutral line of the layer before softening based on the given 
    compressive modulus (E_c), tensile modulus (E_t), and thickness (t).

    Parameters:
    -----------
    E_c : float
        Compressive modulus of the material [Pa].
    
    E_t : float
        Tensile modulus of the material [Pa].
    
    t : float
        Thickness of the ply [m].
    
    Returns:
    --------
    float
        Distance of the neutral line before softening [m].

    References:
    -----------
    Equation (6) in [1]
    """
    z0 = t * (E_c + E_t - 2*np.sqrt(E_c*E_t))/(2*(E_c - E_t))
    return z0

def func_diff_b(E_c, E_t, t, B):
    """
    Calculate the function value used to find the compressive modulus (E_c) 
    when the constant bending stiffness (B), tensile modulus (E_t), and 
    thickness of the ply (t) is known.

    Parameters:
    -----------
    E_c : float
        Compressive modulus of the material [Pa].
    
    E_t : float
        Tensile modulus of the material [Pa].
    
    t : float
        Thickness of the ply [m].

    B : float
        Known linear bending stiffness [Nm].
    
    Returns:
    --------
    float
        The difference between the calculated bending stiffness and the known 
        linear bending stiffness [Nm].

    References:
    -----------
    Equation (7) in [1]
    """
    z0 = calc_z0(E_c, E_t, t)
    Bc = 1./24 * (E_t*(t**3 + 3*t**2*z0 - 4*z0**3) + E_c*(t**3 - 3*t**2*z0 +   \
                                                          4*z0**3)) 
    return Bc - B 

def calc_ksoft(E_c, E_t, t, esoft):
    """
    Calculate the softening curvature given the compressive modulus (E_c), 
    tensile modulus (E_t), thickness of the ply (t), and softening strain 
    (esoft).

    Parameters:
    -----------
    E_c : float
        Compressive modulus of the material [Pa].
    
    E_t : float
        Tensile modulus of the material [Pa].
    
    t : float
        Thickness of the ply [m].

    esoft : float
        Softening strain of the material [m].
    
    Returns:
    --------
    float
        Softening curvature [1/m].
    
    References:
    -----------
    Equation (5) in [1]
    """
    z0 = calc_z0(E_c, E_t, t)
    ksoft =  esoft / (t/2 - z0)
    return ksoft

def integrand(z, z0, E_t):
    """
    Compute the integrand function used in the non-linear moment calculation.

    Parameters:
    -----------
    z : float
        Position along the thickness of the ply [m].
    
    z0 : float
        Location of neutral axis before softening [m].
    
    E_t : float
        Tensile modulus of the material [Pa].
    
    Returns:
    --------
    float
        Value of the integrand function at the given position [Pa m].


    References:
    -----------
    Auxillary function to solve equation (2) in [1]

    """
    return E_t*z*(z-z0)

def M_nl(k, E_t, E_c, E_s, e_soft, t):
    """
    Calculate the non-linear moment given curvature based on the multi-linear 
    asymmetric constitutive law.

    Parameters:
    -----------
    k : float
        Curvature of the neutral axis [1/m].
    
    E_t : float
        Tensile modulus of the material [Pa].
    
    E_c : float
        Compressive modulus of the material [Pa].
    
    E_s : float
        Softening modulus of the material [Pa].
    
    e_soft : float
        Softening strain of the material.
    
    t : float
        Thickness of the ply [m].

    Returns:
    --------
    float
        Non-linear moment based on the given curvature [N m].

    References:
    -----------
    Equation (2) in [1]

    """
    z0,zp =   z0_and_zp(k, E_t, E_c, E_s, e_soft, t)  

    int1 = quad(integrand, -t/2, z0, args=(z0,E_t))[0]
    int2 = quad(integrand, z0, zp, args=(z0,E_c))[0]
    int3 = quad(integrand, zp, t/2, args=(z0,E_s))[0]
    
    M = k*(int1+int2+int3) + (E_c-E_s)*e_soft/8 *(t**2 - 4*zp**2)
    return M

def inverse_z0(z0, Et, Ec, Es, epsilons, t, k):
    """
    Calculate the force equilibrium of the cross-section and solve for the 
    neutral axis position (z0) using an inverse approach. 

    Parameters:
    -----------
    z0 : float
        Location of neutral line  [m].

    Et : float
        Tensile modulus of the material [Pa].
    
    Ec : float
        Compressive modulus of the material [Pa].
    
    Es : float
        Softening modulus of the material [Pa].

    epsilons : float
        Softening strain of the material.
    
    t : float
        Thickness of the ply [m].

    k : float
        Curvature of the layer [1/m].
    
    Returns:
    --------
    float
        Sum of forces in the cross-section [N].

    References:
    -----------
    Equation (4) in [1]
    
    """
    epsilont = (t/2 + z0)*k
    epsilonc = (t/2-z0)*k
    
    zp = epsilons/k + z0
    
    sigmat = Et*epsilont
    sigmac = Ec*epsilons
    sigmas = (epsilonc-epsilons)*Es + sigmac
    
    At = t/2 + z0
    Ac = zp - z0
    As = t/2 - zp 
    
    Ft = 0.5*sigmat*At
    Fc = 0.5*sigmac*Ac
    Fs = 0.5*(sigmac+sigmas)*As
    
    Fsum = Ft - Fc - Fs
    return Fsum 

def z0_and_zp(k, E_t, E_c, E_s, esoft, t):
    """
    Determine the neutral and softening line positions given the curvature,
    using an inverse approach.

    Parameters:
    -----------
    k : float
        Curvature of the neutral axis [1/m].
    
    E_t : float
        Tensile modulus of the material [Pa].
    
    E_c : float
        Compressive modulus of the material [Pa].
    
    E_s : float
        Softening modulus of the material [Pa].
    
    esoft : float
        Softening strain of the material.
    
    t : float
        Thickness of the ply [m].

    Returns:
    --------
    tuple
        Tuple containing two floats: the location of the neutral line
        (z0, [m]) and the location of the softening line (zp, [m]).

    References:
    -----------
    Equation (3) in [1]
    
    """
    z0_guess = [0.0]
    root = op.fsolve(inverse_z0,z0_guess, args=(E_t, E_c, E_s, esoft, t, k))
    z0 = root[0]
    zp = esoft/k + z0
    return z0, zp

def find_compressive_mod(E_t, t, B):
    """
    Determine the compressive modulus (E_c) given the linear bending stiffness 
    (B) of the material.

    Parameters:
    -----------
    E_t : float
        Tensile modulus of the material [Pa].
    
    t : float
        Thickness of the ply [m].
    
    B : float
        Known linear bending stiffness [Nm].

    Returns:
    --------
    float
        Compressive modulus of the material [Pa].

    References:
    -----------
    Solve equation (7) in [1]

    """
    x0 = [1.e9] # Starting guess
    root = op.fsolve(func_diff_b, x0, args=(E_t, t, B))
    return root[0]

def ls_min(x, moment_exp, curvature_exp, E_t, E_c, t, scale_stiff,            \
           scale_strain):
    """
    Least squares objective function for numerical optimisation to fit the 
    model to experimental moment-curvature data and find parameters.

    Parameters:
    -----------
    x : list
        List containing two parameters: [E_s, e_soft], representing the 
        softening modulus [Pa] and softening strain, respectively.
    
    moment_exp : array-like
        Experimental moment data [Nm]
    
    curvature_exp : array-like
        Experimental curvature data [1/m]
    
    E_t : float
        Tensile modulus of the material [Pa].
    
    E_c : float
        Compressive modulus of the material [Pa].
    
    t : float
        Thickness of the ply [m].
    
    scale_stiff : float
        Scaling factor for the softening modulus (dimensionless).
    
    scale_strain : float
        Scaling factor for the softening strain (dimensionless).

    Returns:
    --------
    float
        Objective function value representing the least squares difference 
        between model-predicted and experimental moments.

    References:
    -----------
    Solve equation (2) in [1]
    """
    M_bl = np.zeros(len(moment_exp))
    for i in range(len(M_bl)):
        M_bl[i] = M_nl(curvature_exp[i], E_t, E_c, x[0]*scale_stiff,          \
                       x[1]*scale_strain, t)
    
    diff = moment_exp - M_bl
    diff = np.einsum('...i,...i', diff, diff)
    
    return np.abs(diff)


def f_bilinear(d, K, d_0, d_f):
    """
    Calculate the traction (t) based on a bilinear traction-separation law.

    Parameters:
    -----------
    d : float
        Separation [m].
    
    K : float
        Initial stiffness of the interface [N/m].
    
    d_0 : float
        Displacement at which damage onset occurs [m].
    
    d_f : float
        Displacement at final failure [m].

    Returns:
    --------
    float
        Traction [m] based on the bilinear traction-separation law.
    """
    if d <= d_0:
        t = K*d
    if d > d_0:
        t_0 = K*d_0
        t = np.interp(d, [d_0, d_f], [t_0, 0])
    return t

def make_stress_disp_curve(gamma, t):
    """
    Convert shear strain angles to corresponding displacement values for a 
    material layer with a given thickness.

    Parameters:
    -----------
    gamma : float
        Shear strain angle [degrees].
    
    t : float
        Thickness of the material layer [m].

    Returns:
    --------
    float
        Displacement [m] calculated based on the given shear strain angle and 
        layer thickness.
    """
    gamma_rad = gamma / 360 * 2 * np.pi
    D = np.tan(gamma_rad)*t
    return D


def adjust_curve(slip, traction, C_shear, delta_t, lin_start, n_lin, slip_end):
    """
    Adjust the given shear traction and displacement curves to ensure 
    smoothness and proper behavior for cohesive material modelling. See [1] 
    Sec. 2.4 for detailed description of this

    Parameters:
    -----------
    slip : array-like
        Array containing the original transverse separation data points.

    traction : array-like
        Array containing the original shear traction data points.

    delta_t : float
        Displacement value used to create the linear portion of the curve.

    lin_start : float
        Start displacement for the linear portion of the curve.

    n_lin : int
        Number of points for the linear portion of the curve.

    slip_end : float
        End displacement for the curve.

    Returns:
    --------
    slip_adjusted : array-like
        Array containing adjusted transverse displacement data points.

    traction_adjusted : array-like
        Array containing adjusted shear traction data points.
    """
    #### ADD truncated PART to shear traction curve - to make it go to 0
    rem = 5
    
    slip_adjusted = np.delete(slip, [i for i in range(1,rem)])
    traction_adjusted = C_shear * np.delete(traction, [i for i in range(1,rem)])
    
    
    n_fit = 10
    
    sep_t = slip_adjusted[-1] - delta_t
    
    idx_t = np.where(slip_adjusted > sep_t)[0][0]
    
    x_fit = slip_adjusted[idx_t-int(n_fit/2):idx_t+int(n_fit/2)]
    y_fit = traction_adjusted[idx_t-int(n_fit/2):idx_t+int(n_fit/2)]
    
    p_fit = np.polyfit(x_fit, y_fit, 1)
    
    
    stiff = p_fit[0]
    slip_adjusted = slip_adjusted[:idx_t]
    traction_adjusted = traction_adjusted[:idx_t]
    
    inc_slip = np.linspace(slip_adjusted[-1], lin_start, int(n_lin/2))
    inc_trac = np.linspace(traction_adjusted[-1], traction_adjusted[-1] +     \
                           stiff*(lin_start-slip_adjusted[-1]), int(n_lin/2))
    inc_comp = np.linspace(0, 1e-5, int(n_lin/2)) # compensation 
    inc_trac = inc_trac - inc_comp
    
    slip_adjusted = np.append(slip_adjusted, inc_slip[1:])
    traction_adjusted = np.append(traction_adjusted, inc_trac[1:])
    
    lin_slip = np.linspace(slip_adjusted[-1], slip_end,  int(n_lin/2))
    lin_trac = np.linspace(traction_adjusted[-1], 0,  int(n_lin/2))
    
    slip_adjusted = np.append(slip_adjusted, lin_slip[1:])
    traction_adjusted = np.append(traction_adjusted, lin_trac[1:])
    return slip_adjusted, traction_adjusted 

