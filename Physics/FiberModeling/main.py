
"""

Written by Robert Wolle for Raytum Photonics

"""
import time
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']=200 # setting the dpi of outputted graphs.
from numba import njit
from scipy import sparse


"""
Light travels as a plane wave. It has a complex amplitude described by: 
    A(z) = A(0)e^ikz, with z being the propagation direction
    k = k_r + ik_i: the real and imaginary parts of the wavenumber
    The real part k_r ends up being imaginary in the exp(), meaning that it contributes a phase.
    Thus, k_r = beta, the phase constant.
    The imginary part ends up being a negative real number: an exponential decay
    Thus, k_i = alpha/2, with alpha being the absorption coefficient.


    Time evolution of the complex amplitude can be given by:
    e^(-iwt)
    w = omega = 2pi*f

    The propagation constant k depends on the optical frequency of light: k(w) = 2pi*k(f).
    Optical intensity (power per unit area) is proportional to |A(z)|^2

"""
alpha = 1                       # absorption coefficient

beta = 1                        # phase constant

n_core = 1.44                   # index of refraction
n_cladding = 1.45

lambda_vacuum = 1.5*10**(-3)    # 1500 nm wavelength

a = 10*10**(-6)                 # 10 um radius of fiber core

NA = np.sqrt(n_core**2 - n_cladding**2)

V = 2*np.pi()/lambda_vacuum*a*NA




