
"""

Written by Robert Wolle for Raytum Photonics

"""
import time
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']=200 # setting the dpi of outputted graphs.
#from numba import njit
#from scipy import sparse


"""
WAVE AMPLITUDE

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

n_core = 1.45                   # index of refraction
n_cladding = 1.44

lambda_vacuum = 1.5*10**(-3)    # 1500 nm wavelength

a = 10*10**(-6)                 # 10 um radius of fiber core

A_0 = 1                         # Amplitude maximum

NA = np.sqrt(n_core**2 - n_cladding**2)     # Numerical Aperture

V = 2*np.pi/lambda_vacuum*a*NA            # 

nx = 101

x_min = 0
x_max = 2.2


hx = (x_max - x_min)/(nx-1)


x = np.linspace(x_min,x_max,nx)
A = np.zeros(nx)

for  i in range(nx):
    A[i] = A_0*np.exp(-alpha/2*x[i])

plt.plot(x,A)

plt.savefig('/workspaces/Projects/Physics/FiberModeling/output.png')