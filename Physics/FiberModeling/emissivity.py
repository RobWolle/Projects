
"""

Written by Robert Wolle for Raytum Photonics

"""
import time
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']=200 # setting the dpi of outputted graphs.


"""
Constants and conversion functions
"""
c = 299792458
h = 6.62607015*10**-34  # Planck Constant in Si
k_B = 1.380649*10**-23  # Boltzmann Constant in Si


def get_Kelvin(temp):
    # Input temperature in Celsius
    K = temp + 273.15
    return K

def get_lambda_n(f, n):
    v = c/n
    lambda_n = v/f
    return lambda_n

def get_BoltzmannLaw(wavelength,T):
    # This function returns the radiance of a black body at a given wavelength and temperature
    B = (2*h*c**2/wavelength**5)*1/(np.exp(h*c/(wavelength*k_B*T))-1)
    return B



wavelength_min = 0
wavelength_max = 40*10**-6
wavelength_resolution = 101

# Bookkeeping and avoiding dividing by 0:
lx = (wavelength_max - wavelength_min)/(wavelength_resolution-1)
wavelength_min = wavelength_min+lx
wavelength_max = wavelength_max+lx
wavelengths = np.linspace(wavelength_min,wavelength_max, wavelength_resolution)


placeholder_distribution = np.zeros(wavelength_resolution)
for k in range(wavelength_resolution):
    x = wavelengths[k]*10**6
    f1 = 0.7*np.cos(0.2*x+0.8)
    f2 = 0.4*np.cos(0.5*x+0.2)
    f3 = 1.4*np.cos(0.07*x+1.7)
    placeholder_distribution[k] = f1**2 + f2**2 + f3**2

n_lambda = placeholder_distribution

k_lambda = np.zeros(wavelength_resolution)

for k in range(wavelength_resolution):
    integral = 0
    for s in range(wavelength_resolution):
        if s==k:
            y = (n_lambda[s]-1)/(1-(wavelengths[s]/wavelengths[k])**2+lx**2)
        else:
            y = (n_lambda[s]-1)/(1-(wavelengths[s]/wavelengths[k])**2)
        integral = integral + y*lx
    P = 1/c
    k_lambda[k] = -2*integral/(np.pi*wavelengths[s])*P
    print(k_lambda[k])
    

"""
value of k at lambda = integral through all wavelengths = sum of n(x)dx/(x-lambda)
"""

subportion_k = round(wavelength_resolution*.7)
plt.plot(wavelengths[:subportion_k],n_lambda[:subportion_k])
plt.plot(wavelengths[:subportion_k],k_lambda[:subportion_k])

plt.savefig('/workspaces/Projects/Physics/FiberModeling/outputEmissivity.png')
plt.clf()

