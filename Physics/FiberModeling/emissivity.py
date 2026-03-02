
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


"""
This next section attempts to calculate the Kramers-Kronig relation for k(lambda).
k(lambda) = -2.lambda/pi * int_0^inf[(n(lambda')-1)/(lambda^2 - lambda'^2)dlambda']
    The residue of the pole at lambda = lambda' is -[n(lambda) - 1]/(2.lambda)
This means k(lambda) = [n(lambda) - 1]/pi * int_0^inf[(n(lambda')-1)/(lambda^2 - lambda'^2)dlambda']
This isn't normalized properly. I think this is because in the conversion from frequency to wavelength, I'm missing a factor of c
    It's also not normalized because n is only known on a small region. 
        n --> 1 at lambda = infinity, so n - 1 should converge to 0.
        This fact might be aided by something like an exponential decay. 
        Depending on the tail of the decay, it seems like the integral should converge to a number around 100, which contributes a larger integral factor
"""

k_lambda = np.zeros(wavelength_resolution)

for k in range(wavelength_resolution):
    integral = 0
    for s in range(wavelength_resolution):
        if s==k:
            y = (n_lambda[s]-1)/(1-(wavelengths[s]/wavelengths[k])**2+lx**2)
        else:
            y = (n_lambda[s]-1)/(1-(wavelengths[s]/wavelengths[k])**2)
        integral = integral + y*lx
    Correction = 200/c
    k_lambda[k] = (n_lambda[k]-1)/np.pi*integral*Correction       
    
    

"""
value of k at lambda = integral through all wavelengths = sum of n(x)dx/(x-lambda)
"""

subportion_k = round(wavelength_resolution*.7)
plt.plot(wavelengths[:subportion_k],n_lambda[:subportion_k])
plt.plot(wavelengths[:subportion_k],k_lambda[:subportion_k])

plt.plot(wavelengths,n_lambda)
plt.plot(wavelengths,k_lambda)

plt.savefig('/workspaces/Projects/Physics/FiberModeling/outputEmissivity.png')
plt.clf()


"""
Calculating Reflectance and Transmission
"""

def get_R(n, k):
    R =  ((n-1)**2 + k**2)/((n+1)**2+k**2)
    return R


def get_T(k, wavelengths_array, thickness):
    T = np.exp(-4*np.pi*k*thickness/wavelengths_array)
    return T

def get_emissivity(T, R):
    emissivity = (1-R)*(1-T)/(1-R*T)
    return emissivity


thickness = 0.1*10**-6

R_lambda = get_R(n_lambda,k_lambda)
T_lambda = get_T(k_lambda, wavelengths, thickness)

emissivity_lambda = get_emissivity(T_lambda,R_lambda)

plt.plot(wavelengths,emissivity_lambda)

plt.savefig('/workspaces/Projects/Physics/FiberModeling/outputEmissivity2.png')
plt.clf()
