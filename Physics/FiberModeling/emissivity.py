"""

Written by Robert Wolle for Raytum Photonics

"""
import time
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']=200 # setting the dpi of outputted graphs.

"""
Printing and utility functions
"""
def print_time(section_name, begin_time,stop_time, decimals):
    print(section_name + ' took ' + str(round(stop_time-begin_time,decimals)) + ' s')

"""
Constants and conversion functions
"""
c = 299792458
h = 6.62607015*10**-34  # Planck Constant in Si
k_B = 1.380649*10**-23  # Boltzmann Constant in Si


# k = 2pi/lambda, lamda changes in a material if f stays constant
# f = v/lambda
# v = c/n

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

# Intensity loss:
def get_I_absorption(k_n, lambda_vac, x_i):  # This is a different k (extinction coeff)
    alpha = 4*np.pi*k_n/lambda_vac
    I_absorption = np.exp(-alpha*x_i)
    return I_absorption

"""
Material Constants
"""
#Initialize time:
start = time.time()

Air_n = 1.00027417
Air_k = 0.005

Al2O3_1um_n = 1.6171
Al2O3_1um_k = 0

TiN_1um_n = 1.7547
TiN_1um_k = 3.4512
#TiN_1um_k = 0.4512


# Define regions of material by defining only:
# [0] = their widths <m> (to then call them in order to fill the space)
# [1] = their index of refraction n
# [2] = their extinction coefficient k

Domains = []

Air_L = .1*10**-6
Air = [Air_L, Air_n, Air_k]
Domains.append(Air)

Film1_L = .08*10**-6
Film1 = [Film1_L, Al2O3_1um_n, Al2O3_1um_k]

Film2_L = .02*10**-6
Film2 = [Film2_L, TiN_1um_n, TiN_1um_k]

layers = 5

for i in range (layers):
    Domains.append(Film1)
    Domains.append(Film2)


amplitude_time = time.time()
print_time("Initializing", start, amplitude_time, 4)

"""
Input Conditions
"""


Amplitude = 1               # Constant amplitude of intitial signal. 
lambda_vacuum = 1*10**-6    # Desired vacuum wavelength of light.
resolution = 1001                    # Resolution for each domain.


def E_at_lambda(Domains, Amplitude, lambda_vacuum, resolution):
    nx = resolution
    E_0 = Amplitude
    f = c/lambda_vacuum

    x_all = np.zeros(1)
    E_all = np.zeros(1)
    I_all = np.zeros(1)
    x_starts = [0]
    phi_starts = [0]
    d=0
    for domain in Domains:
        lambda_n = get_lambda_n(f, domain[1])
        x_min = 0
        x_max = domain[0]
        x = np.linspace(x_min,x_max, nx)

        hx = (x_max - x_min)/(nx-1)
        x_starts.append(x_max+x_starts[d]) # sets the start of the next values of x in the x_all list
        
        x_all = np.append(x_all, x+x_starts[d])
        
        E = np.zeros(nx)
        for  i in range(nx):
            E_vac = E_0*np.cos(2*np.pi/lambda_n*x[i]+2*np.pi*phi_starts[d])         # E(x) = E_0.cos(kx) = sqrt(I_0)
            E[i] = E_vac*np.sqrt(get_I_absorption(domain[2], lambda_vacuum, x[i]))  # I(x) = I_0.e^(-alpha.x)
        
        
        E_all = np.append(E_all, E)
        E_0 = E_0*np.sqrt(get_I_absorption(domain[2], lambda_vacuum, x[nx-1]))
        phi_starts.append(phi_starts[d]+x_max/lambda_n)
        
        d=d+1

    E_all = np.delete(E_all, 0)
    I_all = E_all**2
    x_all = np.delete(x_all, 0)
    return E_all, I_all, x_all


"""
Wavelength-dependent intensity input parameters
"""
Wavelength_dependent = True         # Set to false if you want intensity at just one input wavelength
Wavelength_single = 1*10**-6        # Set to the desired input wavelength if you only want intensity at a single wavelength.

Thermal_source = True               # Set to false if you do not want to approximate a black body light source
#^ This may need to be switched to a selector for several standard wavelength-dependent light sources
Temperature = get_Kelvin(300)       # temperature in Kelvin

# Input wavelength range and resolution desired:
wavelength_min = 0
wavelength_max = 10*10**-6
wavelength_resolution = 101

# Bookkeeping and avoiding dividing by 0:
lx = (wavelength_max - wavelength_min)/(wavelength_resolution-1)
wavelength_min = wavelength_min+lx
wavelength_max = wavelength_max+lx
wavelengths = np.linspace(wavelength_min,wavelength_max, wavelength_resolution)
Amplitudes = np.zeros(wavelength_resolution)+Amplitude

# Setting Amplitude over wavelength range:
if Thermal_source == True:
    for j in range(wavelength_resolution):
        Amplitudes[j] = get_BoltzmannLaw(wavelengths[j], Temperature)
#^ This is not as efficient as calling this in the same loop that calculates the E for each wavelength below.
#^ It instead makes a single conditional check instead of once per loop.

wavelengthdependent_time = time.time()
print_time("Calculating amplitudes", amplitude_time,wavelengthdependent_time, 4)

# For single wavelengths:
if Wavelength_dependent == False:
    wavelengths[0] = Wavelength_single

# Performing the calculation once to initiate the arrays:
lambda_vacuum = wavelengths[0]
E_total, I_total, x_total = E_at_lambda(Domains, Amplitudes[0], lambda_vacuum, resolution)

# Iterating for each wavelength in range:
if Wavelength_dependent == True:
    for j in range(1, wavelength_resolution):
        E_j, I_j, dump = E_at_lambda(Domains, Amplitudes[j], wavelengths[j], resolution)
        E_total = E_total + E_j
        I_total = I_total + I_j

plotting_time = time.time()
print_time("Calculating wavelength-dependent intensities", wavelengthdependent_time,plotting_time, 4)


plt.plot(x_total,E_total)
plt.savefig('/workspaces/Projects/Physics/FiberModeling/outputE.png')
plt.clf()

plt.plot(x_total,I_total)
plt.savefig('/workspaces/Projects/Physics/FiberModeling/outputI.png')
plt.clf()

if Thermal_source == True:
    plt.plot(wavelengths, Amplitudes)
    plt.savefig('/workspaces/Projects/Physics/FiberModeling/outputThermalSource.png')

end_time = time.time()
print_time("Plotting", plotting_time, end_time, 4)

