
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
    """v = c/n
    lambda_n = v/f
    return lambda_n"""
    return c/(f*n)

def get_BoltzmannLaw(wavelength,T):
    # This function returns the radiance of a black body at a given wavelength and temperature
    B = (2*h*c**2/wavelength**5)*1/(np.exp(h*c/(wavelength*k_B*T))-1)
    return B

def get_KramersKronigRelation_k_from_n(n_lambda,wavelengths, wavelength_resolution, correction):
    wavelength_resolution = len(n_lambda)
    lx = wavelengths[1]-wavelengths[0]

    k_lambda = np.zeros(wavelength_resolution)

    for k in range(wavelength_resolution):
        integral = 0
        for s in range(wavelength_resolution):
            if s==k:
                y = (n_lambda[s]-1)/(1-(wavelengths[s]/wavelengths[k])**2+lx**2)
            else:
                y = (n_lambda[s]-1)/(1-(wavelengths[s]/wavelengths[k])**2)
            integral = integral + y*lx
        C = correction/c
        k_lambda[k] = (n_lambda[k]-1)/np.pi*integral*C       


"""
AR Coatings Work

Task:

AR Coating at 980 nm +/- 10 nm with <0.2% reflection

get_R(n, k):
    R =  ((n-1)**2 + k**2)/((n+1)**2+k**2)

get_T(k, wavelengths_array, thickness):
    T = np.exp(-4*np.pi*k*thickness/wavelengths_array)

"""

def get_R(n1,n2):
    R = (n2-n1)/(n2+n1)
    return R**2


def getCoatingIndex_range(n1,n3, threshold):
    
    resolution = 100000
    x = np.linspace(0,5,resolution)
    Indexes = getSingleLayer_R(n1,x,n3)
    valid_ranges = []

    #default condition is that it starts in a range under the threshold.
    last_index = -5
    for i in range(resolution):
        if Indexes[i] < threshold:
            if last_index != i-1:
                range_start = i
            last_index = i
        else:
            if last_index == i-1:
                range_end = i-1
                min_value = min(Indexes[range_start:range_end])
                valid_ranges.append([x[range_start],x[range_end],min_value])
    if last_index == resolution-1:
        range_end = resolution-1
        min_value = min(Indexes[range_start:range_end])
        valid_ranges.append([x[range_start],x[range_end],min_value])
        valid_ranges.append([range_start,range_end])
    
    #dropping the 0 ranges:
    if len(valid_ranges)>0:
        if valid_ranges[0][1] < 0.5:
            valid_ranges.pop(0)

    if len(valid_ranges) == 1:
        print("The appropriate index of refraction of the coating layer falls within the range:")
        print(f" > {valid_ranges[0][0]:.3}-{valid_ranges[0][1]:.3}")
    elif len(valid_ranges) > 1:
        print("The appropriate index of refraction of the coating layer falls within these ranges:")
        for valid_range in valid_ranges:
                    print(f" > {valid_range[0]:.3}-{valid_range[1]:.3}")
    else:
        print("No valid ranges of index of refraction were found meeting that coating criteria.")

    return x, Indexes

def getSingleLayer_R(n1,x,n3,printing=False):
    R_1 = get_R(n1,x)
    R_2 = get_R(x,n3)
    Reflectivity = abs(R_1 - R_2)

    if printing==True:
        print(f"The reflectivity of the resulting coating is {Reflectivity:.3}.")
    return Reflectivity
#---------------------
# Constants at 980 nm

Air_n   = 1.00028526

# Films
Al2O3_n = 1.6670
AlN_n   = 2.1320
TiN_n   = 2.1429
TiO2_n  = 2.4880
ZnO_n   = 1.9444
HfO2_n  = 1.8822
MgO_n   = 1.7232
MoO3_n  = 2.0850
MgF2_n  = 1.3738
NaF_n   = 1.3215
CAF_n   = 1.4291
ZrO2_n  = 2.1256
GeO2_n  = 1.5958
Al_n    = 1.4682
CaCO3_n = 1.6442

#Substrates
K108_Glass_n = 1.5048
GaAs_n = 3.4824
Silica_n = 1.4507
Silicon_n = 3.5780

#--------------------

#Single Film Range for R<0.02
K108_Glass_n = 1.5048
GaAs_n = 3.4824
Silica_n = 1.4507
Silicon_range= [3.27167,3.91301]

#--------------------
incident_wavelength = 0.980*10**-6
quarter_thick = 0.25*incident_wavelength



threshold = 0.002
incident_layer = Air_n
substrate_layer = K108_Glass_n
x,y = getCoatingIndex_range(incident_layer,substrate_layer, threshold)

coating_layer = HfO2_n
getSingleLayer_R(incident_layer, coating_layer, substrate_layer, printing=True)

plt.plot(x,y)
plt.axhline(threshold,color='k')
plt.xlabel('index of refraction n')
plt.ylabel('Reflectivity R')
plt.savefig('/workspaces/Projects/Physics/FiberModeling/outputARcoating.png')


"""
To make this general, I would have to add the appropriate sine functions for the waves:

Modify the emissivityClasses file to use sin(k(2d-x)+pi).
I'd have to check the phase situation with multiple reflections

"""