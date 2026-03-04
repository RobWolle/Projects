
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


"""
Intensity making it to interface = T(k,lambda)
Proportion reflected back at interface = R(n,k,lambda)

Only thing that matters is total emissivity,
    which depends on intensity of light leaving and also how much light makes it through.

R,T depends on the properties of the material on the other side of a boundary.
Data structure:
Array of Intensity, propagation direction, and layer.

Only branches that matter are ones that make it to the surface or through,
    but for simulation purposes we should just do raw intensity at all positions.

Layer object tells you to add this to reflected total or transmitted total
    (in layer [0], out layer [max-1])
    This works as each object in the Domains structure, which also stores thickness and mat. properties

Still ends up being recursive:
    For each layer, call the "calculate E" function.
    If on T path:
        Call calc_E(E_0 = transmitted intensity)
        Iterate layer location.
        Break if layer location ends.
    If intensity > 0 threshold AND on R path:
        Call calc_E(E_0 = Reflected intensity)
        Iterate layer location
    

Calculate intensity transmitted continuously through the material. 
Once a boundary is reached, calculate the amount reflected.
Call the "calculate intensity continuously through the material" function with the reflected intensity.
    Once the reflected portion is done calculating, move on the portion transmitted.
        
"""




def make_DataArray(delimiter, data_columns, data_start, file_path):
    with open(file_path, "r") as file:
        data = file.read()

    data_values = []
    last_break = data_start
    for char in range(data_start,len(data)):
        if data[char] == delimiter or data[char] == '\n':
            if data[last_break:char] == '':
                data_values.append(0)
            else:
                data_values.append(float(data[last_break:char]))
            last_break = char+1

    # We can use the number of data points divided by number of columns to get the number of rows of data
    data_rows = int(len(data_values)/data_columns)

    data_organized = []
    for row in range(data_rows):
        data_organized.append(data_values[row*data_columns:(row+1)*data_columns])
    return data_organized


# User inputs:
delimiter = ','
data_columns = 6
data_start = 0    #character position at which the data begins
file_path = "/workspaces/Projects/Physics/FiberModeling/n and k fits - TiN and Al2O3.csv"
DataArray = np.array(make_DataArray(delimiter, data_columns, data_start, file_path))

wavelengths = DataArray[:,0]*10**-6

TiN_n = DataArray[:,1]
TiN_k = DataArray[:,2]

Al2O3_n = DataArray[:,4]
Al2O3_k = DataArray[:,5]



# Material n and k must be arrays which contain values at every wavelength in wavelengths[]

Air_n = np.ones(wavelength_resolution)*1.000273     # This is actually a measured thing that isn't exactly constant, but is almost constant in IR.
Air_k = np.ones(wavelength_resolution)*0.0005       # Technically this should be calculated through Kramers-Kronig, but it's close enough to zero.

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


Domains = []

Air_t = .1*10**-6
Air = [Air_t, Air_n, Air_k, get_R(Air_n, Air_k)]

#print(Air)

Film1_t = .08*10**-6
Film1 = [Film1_t, Al2O3_n, Al2O3_k, get_R(Al2O3_n, Al2O3_k)]

Film2_t = .02*10**-6
Film2 = [Film2_t, TiN_n, TiN_k, get_R(TiN_n, TiN_k)]

Domains.append(Air)
Domains.append(Film1)
Domains.append(Film2)


Amplitude = 1               # Constant amplitude of intitial signal. 
position_resolution = 101 


# This function calculates the transmission-only branch.
# It also produces a position array for the full domain.
def E_at_lambda(Domains, wavelengths, wavelength_index, Amplitude, resolution):
    lambda_vacuum = wavelengths[wavelength_index]


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
        n = domain[1][wavelength_index]
        k = domain[2][wavelength_index]

        lambda_n = get_lambda_n(f, n)
        x_min = 0
        x_max = domain[0]
        x = np.linspace(x_min,x_max, nx)

        hx = (x_max - x_min)/(nx-1)
        x_starts.append(x_max+x_starts[d]) # sets the start of the next values of x in the x_all list
        
        x_all = np.append(x_all, x+x_starts[d])
        E = np.zeros(nx)
        for  i in range(nx):
            E_vac = E_0*np.cos(2*np.pi/lambda_n*x[i]+2*np.pi*phi_starts[d])         # E(x) = E_0.cos(kx) = sqrt(I_0)
            E[i] = E_vac*np.sqrt(get_I_absorption(k, lambda_vacuum, x[i]))  # I(x) = I_0.e^(-alpha.x)
        
        E_all = np.append(E_all, E)
        E_0 = E_0*np.sqrt(get_I_absorption(k, lambda_vacuum, x[nx-1]))
        phi_starts.append(phi_starts[d]+x_max/lambda_n)
        
        d=d+1

    E_all = np.delete(E_all, 0)
    I_all = E_all**2
    x_all = np.delete(x_all, 0)
    return E_all, I_all, x_all

E_total, I_total, x_total = E_at_lambda(Domains, wavelengths, 0, Amplitude, position_resolution)

plt.plot(x_total,E_total)
# list[::-1] reverses the list

# indices of boundaries = (position_resolution)*d, d = domain in Domains

#Reflected_at_Film1(I_total,x_total):
print(I_total[position_resolution-1])
print(Domains[1][3][0])
reflected_intensity = I_total[position_resolution-1]*Domains[1][3][0]
print(reflected_intensity)

# E_total = E_total + add_reflected(reflected_intensity, positions, domain)
# def add_reflected():
def calc_E_through_layer(E_all, resolution, Amplitude, direction, phase, Domains, layer, wavelengths, wavelength_index):
    lambda_vacuum = wavelengths[wavelength_index]
    
    nx = resolution
    E_0 = Amplitude   
    E = np.zeros(nx)
    f = c/lambda_vacuum

    x_max = Domains[layer][0]
    n = Domains[layer][1][wavelength_index]
    k = Domains[layer][2][wavelength_index]
    
    lambda_n = get_lambda_n(f, n)

    x_min = 0
    x = np.linspace(x_min,x_max, nx)
    for  i in range(nx):
        E_vac = E_0*np.cos(2*np.pi/lambda_n*x[i]+2*np.pi*phase)         # E(x) = E_0.cos(kx) = sqrt(I_0)
        E[i] = E_vac*np.sqrt(get_I_absorption(k, lambda_vacuum, x[i]))  # I(x) = I_0.e^(-alpha.x)

    E_resulting = np.ndarray.copy(E_all)
    if direction == 1:
        E_resulting[layer*nx:(layer+1)*nx] = E_all[layer*nx:(layer+1)*nx] + E
    if direction == -1:
        E_reversed = E[::-1]
        E_resulting[layer*nx:(layer+1)*nx] = E_all[layer*nx:(layer+1)*nx] + E_reversed

    resulting_Amplitude = E_0*np.sqrt(get_I_absorption(k, lambda_vacuum, x[nx-1]))
    resulting_phase = phase + x_max/lambda_n

    return E_resulting, resulting_Amplitude, resulting_phase
   



"""
Need phases, which means reflections have to be done "in-line" as soon as a boundary is reached.
This means we should pre-allocate the position array:
"""
def make_x_axis(resolution):
    nx = resolution
    x_all = np.zeros(1)
    x_starts = [0]

    d=0
    for domain in Domains:
        x_min = 0
        x_max = domain[0]
        x = np.linspace(x_min,x_max, nx)

        hx = (x_max - x_min)/(nx-1)
        x_starts.append(x_max+x_starts[d]) # sets the start of the next values of x in the x_all list
        
        x_all = np.append(x_all, x+x_starts[d])
        d=d+1

    x_all = np.delete(x_all, 0)
    return x_all

x_domains = make_x_axis(position_resolution)
E_domains = np.zeros(len(x_domains))


propagation_direction = 1   # +/- 1 depending on the direction of the light propagation
propagation_phase = 0
current_layer = 0
sensitivity = 10**-6
E_domains2, resulting_Amplitude, resulting_phase = calc_E_through_layer(E_domains, position_resolution, Amplitude, propagation_direction,propagation_phase,Domains,current_layer,wavelengths,0)

"""
while current_layer + propagation_direction < len(Domains) and current_layer + propagation_direction > 0:   # This should be called within a branch
"""
def at_boundary(sensitivity, E_domains, position_resolution, boundary_Amplitude, propagation_direction, boundary_phase, Domains, current_layer, wavelengths,wavelength_index):
    current_layer = current_layer + propagation_direction
    # if (condition that allows transmission):
    if current_layer >= 0 and current_layer < len(Domains) and boundary_Amplitude > sensitivity:
    #   transmission
        print('T')
        E_domains, transmit_Amplitude, transmit_phase = calc_E_through_layer(E_domains, position_resolution, boundary_Amplitude, propagation_direction,boundary_phase,Domains,current_layer,wavelengths,wavelength_index)
        E_domains = at_boundary(sensitivity, E_domains, position_resolution, transmit_Amplitude, propagation_direction, transmit_phase, Domains, current_layer, wavelengths,wavelength_index)
    #   reflection
        print('R')
        reflected_Amplitude = boundary_Amplitude*Domains[current_layer][3][wavelength_index]  #first, set intensity of reflection when hitting the layer boundary
        propagation_direction = -propagation_direction                          #then, reverse direction
        current_layer = current_layer + propagation_direction                   #then, iterate direction backwards
        E_domains, reflect_Amplitude, reflect_phase = calc_E_through_layer(E_domains, position_resolution, reflected_Amplitude, propagation_direction,boundary_phase,Domains,current_layer,wavelengths,wavelength_index)
        E_domains = at_boundary(sensitivity, E_domains, position_resolution, reflect_Amplitude, propagation_direction, reflect_phase, Domains, current_layer, wavelengths,wavelength_index)
    return E_domains
# ^ this is not correct


E_domains_recursion = at_boundary(sensitivity, E_domains2, position_resolution, resulting_Amplitude, propagation_direction, resulting_phase, Domains, current_layer, wavelengths,0)
plt.plot(x_domains,E_domains_recursion)


plt.savefig('/workspaces/Projects/Physics/FiberModeling/outputEmissivity2.png')
plt.clf()

"""
basic structure is:
    1. initial portion in air
    2. test at_boundary
    3. On success, transmit and recur
        a. test at_boundary
        b. On success, transmit and recur
            i. test at_boundary
            ii. On failure, break
        c. After transmit, reflect and recur
            i. test at_boundary
            ii. On success, transmit and recur
                1. On failure, break
            iii. After transmit, reflect and recur
                1. On failure, break
    4. After transmit, reflect and recur
        a. On failure, break

At each at_boundary, carry forward the E calculated from the previous step.
The reflected branch can use the E calculated from the transmitted branch as the input E calculated
at_boundary returns the E_calc added from both transmitted and reflected paths
    Tests first,
    Transmits second,
    Reflects third,
    returns fourth

"""


"""
current_layer = current_layer + propagation_direction
#transmission1:
E_domains3, transmit_Amplitude, transmit_phase = calc_E_through_layer(E_domains2, position_resolution, resulting_Amplitude, propagation_direction,resulting_phase,Domains,current_layer,wavelengths,0)
#reflection1:
reflected_Amplitude = resulting_Amplitude*Domains[current_layer][3][0]  #first, set intensity of reflection when hitting the layer boundary
propagation_direction = -propagation_direction                          #then, reverse direction
current_layer = current_layer + propagation_direction                   #then, iterate direction backwards
#E_domains6, reflect_Amplitude, reflect_phase = calc_E_through_layer(E_domains5, position_resolution, reflected_Amplitude, propagation_direction,resulting_phase,Domains,current_layer,wavelengths,0)

#transmission2:
current_layer = current_layer + propagation_direction
E_domains4, dump, dump = calc_E_through_layer(E_domains3, position_resolution, transmit_Amplitude, propagation_direction,transmit_phase,Domains,current_layer,wavelengths,0)
#reflection2:
reflected_Amplitude = transmit_Amplitude*Domains[current_layer][3][0]  #first, set intensity of reflection when hitting the layer boundary
propagation_direction = -propagation_direction                          #then, reverse direction
current_layer = current_layer + propagation_direction
E_domains5, transmit_Amplitude, transmit_phase = calc_E_through_layer(E_domains4, position_resolution, reflected_Amplitude, propagation_direction,transmit_phase,Domains,current_layer,wavelengths,0)


#reflection1:
reflected_Amplitude = resulting_Amplitude*Domains[current_layer][3][0]  #first, set intensity of reflection when hitting the layer boundary
propagation_direction = -propagation_direction                          #then, reverse direction
current_layer = current_layer + propagation_direction                   #then, iterate direction backwards
E_domains6, reflect_Amplitude, reflect_phase = calc_E_through_layer(E_domains5, position_resolution, reflected_Amplitude, propagation_direction,resulting_phase,Domains,current_layer,wavelengths,0)

#plt.plot(x_domains,E_domains2)
#plt.plot(x_domains,E_domains3)
plt.plot(x_domains,E_domains4)
plt.plot(x_domains,E_domains5)
plt.plot(x_domains,E_domains6)
plt.legend(["air","transmitted1","transmitted2","reflected2","reflected1"])


plt.savefig('/workspaces/Projects/Physics/FiberModeling/outputEmissivity2.png')
plt.clf()"""

"""
Calculate intensity transmitted continuously through the material. 
Once a boundary is reached, calculate the amount reflected.
Call the "calculate intensity continuously through the material" function with the reflected intensity.
    Once the reflected portion is done calculating, move on the portion transmitted.

"""

