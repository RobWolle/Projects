
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

def printPlural(number):
    if number == 1:
        return
    else:
        return 's'

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

"""
Reading input n and k data files:
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

def LERP(array):
    # Creates an array where new values are added between the input array values via averages.
    new_length = len(array)*2-1
    new_array = np.zeros(new_length)
    i = 0

    while i < len(array)-1:
        new_array[i*2] = array[i]
        new_array[i*2+1] = (array[i]+array[i+1])/2
        i = i+1
    new_array[new_length-1] = array[len(array)-1]
    return new_array

# User inputs:
delimiter = ','
data_columns = 6
data_start = 0    #character position at which the data begins
file_path = "/workspaces/Projects/Physics/FiberModeling/n and k fits - TiN and Al2O3.csv"
DataArray = np.array(make_DataArray(delimiter, data_columns, data_start, file_path))

wavelengths_original = DataArray[:,0]*10**-6
wavelengths = LERP(wavelengths_original)
wavelength_resolution = len(wavelengths)



TiN_n_1 = DataArray[:,1]
TiN_n = LERP(TiN_n_1)
TiN_k_1 = DataArray[:,2]
TiN_k = LERP(TiN_k_1)

Al2O3_n_1 = DataArray[:,4]
Al2O3_n = LERP(Al2O3_n_1)
Al2O3_k_1 = DataArray[:,5]
Al2O3_k = LERP(Al2O3_k_1)


# Material n and k must be arrays which contain values at every wavelength in wavelengths[]

Air_n = np.ones(wavelength_resolution)*1.000273     # This is actually a measured thing that isn't exactly constant, but is almost constant in IR.
Air_k = np.ones(wavelength_resolution)*0.0005       # Technically this should be calculated through Kramers-Kronig, but it's close enough to zero.
Air_n = np.ones(wavelength_resolution)
Air_k = np.zeros(wavelength_resolution)

"""
Ideally, the data input would all be functionalized.
Users would input:
 1. the delimiter
 2. the number of data columns
 3. character start position
 4. filepath 

 5. wavelength column
 6. number of materials
 7. name of each material
 8. n and k columns of each material (in a list) of [column # of n(lambda), column # of k(lambda)]
    i. Optional Feature: "Calculate k from n" option


Applied to a Domain class which contains the material dictionary,
Or possibly just call Materials dictionary in the function

Then, the user would input:
 1. the signal amplitude of the lightsource (need to add wavelength dependence functionality)
 2. Position resolution (with reasonable bounds and time estimates given for the number of input wavelengths and layers)
 3. Sensitivity (with reasonable bounds as above)
"""

class Settings:
    def __init__(self, initial_Amplitude=1, position_resolution = 51, sensitivity=10**-4, Domains = 0):
        self.initial_Amplitude = initial_Amplitude
        self.position_resolution = position_resolution
        self.sensitivity = initial_Amplitude*sensitivity
        if Domains != 0:
            Domains.makeXaxis(position_resolution)

    def help(self):
        print("Settings has three user inputs:")
        print(" > initial_amplitude.. = The initial amplitude of the incoming light wave.")
        print("         Use 'Settings.initial_amplitude' to access or set this value.")
        print("         Default value is 1.")
        print(" > position_resolution = The number of points used to calculate the electric field between material boundaries.")
        print("         use 'Settings.position_resolution' to access or set this value.")
        print("         Default value is 51.")
        print(" > sensitivity........ = The proportion of initial wave amplitude at which the wave stops propagating once it hits a boundary between materials.")
        print("         use 'Settings.sensitivity' to access or set this value.")
        print("         Default value is 10^-4 (0.01*initial_amplitude).")

class Material:
    def __init__(material, name, n, k):
        material.name = name
        material.n = n
        material.k = k
        material.R = ((n-1)**2 + k**2)/((n+1)**2+k**2)

class Layer:
    def __init__(layer, thickness, material):
        layer.thickness = thickness
        layer.material = material

class Domains:
    def __init__(self):
        self.layer_count = 0
        self.domains = []
        self.blocks = [[0,0]]
        self.latest_block = 0
        self.measure_layer_count = 0
        self.measure_thickness = 0
        self.bulk_thickness = 0
        self.total_thickness = 0
        

    def info(self):
        print("This set of domains has " + str(self.layer_count) + " total layer" + printPlural(self.layer_count) + ":")
        if self.measure_layer_count > 0:
            print("     > " + str(self.measure_layer_count) + " measure layer" + printPlural(self.measure_layer_count) + " with thickness " + f"{self.measure_thickness*10**6:.2e}" + " um")
        if self.layer_count - self.measure_layer_count > 0:
            print("     > " + str(self.layer_count-self.measure_layer_count) + " bulk layer" + printPlural(self.measure_layer_count)+ " with total bulk thickness " + f"{(self.bulk_thickness)*10**6:.2e}" + " um")

    def AddLayer(self, thickness, Material):
        layer = Layer(thickness, Material)
        self.domains.append(layer)
        self.blocks[self.latest_block][1] = self.layer_count
        self.layer_count += 1
        self.total_thickness += thickness
        self.bulk_thickness += thickness

    def AddMeasureLayer(self, wavelengths, Material):
        thickness = 0.75*wavelengths[len(wavelengths)-1]
        layer = Layer(thickness, Material)
        self.domains.append(layer)
        self.blocks[self.latest_block][1] = self.layer_count
        self.layer_count += 1
        self.blocks[self.latest_block]=[self.layer_count,self.layer_count]
        self.measure_layer_count += 1
        self.measure_thickness = thickness
        self.total_thickness += thickness
    
    def RepeatAdd(self, number,input_block = [-1,-1]):
        # If no input_block is entered by the user, use the default:
        if input_block == [-1,-1]:
            input_block = [self.blocks[self.latest_block][0],self.blocks[self.latest_block][1]]
        
        block_start = input_block[0]
        block_end = input_block[1]+1

        block_size = block_end - block_start
        print("Repeating the following " + str(block_size) + " layers " + str(number) + " more times:")
        for j in range(block_start,block_end):
            self.printLayer(j)

        for i in range(number):
            for j in range(block_start,block_end):
                layer = self.domains[j]
                self.AddLayer(layer.thickness,layer.material)
    
        self.blocks.append([self.layer_count,self.layer_count])
        self.latest_block += 1

    def getLayer(self, layer_number):
        return self.domains[layer_number]
    
    def printLayer(self, layer_number):
        printed_layer = self.domains[layer_number]
        printed_material = printed_layer.material
        print("Layer "+str(layer_number) + " out of " + str(self.layer_count-1)+ ": Made of " + printed_material.name + " with thickness " + f"{printed_layer.thickness*10**6:.2e}" + " um")
    
    def printAllLayers(self):
        for i in range(self.layer_count):
            self.printLayer(i)

    def getLayer_property(self, layer_number, property, wavelength_index=0):
        desired_layer = self.domains[layer_number]
        desired_layer_material = desired_layer.material
        if property == 'thickness':
            return desired_layer.thickness
        elif property == 'n':
            desired_layer_n = desired_layer_material.n
            return desired_layer_n[wavelength_index]
        elif property == 'k':
            desired_layer_k = desired_layer_material.k
            return desired_layer_k[wavelength_index]
        elif property == 'R':
            desired_layer_R = desired_layer_material.R
            return desired_layer_R[wavelength_index]
        else:
            print("Domains.getLayer_property() failed. The only valid property names are 'thickness', 'n', 'k', and 'R'.")
    
    def makeXaxis(self, resolution):
        nx = resolution
        x_all = np.zeros(1)
        x_starts = [0]

        for d in range(self.layer_count-1):
            x_min = 0
            x_max = self.getLayer_property(d,'thickness')
            x = np.linspace(x_min,x_max, nx)
            x_starts.append(x_max+x_starts[d]) # sets the start of the next values of x in the x_all list
            
            x_all = np.append(x_all, x+x_starts[d])

        x_all = np.delete(x_all, 0)
        self.Xaxis = x_all
        return self.Xaxis


Air = Material('Air',Air_n,Air_k)
Al2O3 = Material('Al2O3',Al2O3_n,Al2O3_k)
TiN = Material('TiN',TiN_n,TiN_k)

Domains = Domains()
Domains.AddMeasureLayer(wavelengths, Air)
Domains.AddLayer(.8*10**-7, Al2O3)
Domains.AddLayer(.2*10**-7, TiN)
Domains.RepeatAdd(7)
Domains.AddLayer(.7*10**-7, Al2O3)
Domains.AddLayer(.3*10**-7, TiN)
Domains.RepeatAdd(4)
Domains.AddLayer(.5*10**-7, Al2O3)
Domains.AddLayer(.5*10**-7, TiN)
Domains.RepeatAdd(9)
Domains.AddMeasureLayer(wavelengths, Air)

print(Domains.blocks)
print(Domains.getLayer_property(0,'n',0))
print(Domains.layer_count)
Domains.printLayer(Domains.layer_count-2)
Domains.printLayer(Domains.layer_count-1)
Domains.info()






"""
Helpful Functions
"""

# k = 2pi/lambda, lamda changes in a material if f stays constant
# f = v/lambda
# v = c/n

# Intensity loss:
def get_I_absorption(k_n, lambda_vac, x_i):  # This is a different k (extinction coeff)
    alpha = 4*np.pi*k_n/lambda_vac
    I_absorption = np.exp(-alpha*x_i)
    return I_absorption


"""
SETTING UP MATERIAL DOMAINS
"""

Domains = []


# The air layers are meant to determine the reflected intensity vs transmitted intensity, 
#   so they should be thick enough to contain one peak of any input wavelength.
#   This can be calculated as half of the longest wavelength in the wavelenths[] array,
#   since this will certainly contain 1 peak.
last_wavelength_index = wavelength_resolution-1
Air_t = 7.5*wavelengths[last_wavelength_index]
Air = [Air_t, Air_n, Air_k, get_R(Air_n, Air_k)]


Film1_t = .05*10**-6
Film1 = [Film1_t, Al2O3_n, Al2O3_k, get_R(Al2O3_n, Al2O3_k)]

Film2_t = .05*10**-6
Film2 = [Film2_t, TiN_n, TiN_k, get_R(TiN_n, TiN_k)]

Film3_t = .08*10**-6
Film3 = [Film3_t, Al2O3_n, Al2O3_k, get_R(Al2O3_n, Al2O3_k)]

Film4_t = .02*10**-6
Film4 = [Film4_t, TiN_n, TiN_k, get_R(TiN_n, TiN_k)]

Film5_t = .07*10**-6
Film5 = [Film5_t, Al2O3_n, Al2O3_k, get_R(Al2O3_n, Al2O3_k)]

Film6_t = .03*10**-6
Film6 = [Film6_t, TiN_n, TiN_k, get_R(TiN_n, TiN_k)]

# Begin with an Air layer:
Domains.append(Air)    

# Add layers of alternating Film1 and Film2 pairs:
layers = 8             
for i in range(layers):

    Domains.append(Film3)
    Domains.append(Film4)

# Add layers of alternating Film3 and Film4 pairs:
layers = 5             
for i in range(layers):

    Domains.append(Film5)
    Domains.append(Film6)

layers = 10          
for i in range(layers):

    Domains.append(Film1)
    Domains.append(Film2)

# End with a layer of Air:
Domains.append(Air)         # This layer allows us to calculate the total transmission


Amplitude = 1               # Constant amplitude of intitial signal. 
position_resolution = 101 



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
   
def at_boundary(sensitivity, E_domains, position_resolution, boundary_Amplitude, propagation_direction, boundary_phase, Domains, current_layer, wavelengths,wavelength_index,branches):
    n_from = Domains[current_layer][3][wavelength_index]
    current_layer = current_layer + propagation_direction
    # if (condition that allows transmission):
    if current_layer >= 0 and current_layer < len(Domains) and boundary_Amplitude > sensitivity:
        n_to = Domains[current_layer][3][wavelength_index]
    #   transmission
        branches = branches + 1
        # As far as I can tell, T = (boundary intensity - R)
        reflected_Amplitude = boundary_Amplitude*Domains[current_layer][3][wavelength_index]  #first, set intensity of reflection when hitting the layer boundary
        transmitted_Amplitude = boundary_Amplitude - reflected_Amplitude
        E_domains, transmit_Amplitude, transmit_phase = calc_E_through_layer(E_domains, position_resolution, transmitted_Amplitude, propagation_direction,boundary_phase,Domains,current_layer,wavelengths,wavelength_index)
        E_domains,branches = at_boundary(sensitivity, E_domains, position_resolution, transmit_Amplitude, propagation_direction, transmit_phase, Domains, current_layer, wavelengths,wavelength_index,branches)
    #   reflection
        branches = branches + 1
        propagation_direction = -propagation_direction                          #then, reverse direction
        current_layer = current_layer + propagation_direction                   #then, iterate direction backwards
        if n_from < n_to:
            boundary_phase = boundary_phase + 0.5
        E_domains, reflect_Amplitude, reflect_phase = calc_E_through_layer(E_domains, position_resolution, reflected_Amplitude, propagation_direction,boundary_phase,Domains,current_layer,wavelengths,wavelength_index)
        E_domains,branches = at_boundary(sensitivity, E_domains, position_resolution, reflect_Amplitude, propagation_direction, reflect_phase, Domains, current_layer, wavelengths,wavelength_index,branches)
    return E_domains,branches

def get_I_T_R(I_domains, position_resolution, Domains):
    nx = position_resolution
    I_Reflected = max(I_domains[0:nx])
    I_Transmitted = max(I_domains[(len(Domains)-1)*nx:len(Domains)*nx])
    return I_Transmitted,I_Reflected

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
emissivity_lambda = np.ones(wavelength_resolution)

sensitivity = 10**-3

for k in range(wavelength_resolution):
    E_domains = np.zeros(len(x_domains))
    wavelength_index = k

    propagation_direction = 1       # +/- 1 depending on the direction of the light propagation
    propagation_phase = 0           # Initial phase of the incident wave where the air boundary begins
    current_layer = 0               # Initial layer of the incident wave

    # Sensitivity should depend on the input amplitude, and thus is actually a measured of % accuracy:
    Amplitude_sensitivity = Amplitude*sensitivity

    branches = 0                    # Testing variable that allows you to count the number of times E is calculated through a material.
                                    # Along with position_resolution, this determines the computational complexity of the program.


    # Starting with transmission through air alone actually allows me to get the phase and amplitude of the wave at the material without
    #   having to inlcude it in the final intensity.
    E_domains2, resulting_Amplitude, resulting_phase = calc_E_through_layer(E_domains, position_resolution, Amplitude, propagation_direction,propagation_phase,Domains,current_layer,wavelengths,wavelength_index)
    E_domains_recursion,branches = at_boundary(Amplitude_sensitivity, E_domains, position_resolution, resulting_Amplitude, propagation_direction, resulting_phase, Domains, current_layer, wavelengths,wavelength_index,branches)
    I_domains_recursion = E_domains_recursion**2
    T_total, R_total = get_I_T_R(I_domains_recursion,position_resolution, Domains)
    emissivity_lambda[k] = get_emissivity(T_total, R_total)



def get_emissivity_T(emissivity_lambda, wavelengths,T):
    dlambda = wavelengths[1]-wavelengths[0]
    integral1 = 0
    integral2 = 0
    for k in range(len(wavelengths)):
        Blackbody = get_BoltzmannLaw(wavelengths[k],T)
        integral1 = integral1 + Blackbody*emissivity_lambda[k]*dlambda
        integral2 = integral2 + Blackbody*dlambda
    emissivity_T = integral1/integral2
    return emissivity_T


Temperature = get_Kelvin(300) # in Celsius
emissivity_T = get_emissivity_T(emissivity_lambda,wavelengths,Temperature)
print(round(emissivity_T,5))

def axis_lim(array,min,max,second_array):
    # Inputs:
    #        array = axis being limited based on parameters
    #          min = value you want to be the minimum
    #          max = value you want to be the maximum
    # second_array = the other axis that needs to be trimmed

    index_length=len(array)

    true_max = array[index_length-1]
    true_min = array[0]
    max_ranged = max - true_min
    min_ranged = min - true_min

    max_proportion = max_ranged/(true_max-true_min)
    max_index = round(index_length*max_proportion)
    
    min_proportion = min_ranged/(true_max-true_min)
    min_index = round(index_length*min_proportion)

    output_array = array[min_index:max_index]
    output_second_array = second_array[min_index:max_index]
    return output_array,output_second_array

x,y = axis_lim(wavelengths,wavelengths[0],20*10**-6,emissivity_lambda)
plt.plot(x*10**6,y)
plt.ylim(0,1)
plt.savefig('/workspaces/Projects/Physics/FiberModeling/outputEmissivity2.png')
plt.clf()

"""plt.plot(wavelengths,Al2O3_k)
plt.savefig('/workspaces/Projects/Physics/FiberModeling/outputEmissivity2.png')
plt.clf()"""

