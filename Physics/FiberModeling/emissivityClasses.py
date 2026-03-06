
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
Printing and utility functions
"""

class Timer:
    def __init__(self):
        self.initial_time = time.time()
        self.event_times = [self.initial_time]
        self.event_names = ['start']
        self.current_event = 0
        self.decimals = '5'

    def timeTaken(self, event_name='', do_print = True,time_since='last', decimals=5):
        self.decimals = decimals
        self.current_event +=1
        self.event_times.append(time.time())
        self.event_names.append(event_name)
        current_time = self.event_times[self.current_event]
        
        prefix = event_name
        if do_print == True:
            if  prefix == '':
                prefix = "Event " + str(self.current_event)
            if time_since == 'start':
                last_time = self.initial_time
            else:
                last_time = self.event_times[self.current_event-1]
            print(prefix + " took " + str(round(current_time-last_time,self.decimals)) + ' s')
    
    def printTimes(self, time_since='last', decimals=5):
        self.decimals = decimals
        print("\nThere were " + str(self.current_event) + " time"+ printPlural(self.current_event) + " taken:")
        for t in range(1,self.current_event+1):
            current_time = self.event_times[t]
            prefix = self.event_names[t]
            if  prefix == '':
                prefix = "Event " + str(t)
            if time_since == 'start':
                last_time = self.initial_time
            else:
                last_time = self.event_times[t-1]
            print(prefix + " took " + str(round(current_time-last_time,self.decimals)) + ' s')
        print('')

    def help(self):
        print("To use the timer, use the timeTaken() method.")
        print("timeTaken() has a few options for printing the time ellapsed:")
        print(" > do_print is a boolean that enables or disables printing the time.")
        print("     do_print is True by default.")
        print(" > time_since is a setting that determines what time is printed by the method.")
        print("     time_since = 'last' (default setting) means the printed time is the time elapsed since the last timeTaken() call.")
        print("     time_since = 'start' means the printed time is the time elapsed since the Timer object was created.")
        print(" > event_name is a name given to the current time event. This is an optional setting that determines how the events are printed.")
        print(" > decimals determines the number of significant figures displayed.")
        print("     decimals = 5 by default.")

def printPlural(number):
    if number == 1:
        return ''
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

Timer = Timer()

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

wavelength_resolution = len(wavelengths)






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

    def plot_material_property(material, property, wavelengths):
        if property == 'n':
            y = material.n
        elif property == 'k':
            y = material.k
        elif property == 'R':
            y = material.R
        else:
            print("Material.plot_material_property() failed. The only valid property names are 'n', 'k', and 'R'.")
            return
        plt.plot(wavelengths,y)
        plt.savefig('/workspaces/Projects/Physics/FiberModeling/output' + material.name + property + '.png')
        plt.clf()

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
        self.Xaxis = 'The Xaxis has not been made yet. Create a Settings object before retrieving this array.'

    def info(self):
        print("This set of domains has " + str(self.layer_count) + " total layer" + printPlural(self.layer_count) + ":")
        if self.measure_layer_count > 0:
            print("     > " + str(self.measure_layer_count) + " measure layer" + printPlural(self.measure_layer_count) + " with thickness " + f"{self.measure_thickness*10**6:.2e}" + " um")
        if self.layer_count - self.measure_layer_count > 0:
            print("     > " + str(self.layer_count-self.measure_layer_count) + " bulk layer" + printPlural(self.layer_count-self.measure_layer_count)+ " with total bulk thickness " + f"{(self.bulk_thickness)*10**6:.2e}" + " um")

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
        print("Repeating the following " + str(block_size) + " layer" + printPlural(block_size) + " " + str(number) + " more time" + printPlural(number) + ":")
        for j in range(block_start,block_end):
            self.printLayer(j)
        print('')
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

Domains.info()

Timer.timeTaken('Domains')




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
SETTINGS

It seems like wavelength_resolution should be contained in Settings,
    but it's kind of awkward to do that.
    Without that, we need three total things:
        1. Domains
        2. Settings
        3. wavelength_resolution
"""


Amplitude = 1               # Constant amplitude of intitial signal. 
position_resolution = 51 
sensitivity = 10**-4
Settings = Settings(Amplitude, position_resolution, sensitivity, Domains)


def calc_E_through_layer(E_all, Amplitude, direction, phase, Domains, layer, Settings, wavelengths, wavelength_index):
    lambda_vacuum = wavelengths[wavelength_index]
    
    nx = Settings.position_resolution
    E_0 = Amplitude   
    f = c/lambda_vacuum
    
    x_max = Domains.getLayer_property(layer,'thickness')
    n = Domains.getLayer_property(layer,'n',wavelength_index)
    k = Domains.getLayer_property(layer,'k',wavelength_index)
    
    lambda_n = get_lambda_n(f, n)
    x_min = 0
    x = np.linspace(x_min,x_max, nx)
    
    E_vac = E_0*np.cos(2*np.pi/lambda_n*x+2*np.pi*phase)         # E(x) = E_0.cos(kx) = sqrt(I_0)
    E = E_vac*np.sqrt(get_I_absorption(k, lambda_vacuum, x))  # I(x) = I_0.e^(-alpha.x)
    
    E_resulting = np.ndarray.copy(E_all)

    if direction == 1:
        E_resulting[layer*nx:(layer+1)*nx] = E_all[layer*nx:(layer+1)*nx] + E
    elif direction == -1:
        E_reversed = E[::-1]
        E_resulting[layer*nx:(layer+1)*nx] = E_all[layer*nx:(layer+1)*nx] + E_reversed

    resulting_Amplitude = E_0*np.sqrt(get_I_absorption(k, lambda_vacuum, x[nx-1]))
    resulting_phase = phase + x_max/lambda_n

    return E_resulting, resulting_Amplitude, resulting_phase


def at_boundary(E_domains, boundary_Amplitude, propagation_direction, boundary_phase, Settings, Domains, current_layer, wavelengths,wavelength_index,branches):
    n_from = Domains.getLayer_property(current_layer, 'n',wavelength_index)
    current_layer += propagation_direction
    # if (condition that allows transmission):
    if current_layer >= 0 and current_layer < Domains.layer_count-1 and boundary_Amplitude > Settings.sensitivity:
        n_to = Domains.getLayer_property(current_layer, 'n',wavelength_index)
    #   transmission
        branches += 1
        # As far as I can tell, T = (boundary intensity - R)
        reflected_Amplitude = boundary_Amplitude*Domains.getLayer_property(current_layer, 'R',wavelength_index)  #first, set intensity of reflection when hitting the layer boundary
        transmitted_Amplitude = boundary_Amplitude - reflected_Amplitude
        #                                                                   E_all,          Amplitude,          direction,          phase,      Domains,    layer,    Settings, wavelengths, wavelength_index
        E_domains, transmit_Amplitude, transmit_phase = calc_E_through_layer(E_domains, transmitted_Amplitude, propagation_direction,boundary_phase,Domains,current_layer,Settings, wavelengths,wavelength_index)
        #                               (E_domains, boundary_Amplitude, propagation_direction, boundary_phase, Settings, Domains, current_layer, wavelengths,wavelength_index,branches)
        E_domains,branches = at_boundary(E_domains, transmit_Amplitude, propagation_direction, transmit_phase, Settings, Domains, current_layer, wavelengths,wavelength_index,branches)
    #   reflection
        branches += 1
        propagation_direction = -propagation_direction                          #then, reverse direction
        current_layer += propagation_direction                   #then, iterate direction backwards
        if n_from < n_to:
            boundary_phase += 0.5
        #                                                                   E_all,          Amplitude,          direction,          phase,      Domains,    layer,    Settings, wavelengths, wavelength_index
        E_domains, reflect_Amplitude, reflect_phase = calc_E_through_layer(E_domains, reflected_Amplitude, propagation_direction,boundary_phase,Domains,current_layer,Settings, wavelengths, wavelength_index)
        #                               (E_domains, boundary_Amplitude, propagation_direction, boundary_phase, Settings, Domains, current_layer, wavelengths,wavelength_index,branches)
        E_domains,branches = at_boundary(E_domains, reflect_Amplitude, propagation_direction, reflect_phase, Settings, Domains, current_layer, wavelengths,wavelength_index,branches)
    return E_domains,branches

def get_I_T_R(I_domains, Settings, Domains):
    nx = Settings.position_resolution
    I_Reflected = max(I_domains[0:nx])
    I_Transmitted = max(I_domains[(Domains.layer_count-2)*nx:(Domains.layer_count-1)*nx])
    return I_Transmitted,I_Reflected

def get_emissivity_lambda(Domains, Settings, wavelengths, wavelength_resolution):
    emissivity_lambda = np.ones(wavelength_resolution)

    x_domains = Domains.Xaxis
    E_domains = np.zeros(len(x_domains))
    for k in range(wavelength_resolution):
        propagation_direction = 1       # +/- 1 depending on the direction of the light propagation
        propagation_phase = 0           # Initial phase of the incident wave where the air boundary begins
        current_layer = 0               # Initial layer of the incident wave

        branches = 0
        #                                                                           E_all,              Amplitude,              direction,          phase,      Domains, layer,      Settings, wavelengths, wavelength_index
        E_domains2, resulting_Amplitude, resulting_phase = calc_E_through_layer(E_domains, Settings.initial_Amplitude, propagation_direction,propagation_phase,Domains,current_layer,Settings, wavelengths,k)
        #                                           (E_domains, boundary_Amplitude, propagation_direction, boundary_phase, Settings, Domains, current_layer, wavelengths,wavelength_index,branches)
        E_domains_recursion,branches = at_boundary(E_domains, resulting_Amplitude, propagation_direction, resulting_phase, Settings, Domains, current_layer, wavelengths,k,branches)
        I_domains_recursion = E_domains_recursion**2
        print("Done " + str(k) + "/" + str(wavelength_resolution) + " wavelength"+ printPlural(wavelength_resolution)+ " with " +str(branches)+" branche"+ printPlural(branches)+".")
        T_total, R_total = get_I_T_R(I_domains_recursion, Settings, Domains)
        emissivity_lambda[k] = get_emissivity(T_total, R_total)
    return emissivity_lambda

def get_emissivity_T(emissivity_lambda, wavelengths,T):
    dlambda = wavelengths[1]-wavelengths[0]
    integral1 = 0
    integral2 = 0
    for k in range(len(wavelengths)):
        Blackbody = get_BoltzmannLaw(wavelengths[k],T)
        integral1 = integral1 + Blackbody*emissivity_lambda[k]*dlambda
        integral2 = integral2 + Blackbody*dlambda
    emissivity_T = integral1/integral2

    print("\nemissivity at " + str(T-get_Kelvin(0)) + " C is " + str(round(emissivity_T,5)))

    return emissivity_T


emissivity_lambda = get_emissivity_lambda(Domains, Settings, wavelengths, wavelength_resolution)

Timer.timeTaken('emissivity_lambda calculation')

Temperature = get_Kelvin(150) # in Celsius
emissivity_T = get_emissivity_T(emissivity_lambda,wavelengths,Temperature)

Temperature = get_Kelvin(300) # in Celsius
emissivity_T = get_emissivity_T(emissivity_lambda,wavelengths,Temperature)

Temperature = get_Kelvin(450) # in Celsius
emissivity_T = get_emissivity_T(emissivity_lambda,wavelengths,Temperature)

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

Al2O3.plot_material_property('k',wavelengths)

Timer.printTimes()

