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

# k = 2pi/lambda, lamda changes in a material if f stays constant
# f = v/lambda
# v = c/n

def get_lambda_n(f, n, c):
    v = c/n
    lambda_n = v/f
    return lambda_n

"""
Material Constants
"""
Air_n = 1.00027417
Air_k = 0.005

Al2O3_1um_n = 1.6171
Al2O3_1um_k = 0

TiN_1um_n = 1.7547
TiN_1um_k = 3.4512
#TiN_1um_k = 0.4512

"""
Input Conditions
"""
Amplitude = 1               # Constant amplitude of intitial signal. 
lambda_vacuum = 1*10**-6    # Desired vacuum wavelength of light.
nx = 1001                    # Resolution for each domain.

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



# Intensity loss:
def get_I_absorption(k_n, lambda_vac, x_i):  # This is a different k (extinction coeff)
    alpha = 4*np.pi*k_n/lambda_vac
    I_absorption = np.exp(-alpha*x_i)
    return I_absorption

E_0 = Amplitude
f = c/lambda_vacuum

x_all = np.zeros(1)
E_all = np.zeros(1)
I_all = np.zeros(1)
x_starts = [0]
phi_starts = [0]
d=0
for domain in Domains:
    lambda_n = get_lambda_n(f, domain[1],c)
    x_min = 0
    x_max = domain[0]
    x = np.linspace(x_min,x_max, nx)

    hx = (x_max - x_min)/(nx-1)
    x_starts.append(x_max+x_starts[d]) # sets the start of the next values of x in the x_all list
    
    x_all = np.append(x_all, x+x_starts[d])
    
    E = np.zeros(nx)
    I = np.zeros(nx)
    for  i in range(nx):
        E_vac = E_0*np.cos(2*np.pi/lambda_n*x[i]+2*np.pi*phi_starts[d])         # E(x) = E_0.cos(kx) = sqrt(I_0)
        E[i] = E_vac*np.sqrt(get_I_absorption(domain[2], lambda_vacuum, x[i]))  # I(x) = I_0.e^(-alpha.x)
        I[i] = E[i]**2
    
    
    E_all = np.append(E_all, E)
    I_all = np.append(I_all, I)
    E_0 = E_0*np.sqrt(get_I_absorption(domain[2], lambda_vacuum, x[nx-1]))
    phi_starts.append(phi_starts[d]+x_max/lambda_n)
    
    print(x_starts[d])
    print(E)
    d=d+1

E_all = np.delete(E_all, 0)
I_all = np.delete(I_all, 0)
x_all = np.delete(x_all, 0)


plt.plot(x_all,E_all)

plt.savefig('/workspaces/Projects/Physics/FiberModeling/outputE.png')

plt.clf()
plt.plot(x_all,I_all)

plt.savefig('/workspaces/Projects/Physics/FiberModeling/outputI.png')

