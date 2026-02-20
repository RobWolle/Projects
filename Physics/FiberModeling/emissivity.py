"""

Written by Robert Wolle for Raytum Photonics

"""
import time
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']=200 # setting the dpi of outputted graphs.


# Field strength E = E_0.sin(k.x - omega.t)

# k = 2pi/lambda, lamda changes in a material if f stays constant
# f = v/lambda
# v = c/n

E_0 = 1

c = 3*10**8

def get_lambda_n(f, n, c):
    v = c/n
    lambda_n = v/f
    return lambda_n

Air_n = 1.0003
Air_k = 0.05

Al2O3_1um_n = 1.6171
Al2O3_1um_k = 0

TiN_1um_n = 1.7547
TiN_1um_k = 3.4512


# Define regions of material by defining only:
# [0] = their widths <m> (to then call them in order to fill the space)
# [1] = their index of refraction n
# [2] = their extinction coefficient k

Domains = []

Air_L = 5*10**-6
Air = [Air_L, Air_n, Air_k]
Domains.append(Air)

Film1_L = 2.5*10**-6
Film1 = [Film1_L, Al2O3_1um_n, Al2O3_1um_k]
Domains.append(Film1)

Film2_L = 2.5*10**-6
Film2 = [Film2_L, TiN_1um_n, TiN_1um_k]
Domains.append(Film2)

print(Domains)


lambda_vacuum = 1*10**-6


f = c/lambda_vacuum



# Intensity loss:
def get_I_n(I_0, k_n, lambda_n, x_i):  # This is a different k (extinction coeff)
    alpha = 4*np.pi*k_n/lambda_n
    I_n = I_0*np.exp(-alpha*x_i)
    return I_n


nx = 101    # Resolution for each domain
x_all = np.zeros(1)
E_all = np.zeros(1)
x_starts = [0]
phi_starts = [0]

d=0
for domain in Domains:
    lambda_n = get_lambda_n(f, domain[1],c)
    x_min = 0
    x_max = domain[0]
    x = np.linspace(x_min,x_max, nx)

    x_starts.append(x_max+x_starts[d]) # sets the start of the next values of x in the x_all list
    #print(x+x_starts[d])
    x_all = np.append(x_all, x+x_starts[d])
    
    phi_starts.append(x_max/lambda_n)
    
    
    
    E = np.zeros(nx)
    for  i in range(nx):
        I_0 = E_0*np.cos(2*np.pi/lambda_n*x[i]-2*np.pi*phi_starts[d])         # E(x) = E_0.cos(kx) = I_0
        E[i] = get_I_n(I_0, domain[2], lambda_n, x[i])  # I(x) = I_0.e^(-alpha.x)
    E_all = np.append(E_all, E)
    E_0 = get_I_n(E_0, domain[2], lambda_n, x[nx-1])
    #print(E_0)

    d=d+1

    
plt.plot(x_all,E_all)

plt.savefig('/workspaces/Projects/Physics/FiberModeling/output.png')

