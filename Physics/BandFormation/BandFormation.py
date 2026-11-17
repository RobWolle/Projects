# -*- coding: utf-8 -*-
"""
This code is intended to help calculate the energy of bonding and antibonding orbitals,
    specifically in solids.
    
Created by Robert Wolle, 10/7/2024
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']=200


# Angstrom rescaling factor, puts everything in the units of Angstroms
Angstrom = 10**-10

# Bohr radius:
a0 = 5.29177210903*10**-11/Angstrom

# Number of protons of interest (Copper -> 29)
Z = 29

# Lattice constant in solid copper
a = 2.56*10**-10/Angstrom

def radial_3d(r,Z,a0):
    # This calculates the radial distribution function:
    R = 8/(15*81**2)*Z**7*a0**6*r**6*np.exp(-2*Z*a0*r/3)
    return R

def wavefunction_3d(r,Z,a0):
    psi = 4/(81*np.sqrt(30))*Z**(7/2)*a0**2*r**2*np.exp(-Z*a0*r/3)
    return psi

def radial_4s(r, Z, a0):
    R = 4*np.pi/96**2*Z**3*a0**2*r**2*(24-18*Z*a0*r + 3*Z**2*a0**2*r**2 - Z**3*a0**3*r**3/8)**2*np.exp(-Z*a0*r/2)
    return R

def wavefunction_4s(r, Z, a0):
    psi = 1/96*Z**(3/2)*(24-18*Z*a0*r + 3*Z**2*a0**2*r**2 - Z**3*a0**3*r**3/8)*np.exp(-Z*a0*r/4)
    return psi

def get_i(target,range_settings):
    r_min = range_settings[1]
    dr = range_settings[3]
    i = round((target-r_min)/dr)
    return i
    
def integral(interval, function, range_settings):
    integral=0
    dr = range_settings[3]
    for i in range(get_i(interval[0],range_settings),get_i(interval[1],range_settings)):
        integral = integral + function[i]*dr
    return integral

"""
#####  USER INPUTS  ######
"""
# Position of each ion site in Angstroms:
proton_positions   =  [-0.5*a,0.5*a]
# The sign of each wavefunction (+/-1):
wavefunction_signs =  [1,-1,1,-1,-1,1]   # The number of entries doesn't have to equal the number of protons


# Region of interest:
r_min = -3*a  # in Angstroms
r_max = 3*a   # in Angstroms

# Number of points to calculate:
resolution = 1000   # 1000 is good enough and fast

dr = (r_max - r_min)/(resolution-1)
range_settings = [r_max,r_min,resolution,dr]

"""
#####  CREATING SPACE  #####
"""

r = np.linspace(r_min,r_max,resolution)
mod_rR3d = np.zeros(resolution)
R3d = np.zeros(resolution)
psi3d = np.zeros(resolution)
R4s = np.zeros(resolution)
mod_rR4s = np.zeros(resolution)
psi4s = np.zeros(resolution)



####################
"""
# This code calculates the radial distribution function directly
for r_p in proton_positions:
    for i in range(resolution):
        dist = abs(r[i]-r_p)
        R3d[i] = R3d[i] + radial_3d(dist, Z,a0)
        R4s[i] = R4s[i] + radial_4s(dist,Z,a0)
plt.plot(r,R3d)
plt.plot(r,R4s)
plt.show()"""
###################

"""
#####  CALCULATING WAVEFUNCTIONS  #####
 1. Calculate radial wavefunction at each position given each proton's position and sign.
 2. Calculate pure square root of radial probability distribution by multiplying the wavefunction at each position by wavefunction sign and distance to protons
 3. Square the radial probability distributions
 4. Plot the radial distribtutions and wavefunctions
"""

s=0
for r_p in proton_positions:
    for i in range(resolution):
        dist = abs(r[i]-r_p)
        psi3d[i] = psi3d[i] + wavefunction_signs[s]*wavefunction_3d(dist, Z,a0)
        psi4s[i] = psi4s[i] + wavefunction_signs[s]*wavefunction_4s(dist,Z,a0)
    s+=1
    
s=0
for r_p in proton_positions:
    for i in range(resolution):
        dist = abs(r[i]-r_p)
        mod_rR3d[i] = mod_rR3d[i] + wavefunction_signs[s]*wavefunction_3d(dist, Z,a0)*dist*a0
        mod_rR4s[i] = mod_rR4s[i] + wavefunction_signs[s]*wavefunction_4s(dist,Z,a0)*dist*a0
    s+=1
    
s=0
for r_p in proton_positions:
    for i in range(resolution):
        dist = abs(r[i]-r_p)
        R3d[i] = R3d[i] + wavefunction_signs[s]*wavefunction_3d(dist, Z,a0)
        R4s[i] = R4s[i] + wavefunction_signs[s]*wavefunction_4s(dist,Z,a0)
    s+=1

mod_R_distribution_3d = mod_rR3d**2
mod_R_distribution_4s = mod_rR4s**2*4*np.pi







"""
#######  CALCULATING AVERAGE VALUE OF RADIAL DISTRIBUTION
 1. Take integral over defined interval, divide by size of the integral.
"""
interval_full = [r_min, r_max]

interval1 = [a, 2*a]

integral_mod_R4s = integral(interval1,mod_rR4s*psi4s, range_settings)

integral_R3d = integral(interval_full, psi3d**2, range_settings)

integral_R4s = integral(interval_full, psi4s**2, range_settings)


avgR4s = integral_mod_R4s/(interval1[1]-interval1[0])
print(avgR4s)


plt.plot(r,psi3d)
plt.plot(r,psi4s)
plt.show()
plt.savefig('/workspaces/Projects/Physics/BandFormation/output1.png')
plt.clf()
plt.plot(r,mod_R_distribution_3d/integral_R3d)
plt.plot(r,mod_R_distribution_4s/integral_R4s)
plt.show()
plt.savefig('/workspaces/Projects/Physics/BandFormation/output2.png')

