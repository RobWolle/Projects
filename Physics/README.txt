This folder contains scripts for various simulations of physical systems.

Light-and-Materials
  Contains programs which simulate interactions between light and user-input materials. These programs can calculate emissivity, 
  help design anti-reflective coatings, and simulate internal reflections in highly complex coating designs.

FiberModeling
  Contains programs which model the modes and features of fiber optic waveguides depending on input parameters.

PeakFitting
  Contains programs which load input data with arbritary delimitters, column, and row counts, plot that data, and plot 
  user-input peaks (either Lorentzian or Gaussian) + background noise found to match the data. Currently working on a program 
  to automatically fit peaks to data with gradient descent. The program was built for Raman spectroscopy and Photoluminescence 
  data fitting.

BandFormation.py
  Calculate the energy of bonding and antibonding orbitals in solids. This currently only works in 1D to calculate the strength
  of overlapping orbitals to visualize the formation of energy bands from the splitting of energy states.

Electrostatic Simulation.py & Electrostatic Modeling Presentation.pdf
  The Python script uses finite difference methods to calculate the electric potential and field in 2D around user-input square
  or triangular geometry held at input voltages. The presentation explains its function and results in the case of a particular 
  problem involving a conical ion trap.

TMDTunnelingProbability.py
  Calculate the probability of a charged particle tunneling through a user-input lattice of 2D material based on the 
  WKB approximation. Supports TMD monolayers, and any number of layers or switching between them.

Tight Binding Contours Midterm prob3.py
  Calculates and plots contours of equal energy in k-space for a square lattice material in the tight binding approximation.

Tight Binding Energy Bands.py
  Calculates and plots the energy bands of graphene's Brillouin Zone
