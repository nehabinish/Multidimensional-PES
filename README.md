# Multidimensional-PES
The goal of this project is to build multidimensional potential energy surface to describe N-N2 interacting with the tungsten surfaces, specifically for the crystallographic plane W(100). The potential energy surfaces are build using an analytical model called the LEPS (London Eyring Polanyi Sato) model. Initially, the required potential energies are extracted from the numerical model called the CRP model that has been provided for the N2/W(100) system. These potential energies are then used to carry out a curve fitting to obtain the various parameters of the LEPS model. The CRP-PES also behaves as a reference model that serves to compare the quality of the LEPS- PES generated. The attained version of the program builds a 3D-PES for the N/W(100) atom- surface

In order to increase readability, the project is made of several files. The curve fitting and the plots are made using python3, whereas the development of the LEPS potential energy surfaces is done using FORTRAN90.
fit.py : Non - linear curve fitting in Python3
main.f90 (main program) : Building of the LEPS potentials sub.f90: Modules required to generate 3D-PES scripts and 1-D cuts savearray.py : Python Program to read text files and save as arrays. 3dplot.py: Plot of the 3-D-PES script
error.f90 : Program for Error Analysis
min.f90 : Program to find the Global Minimum
