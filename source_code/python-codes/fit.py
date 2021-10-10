#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 17:05:19 2021

@author: nehabinish
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def read_file(filename):
    '''
    --------------------------------------------------------------------------
    function to read the file and store data in the file 
    -------------------------------------------------------------------------- 
    '''
    with open(filename,'r') as data:
        i = []
        j = []
        
        for line in data:
            c = line.rstrip( '\n' ) 
            c = line.split()
           
            i.append(float(c[0])) 
            j.append(float(c[1])) 
    return i,j


z_h,p_h = read_file('hollow_final.txt')
z_b1,p_b1 = read_file('b1_f.txt')
z_t,p_t = read_file('top_final.txt')




'''
 ----------------------------------------------------------------------------
 
                MORSE POTENTIAL MODEL FOR CURVE FITTING
 
 ----------------------------------------------------------------------------

'''

def morse(x, q, m, u):
    return (q * (np.exp(-2*m*(x-u))-2*np.exp(-m*(x-u))))

'''
 ----------------------------------------------------------------------------
 
             CURVE FITTING FOR HIGH SYMMETRY SITE - HOLLOW
 
 ----------------------------------------------------------------------------

'''

p0_h = [-6,1,0.1]

parameters_h, pcov_h = curve_fit(morse, z_h, p_h, p0_h, maxfev = 10000)
print('Parameters for hollow symmetry are {}'.format(parameters_h))

plt.figure(figsize=(20,10))

plt.plot(z_h,p_h,'ro-',label='data')
plt.plot(z_h, morse(z_h, *parameters_h),label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(parameters_h))

plt.grid()
plt.legend()
plt.xlabel('Z ($\AA$)')
plt.ylabel('Potential Energy (eV)')
plt.title('Fit for Hollow Symmetry')


'''
 ----------------------------------------------------------------------------
 
             CURVE FITTING FOR HIGH SYMMETRY SITE - BRIDGE
 
 ----------------------------------------------------------------------------

'''

p0_b1 = [-6,1,1]

parameters_b1, pcov_b1 = curve_fit(morse, z_b1, p_b1, p0_b1, maxfev = 10000)
print('Parameters for bridge symmetry are {}'.format(parameters_b1))

plt.figure(figsize=(20,10))

plt.plot(z_b1,p_b1,'ro-',label='data')
plt.plot(z_b1, morse(z_b1, *parameters_b1),label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(parameters_b1))

plt.grid()
plt.legend()
plt.xlabel('Z ($\AA$)')
plt.ylabel('Potential Energy (eV)')
plt.title('Fit for Bridge Symmetry')


'''
 ----------------------------------------------------------------------------
 
             CURVE FITTING FOR HIGH SYMMETRY SITE - TOP
 
 ----------------------------------------------------------------------------

'''


p0_t = [-8,2,3]

parameters_t, pcov_t = curve_fit(morse, z_t, p_t, p0_t, maxfev = 10000000)
print('Parameters for top symmetry are {}'.format(parameters_t))

plt.figure(figsize=(20,10))

plt.plot(z_t,p_t,'ro-',label='data')
plt.plot(z_t, morse(z_t, *parameters_t),label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(parameters_t))

plt.grid()
plt.legend()
plt.xlim(0,6)
plt.ylim(-8,1)
plt.xlabel('Z ($\AA$)') 
plt.ylabel('Potential Energy (eV)')
plt.title('Fit for Top Symmetry')


""
