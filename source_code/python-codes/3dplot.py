#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 23:26:19 2021

@author: nehabinish
"""

import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

def read_file(filename):
    '''
    --------------------------------------------------------------------------
    function to read the file and store data in the file 
    -------------------------------------------------------------------------- 
    '''
    with open(filename,'r') as data:
        
        arr = []
        
        for line in data:
            arr.append(float(line))

    return arr

potential = np.zeros((200,200))

with open('potential3d.txt','r') as data:
        
    arr = []
    
    i = 0 
    count = 0

    for line in data:
        
        if line.split() == []:
            i += 1
            count = 0

        else:
            potential[i,count] = float(line)    
            count += 1    
            

z = np.array(read_file('z_arr3d.txt'))
y = np.array(read_file('y_arr3d.txt'))

Y, Z = np.meshgrid(y, z)

fig = plt.figure(figsize=(10,5))
ax = plt.axes(projection='3d')
surf = ax.contour3D(Y, Z, potential, 50, cmap='hot')
fig.colorbar(surf)
ax.view_init(10,250)
ax.set_xlabel('Y ($\AA$)')
ax.set_ylabel('Z ($\AA$)')
ax.set_zlabel('potential (eV)')
plt.title('Plot generated from the 3D-PES script, X = 0.5')
plt.xlim(0,1.4)


