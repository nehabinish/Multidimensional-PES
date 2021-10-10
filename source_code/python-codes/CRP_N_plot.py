#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 22:29:35 2021

@author: nehabinish
"""
import numpy as np
import matplotlib.pyplot as plt

def read_file(filename):
    '''
    --------------------------------------------------------------------------
    function to read the file and store data in the file into two lists 
    giving the linked nodes.
    -------------------------------------------------------------------------- 
    '''
    with open(filename,'r') as data:
        i = []
        j = []
        
        for line in data:
            c = line.rstrip( '\n' ) 
            c = line.split()
            # creating lists for the linked lists
            i.append(float(c[0])) # out node
            j.append(float(c[1])) # in node

    return i,j


plt.figure(figsize=(20,10))

z_top,p_top = read_file('top.txt')

plt.plot(z_top,p_top,label='Top')

z_h,p_h = read_file('hollow.txt')
plt.plot(z_h,p_h,'ro-',label='Hollow')

z_b1,p_b1 = read_file('b1.txt')
plt.plot(z_b1,p_b1,'k-*',label='Bridge 1')

z_b2,p_b2 = read_file('b2.txt')
plt.plot(z_b2,p_b2,'-',label='Bridge 2')




plt.ylim(-8,1)
plt.grid()
plt.legend()
plt.xlabel('Z ($\AA$)')
plt.ylabel('Potential Energy (eV)')
plt.title('1-D cut N/W(100) system')

for i in p_top:
    print(i)
