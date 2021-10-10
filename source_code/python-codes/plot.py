#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 01:53:19 2021

@author: nehabinish
"""

import numpy as np
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D


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

def read_file_2(filename):
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


zcrp_h,pcrp_h = read_file_2('hollow_final.txt')
zcrp_b,pcrp_b = read_file_2('b1_f.txt')
zcrp_t,pcrp_t = read_file_2('top_final.txt')
zcrp_tbridge,pcrp_tbridge = read_file_2('top_bridge.txt')
zcrp_thollow,pcrp_thollow = read_file_2('top_hollow.txt')
zcrp_bridgehollow,pcrp_bridgehollow = read_file_2('bridge_hollow.txt')


zleps_h = np.array(read_file('z_arr_hollow.txt'))
Vleps_h = np.array(read_file('potential_hollow.txt'))


zleps_b = np.array(read_file('z_arr_bridgef.txt'))
Vleps_b = np.array(read_file('potential_bridgef.txt'))


zleps_t = np.array(read_file('z_arr_top.txt'))
Vleps_t = np.array(read_file('potential_top.txt'))


zleps_bridgehollow = np.array(read_file('z_arr_bridgehollow.txt'))
Vleps_bridgehollow = np.array(read_file('potential_bridgehollow.txt'))


zleps_tbridge = np.array(read_file('z_arr_tbridge.txt'))
Vleps_tbridge = np.array(read_file('potential_tbridge.txt'))


zleps_thollow = np.array(read_file('z_arr_thollow.txt'))
Vleps_thollow = np.array(read_file('potential_thollow.txt'))

'''
 ----------------------------------------------------------------------------
 
          1-D CUT - HIGH SYMMETRY SITE -  HOLLOW SYMMETRY
 
 ----------------------------------------------------------------------------

'''

plt.figure(figsize=(20,10))


plt.plot(zcrp_h,pcrp_h,'ro-',label= 'CRP Curve')


plt.plot(zleps_h,Vleps_h,'k-*', label = 'LEPS Curve')
plt.grid()
plt.xlabel('Z ($\AA$)')
plt.ylabel('Potential Energy (eV)')
plt.legend()
plt.title('1D cut HIGH SYMMETRY SITE - Hollow Symmetry')

'''
 ----------------------------------------------------------------------------
 
          1-D CUT - HIGH SYMMETRY SITE -  BRIDGE SYMMETRY
 
 ----------------------------------------------------------------------------

'''

plt.figure(figsize=(20,10))

plt.plot(zcrp_b,pcrp_b,'ro-',label= 'CRP Curve')

plt.plot(zleps_b,Vleps_b,'k-*', label = 'LEPS Curve')
plt.grid()
plt.xlabel('Z ($\AA$)')
plt.ylabel('Potential Energy (eV)')
plt.legend()
plt.xlim(0,6)
plt.ylim(-8,1)
plt.title('1D cut HIGH SYMMETRY SITE - Bridge Symmetry')


'''
 ----------------------------------------------------------------------------
 
          1-D CUT - HIGH SYMMETRY SITE -  TOP SYMMETRY
 
 ----------------------------------------------------------------------------

'''


plt.figure(figsize=(20,10))

plt.plot(zcrp_t,pcrp_t,'ro-',label= 'CRP Curve')

plt.plot(zleps_t,Vleps_t,'k-*', label = 'LEPS Curve')
plt.grid()
plt.xlabel('Z ($\AA$)')
plt.ylabel('Potential Energy (eV)')
plt.legend()
plt.xlim(0,6)
plt.ylim(-8,1)
plt.title('1D cut HIGH SYMMETRY SITE - Top Symmetry')

'''
 ----------------------------------------------------------------------------
 
          1-D CUT - LOW SYMMETRY SITE - TOP BRIDGE SYMMETRY
 
 ----------------------------------------------------------------------------

'''

plt.figure(figsize=(20,10))

plt.plot(zcrp_tbridge,pcrp_tbridge ,'ro-',label= 'CRP Curve')

plt.plot(zleps_tbridge,Vleps_tbridge,'k-*', label = 'LEPS Curve')
plt.grid()
plt.xlabel('Z ($\AA$)')
plt.ylabel('Potential Energy (eV)')
plt.legend()
plt.xlim(0,6)
plt.ylim(-8,1)
plt.title('1D cut LOW SYMMETRY SITE - Top Bridge Symmetry')


'''
 ----------------------------------------------------------------------------
 
          1-D CUT - LOW SYMMETRY SITE - TOP HOLLOW SYMMETRY
 
 ----------------------------------------------------------------------------

'''

plt.figure(figsize=(20,10))

plt.plot(zcrp_thollow,pcrp_thollow ,'ro-',label= 'CRP Curve')

plt.plot(zleps_thollow,Vleps_thollow,'k-*', label = 'LEPS Curve')
plt.grid()
plt.xlabel('Z ($\AA$)')
plt.ylabel('Potential Energy (eV)')
plt.legend()
plt.xlim(0,6)
plt.ylim(-8,1)
plt.title('1D cut LOW SYMMETRY SITE - Top Hollow Symmetry')

'''
 ----------------------------------------------------------------------------
 
          1-D CUT - LOW SYMMETRY SITE - BRIDGE HOLLOW SYMMETRY
 
 ----------------------------------------------------------------------------

'''

plt.figure(figsize=(20,10))

plt.plot(zcrp_bridgehollow,pcrp_bridgehollow ,'ro-',label= 'CRP Curve')

plt.plot(zleps_bridgehollow,Vleps_bridgehollow,'k-*', label = 'LEPS Curve')
plt.grid()
plt.xlabel('Z ($\AA$)')
plt.ylabel('Potential Energy (eV)')
plt.legend()
plt.xlim(0,6)
plt.ylim(-8,1)
plt.title('1D cut LOW SYMMETRY SITE - Bridge Hollow Symmebtry')

'''
 ----------------------------------------------------------------------------
 
                  COMPARISON OF HIGH SYMMETRY SITES
 
 ----------------------------------------------------------------------------

'''

plt.figure(figsize=(20,10))

plt.plot(zleps_t,Vleps_t,'k-*', label = 'TOP')
plt.plot(zleps_b,Vleps_b,'b-*', label = 'BRIDGE')
plt.plot(zleps_h,Vleps_h,'r-*', label = 'HOLLOW')

plt.grid()
plt.xlabel('Z ($\AA$)')
plt.ylabel('Potential Energy (eV)')
plt.legend()
plt.xlim(0,6)
plt.ylim(-8,1)
plt.title('Comparison of the High Symmetry Sites')



'''
 ----------------------------------------------------------------------------
 
                  COMPARISON OF LOW SYMMETRY SITES
 
 ----------------------------------------------------------------------------

'''


plt.figure(figsize=(20,10))

plt.plot(zleps_tbridge,Vleps_tbridge,'k-*', label = 'TOP-BRIDGE')
plt.plot(zleps_bridgehollow,Vleps_bridgehollow,'b-*', label = 'BRIDGE-HOLLOW')
plt.plot(zleps_thollow,Vleps_thollow,'r-*', label = 'TOP-HOLLOW')

plt.grid()
plt.xlabel('Z ($\AA$)')
plt.ylabel('Potential Energy (eV)')
plt.legend()
plt.xlim(0,6)
plt.ylim(-8,30)
plt.title('Comparison of the Low Symmetry Sites')


'''
 ----------------------------------------------------------------------------
 
               COMPARISON OF ALL THE SYMMETRY SITES
 
 ----------------------------------------------------------------------------

'''

plt.figure(figsize=(20,10))

plt.plot(zleps_tbridge,Vleps_tbridge,'-*', label = 'TOP-BRIDGE')
plt.plot(zleps_bridgehollow,Vleps_bridgehollow,'-*', label = 'BRIDGE-HOLLOW')
plt.plot(zleps_thollow,Vleps_thollow,'-*', label = 'TOP-HOLLOW')
plt.plot(zleps_t,Vleps_t,'-*', label = 'TOP')
plt.plot(zleps_b,Vleps_b,'-*', label = 'BRIDGE')
plt.plot(zleps_h,Vleps_h,'-*', label = 'HOLLOW')

plt.grid()
plt.xlabel('Z ($\AA$)')
plt.ylabel('Potential Energy (eV)')
plt.legend()
plt.xlim(0,6)
plt.ylim(-8,20)
plt.title('Comparison of all Symmetry Sites')



