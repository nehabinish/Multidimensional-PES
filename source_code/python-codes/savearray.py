#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 23:08:59 2021

@author: nehabinish
"""

import numpy as np

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

            if (float(c[0])<7.0):
                i.append(float(c[0])) 
                j.append(float(c[1])) 

    return i,j

zcrp_h,pcrp_h = read_file_2('hollow_final.txt')
zcrp_b,pcrp_b = read_file_2('b1_f.txt')
zcrp_t,pcrp_t = read_file_2('top_final.txt')
zcrp_tbridge,pcrp_tbridge = read_file_2('top_bridge.txt')
zcrp_thollow,pcrp_thollow = read_file_2('top_hollow.txt')
zcrp_bridgehollow,pcrp_bridgehollow = read_file_2('bridge_hollow.txt')


with open('zcrp_h.txt', 'w') as filehandle: 
    for listitem in zcrp_h:
        filehandle.write('%s\n' % listitem)

with open('zcrp_b.txt', 'w') as filehandle:
    for listitem in zcrp_b:
        filehandle.write('%s\n' % listitem)
  
with open('zcrp_t.txt', 'w') as filehandle:
    for listitem in zcrp_t:
        filehandle.write('%s\n' % listitem)   

with open('zcrp_tbridge.txt', 'w') as filehandle:
    for listitem in zcrp_tbridge:
        filehandle.write('%s\n' % listitem)        
        
with open('zcrp_thollow.txt', 'w') as filehandle:
    for listitem in zcrp_thollow:
        filehandle.write('%s\n' % listitem)        

with open('zcrp_bridgehollow.txt', 'w') as filehandle:
    for listitem in zcrp_bridgehollow:
        filehandle.write('%s\n' % listitem)
        

with open('pcrp_h.txt', 'w') as filehandle:
    for listitem in pcrp_h:
        filehandle.write('%s\n' % listitem)

with open('pcrp_b.txt', 'w') as filehandle:
    for listitem in pcrp_b:
        filehandle.write('%s\n' % listitem)
  
with open('pcrp_t.txt', 'w') as filehandle:
    for listitem in pcrp_t:
        filehandle.write('%s\n' % listitem)   

with open('pcrp_tbridge.txt', 'w') as filehandle:
    for listitem in pcrp_tbridge:
        filehandle.write('%s\n' % listitem)        
        
with open('pcrp_thollow.txt', 'w') as filehandle:
    for listitem in pcrp_thollow:
        filehandle.write('%s\n' % listitem)        

with open('pcrp_bridgehollow.txt', 'w') as filehandle:
    for listitem in pcrp_bridgehollow:
        filehandle.write('%s\n' % listitem)     
