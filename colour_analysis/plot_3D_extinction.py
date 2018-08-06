# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 15:04:19 2018

@author: rstreet
"""

from os import path
from sys import argv
import numpy as np
import matplotlib.pyplot as plt

def plot_3D_extinction_data():
    """Function to plot the 3D extinction derived from the maps from Pan-STARRS1
    by Green et al. 2015ApJ...810...25G"""
    
    if len(argv) == 1:
        
        data_file = raw_input('Please enter the path to the data file: ')
        plot_file = raw_input('Please enter the path to the output plot file: ')
        DM_source = float(raw_input('Please enter the distance modulus of the source: '))
        
    else:
        
        data_file = argv[1]
        plot_file = argv[2]
        DM_source = float(argv[3])
        
    (DistMod,EBV) = read_3D_map_data(data_file)
    
    plot_EBV_distance(DistMod, EBV, DM_source, plot_file)
    

def read_3D_map_data(data_file):
    """Function to read the results of a query through
    http://argonaut.skymaps.info/
    """
    
    if path.isfile(data_file) == False:
        print('ERROR: Cannot find input file '+data_file)
        exit()
        
    line_list = open(data_file,'r').readlines()
    
    DistMod = []
    EBV = []
    
    for i in range(0,len(line_list),1):
        
        line = line_list[i]
        
        if line[0:1] != '#':
            
            if 'SampleNo | DistanceModulus' in line:
                
                columns = line_list[i+1]
                
                for item in columns.replace('|','').replace('\n','').split():
                    
                    if len(item) > 0:
                        DistMod.append(float(item))
            
                i += 1
            
            if 'BestFit' in line:
                
                columns = line.replace('|','').replace('\n','').split()
                
                for item in columns[1:]:
                    
                    if len(item) > 0:
                        EBV.append(float(item))
    
    DistMod = np.array(DistMod)
    EBV = np.array(EBV)
    
    return DistMod, EBV

def plot_EBV_distance(DistMod, EBV, DM_source, plot_file):
    """Function to plot the extinction as a function of distance along the
    line of sight"""
    
    fig = plt.figure(1,(10,10))

    plt.rcParams.update({'font.size': 18})
    
    plt.plot(DistMod,EBV,'k-')
    
    [xmin,xmax,ymin,ymax] = plt.axis()
    
    if DM_source < xmax:
        ydata = np.arange(EBV.min(), EBV.max()+0.2, 0.5)    
        xdata = np.zeros(len(ydata))
        xdata.fill(DM_source)
    
        plt.plot(xdata, ydata,'m-.')

    else:
        
        plt.arrow(xmax-1.0, (ymax-ymin)/2.0, 1.0, 0.0)
        plt.text('$source_{m-M}$ = '+str(round(DM_source,1))+' mag', 
                 xmax-1.0, (ymax-ymin)/2.0)
    plt.xlabel('Distance Modulus [mag]')
    plt.ylabel('E(B-V) [mag]')

    plt.grid()
    
    plt.savefig(plot_file)
    
if __name__ == '__main__':

    plot_3D_extinction_data()