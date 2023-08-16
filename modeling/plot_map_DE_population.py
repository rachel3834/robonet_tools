# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 09:07:37 2018

@author: rstreet
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import os
from sys import argv
plt.rc('font', family='DejaVu Sans')
plt.rc('xtick', labelsize='large')
plt.rc('ytick', labelsize='large')


def map_DE_population(input_file):
    """Function to plot maps of the (logs .vs. logq) parameter space mapped
    out by a Differential Evolution (DE) algorithm search.  
    
    The first plot is a 2D histogram of the number of DE samples
    for each pixel in (logs .vs. logq) space, as a proxy for a map of the 
    population of solutions calculated during a Differential Evolution 
    search of parameter space.
    
    The second plot is a 2D histogram of the chi squared values of each pixel
    in the parameter space.
    
    Based on original code by E. Bachelet.
    """
    
    map_data = np.load(input_file)
    order = 'descending'
    # Descending order sort
    if order == 'descending':
        map_data_sort = map_data[map_data[:,-1].argsort()[::-1][:len(map_data[:,-1])]]
    else:
        map_data_sort = map_data[map_data[:, -1].argsort()]

    fig, (ax1, ax2) = plt.subplots(1, 2)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)

    # Down sample the number of map points to speed up plotting
    #index = np.random.randint(0,len(map_data),int(len(map_data)*0.12))
    index = np.arange(0,len(map_data_sort),1, dtype=int)

    h = ax1.hist2d(map_data_sort[index,4],map_data_sort[index,5],
               norm=LogNorm(),bins=(30,30))
    
    ax1.set_title('N DE samples')
    ax1.set_xlabel('$log_{10}(s)$')
    ax1.set_ylabel('$log_{10}(q)$')
    plt.colorbar(h[3], ax=ax1)

    im2 = ax2.scatter(map_data_sort[index,4],map_data_sort[index,5],
                c=np.log10(map_data_sort[index,-1]))

    ax2.set_title('$\log_{10}(\chi^{2})$')
    ax2.set_xlabel('$log_{10}(s)$')
    plt.colorbar(im2, ax=ax2)

    plt.savefig('DE_population_maps.png')
    
    plt.close(1)
    
if __name__ == '__main__':
    
    if len(argv) == 1:
        input_file = raw_input('Please enter the path to the DE population output: ')
        
    else:
        input_file = argv[1]
    
    map_DE_population(input_file)
    