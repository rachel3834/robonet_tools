# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 10:25:19 2019

@author: rstreet
"""

from sys import argv
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='DejaVu Sans')
plt.rc('xtick', labelsize='large')
plt.rc('ytick', labelsize='large')

def plot_chain_distributions(input_file):
    
    mcmc_chains = np.loadtxt(input_file)
        
    fig = plt.figure(1,(7,5))
        
    #plt.hist(mcmc_chains[:,2],bins=50)
    
    plt.plot(mcmc_chains[:,2],mcmc_chains[:,-1]*-2.0,'r.')
    
    plt.xlabel('$t_{E}$ [days]')
    plt.ylabel('$\chi^{2}$')
    
    plt.show()
    #plt.savefig('mcmc_chains_map.png')
    
    plt.close(1)
    

if __name__ == '__main__':
    
    if len(argv) == 1:
        input_file = raw_input('Please enter the path to the file containing the MCMC chains data: ')
        
    else:
        input_file = argv[1]
    
    plot_chain_distributions(input_file)
    
    