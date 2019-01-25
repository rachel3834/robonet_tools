# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 13:26:44 2018

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

def map_mcmc_chains(input_file):
    """Function to plot maps of the (logs .vs. logq) parameter space mapped
    out by a Markov Chain Monte Carlo (MCMC) algorithm search.  
    """
    
    mcmc_chains = np.loadtxt(input_file)
    mcmc_chains = mcmc_chains[mcmc_chains[:,-1].argsort(),]
    
    
    fig = plt.figure(1,(10,5))
        
    plt.subplot(121)
    plt.subplots_adjust(wspace=0.5)
                        
    # Down sample the number of map points to speed up plotting
    #index = np.random.randint(0,len(map_data),int(len(map_data)*0.12))
    index = np.arange(0,len(mcmc_chains),1, dtype=int)
    
    chisq = -2*mcmc_chains[index, -1]
    

    plt.scatter(mcmc_chains[index,4],mcmc_chains[index,5],
                c=mcmc_chains[index,-1],alpha=0.25)
    
    plt.title('$\chi^{2}$')
    plt.xlabel('$log_{10}(s)$')
    plt.ylabel('$log_{10}(q)$')
    plt.colorbar()

    plt.subplot(122)
    
    plt.scatter(mcmc_chains[index,7],mcmc_chains[index,8],
                c=mcmc_chains[index,-1],alpha=0.25)
    
    plt.title('$\chi^{2}$')
    plt.xlabel('$\pi_{E,N}$')
    plt.ylabel('$\pi_{E,E}$')
    plt.colorbar()

    plt.savefig('mcmc_chains_map.png')
    
    plt.close(1)
    
if __name__ == '__main__':
    
    if len(argv) == 1:
        input_file = raw_input('Please enter the path to the MCMC chains data file: ')
        
    else:
        input_file = argv[1]
    
    map_mcmc_chains(input_file)
    