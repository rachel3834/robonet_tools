# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 17:37:03 2018

@author: rstreet
"""

from sys import argv
from os import path
import numpy as np

def extract_best_fit_mcmc_parameters(input_file):
    """Function to identify the best fitting model from an output file of 
    MCMC chain data, estimate the uncertainties at 16% and 84% and calculate
    the model chi squared from the maximum likelihood data.
    
    Based on code from pyLIMA by E. Bachelet.
    """
    
    model_parameters = ['t0', 'u0', 'tE', 'rho', 
                        'logs', 'logq', 'alpha', 'piEN', 'piEE',
                        '1/sdsdt', 'dalpha/dt', '1/sdsz/dt']
    
    if path.isfile(input_file) == False:
        raise IOError('Cannot find input file '+input_file)
        exit()
        
    mcmc_chains = np.loadtxt(input_file)
    
    n_params = len(mcmc_chains[0,:]) - 1
    
    best_model_index = np.argmax(mcmc_chains[:, -1])

    for index, key in enumerate(model_parameters):
        if index < n_params:
            best_param = mcmc_chains[best_model_index, index]
            percent_34 = np.percentile(mcmc_chains[:, index], 16)
            percent_50 = np.percentile(mcmc_chains[:, index], 50)
            percent_84 = np.percentile(mcmc_chains[:, index], 84)
        
            print(key + ' ' + str(best_param) + ' [' + str(percent_34) + ',' + str(percent_50) + ',' + str(
                    percent_84) + ']')
                    
    print('Chi^2 ' + str(-2*mcmc_chains[best_model_index, -1]) + '&0\n')
    
if __name__ == '__main__':
    
    if len(argv) == 1:
        input_file = raw_input('Please enter the path to the file containing the MCMC chains data: ')
        
    else:
        input_file = argv[1]
    
    extract_best_fit_mcmc_parameters(input_file)
    
    