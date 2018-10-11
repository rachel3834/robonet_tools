# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 17:37:03 2018

@author: rstreet
"""

from sys import argv
from os import path
import numpy as np
import calc_sq_with_errors

def extract_best_fit_mcmc_parameters(input_file):
    """Function to identify the best fitting model from an output file of 
    MCMC chain data, estimate the uncertainties at 16% and 84% and calculate
    the model chi squared from the maximum likelihood data.
    
    Based on code from pyLIMA by E. Bachelet.
    """
    
    if path.isfile(input_file) == False:
        raise IOError('Cannot find input file '+input_file)
        exit()
        
    mcmc_chains = np.loadtxt(input_file)
    
    n_params = len(mcmc_chains[0,:]) - 1
    
    if n_params == 7:
        model_parameters = ['t0', 'u0', 'tE', 'rho', 
                            'logs', 'logq', 'alpha']
    elif n_params == 9:
        model_parameters = ['t0', 'u0', 'tE', 'rho', 
                        'logs', 'logq', 'alpha', 'piEN', 'piEE']
    elif n_params == 11:
        model_parameters = ['t0', 'u0', 'tE', 'rho', 
                        'logs', 'logq', 'alpha', 'piEN', 'piEE',
                        'dsdt', 'dalpha/dt']
    elif n_params == 12:
        model_parameters = ['t0', 'u0', 'tE', 'rho', 
                        'logs', 'logq', 'alpha', 'piEN', 'piEE',
                        '1/sdsdt', 'dalpha/dt', '1/sdsz/dt']
    
    best_model_index = np.argmax(mcmc_chains[:, -1])

    for index, key in enumerate(model_parameters):
        if index < n_params:
            best_param = mcmc_chains[best_model_index, index]
            percent_34 = np.percentile(mcmc_chains[:, index], 16)
            percent_50 = np.percentile(mcmc_chains[:, index], 50)
            percent_84 = np.percentile(mcmc_chains[:, index], 84)
            
            lower = percent_50 - percent_34
            upper = percent_84 - percent_50
            
            if 'log' in key:
                (best_param10, lower10, upper10) = calc_sq_with_errors.calc_sq_with_errors(best_param, lower, upper)
                
            if key in ['t0']:
                ndp = 5
            elif key in ['tE']:
                ndp = 3
            else:
                ndp = 6
            
            best_param = round(best_param,ndp)
            percent_34 = round(percent_34,ndp)
            percent_50 = round(percent_50,ndp)
            percent_84 = round(percent_84,ndp)
            lower = round(lower,ndp)
            upper = round(upper,ndp)
            
            if 'log' in key:
                best_param10 = round(best_param10,ndp)
                lower10 = round(lower10,ndp)
                upper10 = round(upper10,ndp)
            
            print(key + ' ' + str(best_param) + ' [' + str(percent_34) + ',' + str(percent_50) + ',' + str(
                    percent_84) + ']')
            print(key + ' ' + str(best_param) + '$_{-' + str(lower) + '}^{+' + str(upper)+'}$')
            if 'log' in key:
                print(key.replace('log','') + ' ' + str(best_param10) + '$_{-' + str(lower10) + '}^{+' + str(upper10)+'}$')
                
    chichi = -2*mcmc_chains[best_model_index, -1]
    print('Chi^2 ' + str(round(chichi,3)) + '\n')
    
if __name__ == '__main__':
    
    if len(argv) == 1:
        input_file = raw_input('Please enter the path to the file containing the MCMC chains data: ')
        
    else:
        input_file = argv[1]
    
    extract_best_fit_mcmc_parameters(input_file)
    
    