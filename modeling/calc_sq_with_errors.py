# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 09:37:11 2018

@author: rstreet
"""
from sys import argv

def calc_sq_with_errors(log_param, lower_range, upper_range):
    """Function to convert a parameter computed in log_10 (e.g. logs, logq)
    to its un-logged value and calculate appropriate positive and negative 
    errorbars from upper and lower range values
    """
    
    best_param10 = 10**log_param
    upper10 = 10**(log_param+upper_range)
    lower10 = 10**(log_param-lower_range)
    upper10 = upper10 - best_param10
    lower10 = best_param10 - lower10
    
    return best_param10, lower10, upper10

if __name__ == '__main__':
    
    if len(argv) == 1:
        log_param = float(raw_input('Please enter the log10 parameter value: '))
        lower_range = float(raw_input('Please enter the lower range value: '))
        upper_range = float(raw_input('Please enter the upper range value: '))
        
    else:
        log_param = float(argv[1])
        lower_range = float(argv[2])
        upper_range = float(argv[3])
    
    (best_param10, lower10, upper10) = calc_sq_with_errors(log_param, lower_range, upper_range)
    
    print(str(round(best_param10,6))+'_-'+str(round(lower10,6))+'^+'+str(round(upper10,6)))
    