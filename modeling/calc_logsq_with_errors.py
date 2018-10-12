# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 16:56:01 2018

@author: rstreet
"""

from sys import argv
import numpy as np

def calc_logsq_with_errors(param, sig_param):
    """Function to convert a parameter computed in log_10 (e.g. logs, logq)
    to its un-logged value and calculate appropriate positive and negative 
    errorbars from upper and lower range values
    """
    
    log_param = np.log10(param)
    
    sig_log_param = np.sqrt((sig_param/param)**2) * abs(log_param)
    
    return log_param, sig_log_param

if __name__ == '__main__':
    
    if len(argv) == 1:
        param = float(raw_input('Please enter the parameter value: '))
        sig_param = float(raw_input('Please enter the uncertainty on the parameter value: '))
        
    else:
        param = float(argv[1])
        sig_param = float(argv[2])
    
    (log_param, sig_log_param) = calc_logsq_with_errors(param, sig_param)
    
    print(str(round(log_param,6))+'+/-'+str(round(sig_log_param,6)))
    