# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 10:09:36 2019

@author: rstreet
"""

from sys import argv
from os import path
import numpy as np
import copy
from pyLIMA import microltoolbox
import matplotlib.pyplot as plt
import data_handling_utils

def plot_event_lightcurve():
    """Function to plot an event lightcurve with a model overplotted from
    pyLIMA"""
    
    params = get_params()

    current_event = generate_event(params)
    
    for i,data_file in enumerate(params['datasets']):
        
        gamma = fetch_gamma(os.path.basename(data_file))
        
        rescaling = fetch_rescaling_factors(data_file,params)
        
        tel = data_handling_utils.read_data_file(data_file,calibrated=False, 
                                                 gamma=gamma, rescaling=rescaling)
        
        current_event.telescopes.append(telescope)
    
    current_event.find_survey(params['survey'])
    current_event.check_event()


def get_params():
    
    params = {}
    
    if len(argv) == 1:
        
        input_file = raw_input('Please enter the path to the parameter file: ')

    else:

        input_file = argv[1]
    
    if path.isfile(input_file) == False:
        
        print('ERROR: Cannot find input parameter file')
        exit()
        
    flines = open(input_file,'r').readlines()
    
    datasets = []
    
    for line in flines:
        (key,value) = line.replace('\n','').split()
        
        key = str(key).lower()
        
        if 'dataset' in key:
            datasets.append(value)
        
        if key in ['output','error_scaling_file','survey']:
            params[str(key).lower()] = value
        
        if key in ['ra', 'dec']:
            params[str(key).lower()] = float(value)
            
    params['datasets'] = datasets
    
    if 'error_scaling_file' in params.keys():

        rescaling_coeffs = read_error_rescaling_factors(params['error_scaling_file'])
        
        for key, value in rescaling_coeffs.items():
            
            params[key] = value
            
    return params
    
def read_error_rescaling_factors(coeffs_file):
    
    rescaling_coeffs = {}
    
    if path.isfile(coeffs_file):
        print('Read error rescaling factors:')
        
        flines = open(coeffs_file,'r').readlines()
        
        for line in flines:
            if line[0:1] != '#' and len(line.replace('\n','')) > 0:
                
                
                (tel_code, a0, a1) = line.replace('\n','').split()
                
                rescaling_coeffs[tel_code] = [float(a0), float(a1)]
                print(tel_code+' a0 = '+a0+', a1='+a1)
        
    elif path.isfile(coeffs_file) == False:
        print('ERROR: Error rescaling is switched on but no '+coeffs_file+' found')
        exit()
        
    return rescaling_coeffs

def fetch_rescaling_factors(data_file,params):
    
    if 'DK-1.54' in data_file:
        name = 'Danish'
    else:
        name = data_file[14:17]+data_file[21]+'_'+data_file[-4:-2]
    
    if name in params.keys():
        factors = params['name']
    else:
        factors = []
        
    return factors
        
def fetch_gamma(lightcurve_file):

    gamma_coeffs = { 'gp': 0.8371,
                     'rp': 0.6445,
                     'ip': 0.503,
                     'Z': 0.4134 }

    gamma = None

    for f in gamma_coeffs.keys():

        if '_'+f in lightcurve_file:
            gamma = gamma_coeffs[f]
            #print(lightcurve_file+' -> filter = '+f+' gamma = '+str(gamma))

    if gamma == None:
        print('ERROR: No gamma value available for lightcurve '+lightcurve_file)
        exit()

    return gamma

def generate_event(params):
    """Function to initialise a pyLIMA event object"""
    
    current_event = event.Event()
    current_event.name = params['name']
    
    current_event.ra = params['ra']
    current_event.dec = params['dec']
    
    return current_event
    
def generate_model_lightcurve(e,ts=None,diagnostics=False):
    """Function to produce a model lightcurve based on a parameter set
    fitted by pyLIMA
    
    Inputs:
    e  Event object, with attributed lightcurve data and model fit(s)
    """
    
    lc = e.telescopes[0].lightcurve_magnitude
    
    fit_params = e.fits[-1].model.compute_pyLIMA_parameters(e.fits[-1].fit_results)
    
    if type(ts) != type(np.zeros(1)):
        ts = np.linspace(lc[:,0].min(), lc[:,0].max(), len(lc[:,0]))

    reference_telescope = copy.copy(e.fits[-1].event.telescopes[0])
    
    reference_telescope.lightcurve_magnitude = np.array([ts, [0] * len(ts), [0] * len(ts)]).T
    
    reference_telescope.lightcurve_flux = reference_telescope.lightcurve_in_flux()

    if e.fits[-1].model.parallax_model[0] != 'None':
        
        reference_telescope.compute_parallax(e.fits[-1].event, e.fits[-1].model.parallax_model)

    flux_model = e.fits[-1].model.compute_the_microlensing_model(reference_telescope, fit_params)[0]
    
    mag_model = microltoolbox.flux_to_magnitude(flux_model)

    if diagnostics:
        fig = plt.figure(1,(10,10))
    
        plt.plot(ts,mag_model,'r-')
    
        plt.xlabel('HJD')
        plt.ylabel('Magnitude')
        
        rev_yaxis = True
        if rev_yaxis:
            [xmin,xmax,ymin,ymax] = plt.axis()
            plt.axis([xmin,xmax,ymax,ymin])
        
        plt.grid()
        plt.savefig('lc_model_test.png')
        
        plt.close(1)
        
    return mag_model
    

if __name__ == '__main__':
    
    plot_event_lightcurve()
    