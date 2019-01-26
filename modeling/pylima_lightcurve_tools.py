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
from pyLIMA import event
from pyLIMA import microlfits
from pyLIMA import microlmodels

def plot_event_lightcurve():
    """Function to plot an event lightcurve with a model overplotted from
    pyLIMA"""
    
    params = get_params()
    
    rescaling_coeffs = read_error_rescaling_factors(params['error_scaling_file'])
    
    params = read_data_files(params, rescaling_coeffs)
    
    current_event = create_event(params)
    
    for d in params['data']:
        
        current_event.telescopes.append(d.tel)
    
    current_event.find_survey(params['survey'])
    current_event.check_event()

    (current_event,params) = create_model(current_event,params)
    
    residual_lcs = generate_residual_lcs(current_event,params)
    
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
    phot_params = {}
    
    for line in flines:
        key = line.replace('\n','').split()[0]
        key = str(key).lower()
        
        if 'dataset' in key:
            entries = line.replace('DATASET','').replace('\n','').split()
            ddict = {'name': entries[0], 'filter': entries[1], 
                     'data_file': entries[2], 'gamma': float(entries[3])}
            datasets.append(ddict)
        
        if key in ['name','output','error_scaling_file','survey','model_type']:
            entries = line.replace(key.upper(),'').replace('\n','').split()
            params[key] = entries[0]
        
        if key in ['ra', 'dec', 'to', 'uo', 'te', 'rho', 'logs', 'logq', \
                    'alpha', 'pien', 'piee', 'dsdt', 'dalphadt']:
            entries = line.replace(key.upper(),'').replace('\n','').split()
            params[key] = float(entries[0])
        
        if 'fs_' in key or 'fb_' in key:
            
            name = key.replace('fs_','').replace('fb_','')
            value = float(line.replace('\n','').split()[-1])
            
            if name in phot_params.keys():
                p = phot_params[name]
            else:
                p = {}
            
            if 'fs_' in key:
                p['fs'] = value
            else:
                p['fb'] = value
            
            phot_params[name] = p
            
    params['datasets'] = datasets
    params['phot_params'] = phot_params
    
    return params

def read_data_files(params, rescaling_coeffs):
    
    datasets = []
    
    for ddict in params['datasets']:
        
        d = data_handling_utils.Dataset(ddict)
        
        if 'gamma' in ddict.keys():
            d.gamma = ddict['gamma']
        
        rescale_factors = fetch_rescaling_factors(d.name,rescaling_coeffs)
        
        d.read_dataset_to_telescope(params['model_type'], 
                                    rescaling=rescale_factors)
        
        datasets.append( d )
        
        print(d.summary())
        
    params['data'] = datasets
    
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

def fetch_rescaling_factors(name,rescaling_coeffs):
    
    if name in rescaling_coeffs.keys():
        factors = rescaling_coeffs[name]
    else:
        factors = []

    return factors
        
def create_event(params):
    """Function to initialise a pyLIMA event object"""
    
    current_event = event.Event()
    current_event.name = params['name']
    
    current_event.ra = params['ra']
    current_event.dec = params['dec']
    
    return current_event

def create_model(current_event,params,diagnostics=False):
    
    f = microlfits.MLFits(current_event)
    
    parallax_params = ['None', params['to']]
    orbital_params = ['None', params['to']]
    model_params = ['to', 'uo', 'tE']
    
    if params['rho'] != None:
        
        model_params.append('rho')
    
    if params['logs'] != None and params['logq'] != None and params['alpha'] != None:
        
        model_params.append('logs')
        model_params.append('logq')
        model_params.append('alpha')
        
    if params['pien'] != None and params['piee'] != None:
        
        parallax_params = ['Full', params['to']]
        model_params.append('piEN')
        model_params.append('piEE')

    if 'sdsdt' in params.keys() and 'dalphadt' in params.keys() and 'sdszdt' in params.keys():
        
        orbital_params = ['3D', params['to']]
        model_params.append('sdsdt')
        model_params.append('dalphadt')
        model_params.append('sdszdt')
    
    elif 'dsdt' in params.keys() and 'dalphadt' in params.keys():
        
        orbital_params = ['2D', params['to']]
        model_params.append('dsdt')
        model_params.append('dalphadt')
        
    params['model_params'] = model_params
    
    model = microlmodels.create_model(params['model_type'], current_event,
                                  parallax=parallax_params,
                                  orbital_motion=orbital_params,
                                  blend_flux_ratio = False)
    
    if 'binary origin' in params.keys():
        model.binary_origin  = params['binary_origin']
        
    model.define_model_parameters()
    
    f.model = model
    
    results = []
    for key in model_params:
        results.append(params[str(key).lower()])
    
    params['fitted_model_params'] = results
    
    results = results + params['phot_params']
    results.append(params['chisq'])
    
    f.fit_results = results
    
    current_event.fits.append(f)
    
    if diagnostics:
        fig = microloutputs.LM_plot_lightcurves(f)
    
        plt.show()
        
        plt.close()
    
    return current_event,params

def generate_residual_lcs(current_event,params):
    
    f = e.fits[-1]
    model = f.model
    
    pylima_params = model.compute_pyLIMA_parameters(params['fitted_model_params'])
    
    residual_lcs = []
    
    for tel in current_event.telescopes:
        
        flux = tel.lightcurve_flux[:, 1]
        flux_model = model.compute_the_microlensing_model(tel, pylima_params)[0]
        
        dmag = 2.5 * np.log10(flux_model/flux)
        
        res_lc = np.copy(tel.lightcurve_magnitude)
        res_lc[:,1] = dmag
        
        residual_lcs.append(res_lc)
    
    return residual_lcs
    
def plot_lcs(lc_data):

    fig = plt.figure(1,(10,10))
    
    for i,lc in enumerate(lc_data):
        plt.errorbar(lc[:,0],lc[:,1],yerr=lc[:,2],fmt='.')
    
    plt.xlabel('HJD')
    plt.ylabel('Magnitude')
    
    rev_yaxis = True
    if rev_yaxis:
        [xmin,xmax,ymin,ymax] = plt.axis()
        plt.axis([xmin,xmax,ymax,ymin])
    
    plt.grid()
    plt.savefig('lc_test.png')
    
    plt.close(1)

 
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
    