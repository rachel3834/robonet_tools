# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 10:41:19 2019

@author: rstreet
"""

from os import path
from sys import argv
import pylima_lightcurve_tools
from pyLIMA import event
from pyLIMA import microlfits
from pyLIMA import telescopes
from pyLIMA import microlmodels
from pyLIMA import microlmagnification
from pyLIMA import microloutputs
import numpy as np
import matplotlib.pyplot as plt

class Dataset():
    
    def __init__(self, params=None):
        self.name = None
        self.filter = None
        self.data_file = None
        self.tel = None
        self.gamma = None
        
        if params != None:
            self.name = params['name']
            self.filter = params['filter']
            self.data_file = params['data_file']
    
    def summary(self):
        if self.gamma != None:
            return self.name+' '+self.filter+' '+self.data_file+' '+str(self.gamma)
        else:
            return self.name+' '+self.filter+' '+self.data_file

    def read_dataset_to_telescope(self):
        
        if 'p.t' in self.data_file:
            lightcurve = np.loadtxt(self.data_file,dtype=str)
            lightcurve = np.c_[lightcurve[:,1],lightcurve[:,6],lightcurve[:,7]].astype(float)
        
        if 'cal.t' in self.data_file:
            lightcurve = np.loadtxt(self.data_file,dtype=str)
            lightcurve = np.c_[lightcurve[:,1],lightcurve[:,8],lightcurve[:,9]].astype(float)
                
        if 'DK-1.54' in self.data_file:
            lightcurve = np.loadtxt(self.data_file,dtype=str)
            lightcurve = np.c_[lightcurve[:,1],lightcurve[:,6],lightcurve[:,7]].astype(float)
        
        self.tel = telescopes.Telescope(name=self.name, camera_filter=self.filter, 
                                 light_curve_magnitude=lightcurve,
                                 light_curve_magnitude_dictionnary={'time': 0, 'mag': 1, 'err_mag': 2})
        
def estimate_error_scaling():
    
    params = get_params()
    
    e = create_event(params)
    
    for d in params['datasets']:
        
        d.read_dataset_to_telescope()
    
        e.telescopes.append(d.tel)
        
    e = create_model(e,params)
    
    
def get_params():
    """Input file format expected:  ASCII, with lines:
    FIT_PARAMETERS: fit parameters as space-separated list in pyLIMA order on one line
    DATASET: Paths to a lightcurve datafile, one per line
    DATASET: Paths to a lightcurve datafile, one per line
    Datasets must be listed in the same order in which their fs, fb values appear
    in the fit parameters list.
    ...
    """
    
    params = {'datasets': []}    
    if len(argv) == 1:
        params['input_file'] = raw_input('Please enter the path to the input file: ')
    else:
        params['input_file'] = argv[1]

    if path.isfile(params['input_file']) == False:
        print('Error: Cannot find input file '+params['input_file'])
        exit()
    
    flines = open(params['input_file']).readlines()
    
    for line in flines:
        
        if '#' not in line[0:1]:
            
            if 'FIT_PARAMETER' in line:
                
                entries = line.replace('FIT_PARAMETER','').replace('\n','').split()
                if 'None' not in line:
                    params[entries[0].strip()] = float(entries[1])
                else:
                    params[entries[0].strip()] = None
            
            if 'MODEL_TYPE' in line:
                params['model_type'] = line.replace('MODEL_TYPE','').replace('\n','').strip()
                
            elif 'PHOT_PARAMETERS' in line:
                
                entries = line.replace('PHOT_PARAMETERS','').replace('\n','').split(',')
                ppars = []
                for value in entries:
                    ppars.append(float(value))
                params['phot_params'] = ppars
                
            elif 'CHISQ' in line:
                
                params['chisq'] = float(line.replace('CHISQ','').replace('\n',''))
                
            elif 'DATASET' in line:
                
                entries = line.replace('DATASET','').replace('\n','').split()
                d = Dataset({'name': entries[0], 'filter': entries[1], 'data_file': entries[2]})
                if len(entries) > 3:
                    d.gamma = float(entries[3])
                params['datasets'].append( d )
                
            elif 'EVENT_PARAMETER' in line:
                
                entries = line.replace('EVENT_PARAMETER','').replace('\n','').split()
                params[str(entries[0].strip()).lower()] = entries[1]
    
    for d in params['datasets']:
        print(d.summary())
        
    return params

def create_event(params):
    
    e = event.Event()
    
    e.name = params['name']
    e.ra = params['ra']
    e.dec = params['dec']    

    return e
    
def create_model(e,params):
    
    f = microlfits.MLFits(e)
    
    parallax_params = ['None', params['to']]
    orbital_params = ['None', params['to']]
    model_params = ['to', 'uo', 'tE']
    
    if params['rho'] != None:
        
        model_params.append('rho')
    
    if params['logs'] != None and params['logq'] != None and params['alpha'] != None:
        
        model_params.append('logs')
        model_params.append('logq')
        model_params.append('alpha')
        
    if params['piEN'] != None and params['piEE'] != None:
        
        parallax_params = ['Full', params['to']]
        model_params.append('piEN')
        model_params.append('piEE')

    if params['sdsdt'] != None and params['dalphadt'] != None and params['sdszdt'] != None:
        
        orbital_params = ['2D', params['to']]
        model_params.append('sdsdt')
        model_params.append('dalphadt')
        model_params.append('sdszdt')
            
    model = microlmodels.create_model(params['model_type'], e,
                                  parallax=parallax_params,
                                  orbital_motion=orbital_params,
                                  blend_flux_ratio = False)
    model.define_model_parameters()
    
    f.model = model
    
    results = []
    for key in model_params:
        results.append(params[key])
    
    results = results + params['phot_params']
    results.append(params['chisq'])
    
    f.fit_results = results
    
    e.fits.append(f)
    
    fig = microloutputs.LM_plot_lightcurves(f)
    
    #pylima_params = f.model.compute_pyLIMA_parameters(f.fit_results)
        
    #A = model.model_magnification(e.telescopes[0],pylima_params)
    
    #lightcurve = e.telescopes[0].lightcurve_magnitude
        
    #fig = plt.figure(1,(10,10))
    
    #plt.plot(lightcurve[:,0],A,'b.')
    
    #plt.xlabel('HJD')
    #plt.ylabel('Magnification')
    
    #rev_yaxis = False
    #if rev_yaxis:
    #    [xmin,xmax,ymin,ymax] = plt.axis()
    #    plt.axis([xmin,xmax,ymax,ymin])
    
    #plt.grid()
    #plt.savefig('lc_test.png')
    
    plt.show()
    
    return e
    
    
if __name__ == '__main__':
    estimate_error_scaling()
    