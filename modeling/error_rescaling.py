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
from pyLIMA import microltoolbox
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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
        
    (e,params) = create_model(e,params)
    
    norm_lcs = generate_residual_lcs(e,params)
    
    model_lc = pylima_lightcurve_tools.generate_model_lightcurve(e,diagnostics=True)
    
    estimate_err_scale_factors(norm_lcs,e,params)
    
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
            
            elif 'OUTPUT' in line:
                
                entries = line.replace('OUTPUT','').replace('\n','').split()
                params['output'] = entries[0]
                
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
    
    params['model_params'] = model_params
    
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

    plt.show()
    
    plt.close()
    
    return e,params

def generate_residual_lcs(e,params):
    
    f = e.fits[-1]
    model = f.model
    
    pylima_params = model.compute_pyLIMA_parameters(params['model_params'])
    
    norm_lcs = []
    
    for tel in e.telescopes:
        nlc = f.model_residuals(tel, pylima_params)
        
        norm_lcs.append(nlc)
    print(norm_lcs)
    
    return norm_lcs
    
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

def straightline(x, a0, a1):
    return a0 + a1*x
    
def estimate_err_scale_factors(norm_lcs,e,params):

    f = open(path.join(params['output'],'fitted_err_rescalings.dat'),'w')
    
    for i in range(0,len(norm_lcs),1):
        
        lc = norm_lcs[i]
        
        model_lc = pylima_lightcurve_tools.generate_model_lightcurve(e,ts=lc[:,0],
                                                                     diagnostics=True)
        
        N = float(len(lc))
        Ndof = N - len(params['model_params'])
        sqrtNdof = np.sqrt(Ndof/N)
        
        delta = (lc[:,1] - model_lc)
        
        idx1 = np.where(abs(delta) <= 0.1)[0]
        idx2 = np.where(lc[:,2] <= 0.05)[0]
        idx = list(set(idx1).intersection(set(idx2)))
        
        sigmas = lc[:,2]**2
        delta = (delta/sqrtNdof)**2
        
        (intercept,slope) = curve_fit(straightline, sigmas[idx], delta[idx])[0]
        
        if intercept > 0:
            a0 = np.sqrt(intercept)
            #a0 = intercept
        else:
            a0 = np.median(delta[idx])
        if slope > 0:
            a1 = np.sqrt(slope)
            #a1 = slope
        else:
            a1 = 0.0
        median_delta = np.median(np.sqrt(delta[idx]))
        
        f.write(params['datasets'][i].name+' a0='+str(a0)+' a1='+str(a1)+'\n')
        print(params['datasets'][i].name+\
            ' intercept='+str(intercept)+' slope='+str(slope)+\
            ' a0='+str(a0)+' a1='+str(a1)+' median delta='+str(median_delta))
        
        x = np.linspace(sigmas[idx].min(),sigmas[idx].max(),20)
        y = straightline(x, intercept, slope)
        
        fig = plt.figure(1,(10,10))
        
        plt.plot(sigmas,delta,'b.')
        
        plt.plot(x,y,'k-')
        
        [xmin,xmax,ymin,ymax] = plt.axis()
        plt.axis([xmin,xmax,-0.01,0.25])
        plt.xlabel('$\sigma^{2}$')
        plt.ylabel('$[(Data-model)/\sqrt{N_{dof}/N}]^{2}$')
        
        plt.savefig(path.join(params['output'],
                    'err_factor_fit_'+params['datasets'][i].name+'.png'))
    
        plt.close(1)
        
        err_new = np.sqrt(median_delta**2 + sigmas[idx])
        err_factor = np.median(err_new/np.sqrt(sigmas[idx]))
        f.write('Error factor = '+str(err_factor)+'\n')
        f.write('Median old error = '+str(np.median(np.sqrt(sigmas[idx])))+'\n')
        f.write('Median new error = '+str(np.median(err_new))+'\n')
        
    f.close()
    
if __name__ == '__main__':
    estimate_error_scaling()
    