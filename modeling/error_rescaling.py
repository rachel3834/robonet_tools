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
from scipy import optimize

class Dataset():
    
    def __init__(self, params=None):
        self.name = None
        self.filter = None
        self.data_file = None
        self.tel = None
        self.gamma = None
        self.k = None
        self.k_err = None
        self.emin = None
        self.emin_err = None
        self.a0 = None
        self.a0_err = None
        self.a1 = None
        self.a1_err = None
        
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
        
    def k_emin_results(self,latex=False):
        if latex and self.k_err != None and self.emin_err != None:
            return self.name+' k='+str(round(self.k,3))+'$\pm$'+str(round(self.k_err,3))+\
                        ', emin='+str(round(self.emin,3))+'$\pm$'+str(round(self.emin_err,3))
        else:
            return self.name+' k='+str(self.k)+'+/-'+str(self.k_err)+', emin='+str(self.emin)+'+/-'+str(self.emin_err)
    
    def slope_results(self,latex=False):
        if latex and self.a0_err != None and self.a1_err != None:
            return self.name+' a0='+str(round(self.a0,3))+'$\pm$'+str(round(self.a0_err,3))+\
                    ', a1='+str(round(self.a1,3))+'$\pm$'+str(round(self.a1_err,3))
        else:
            return self.name+' a0='+str(self.a0)+'+/-'+str(self.a0_err)+', a1='+str(self.a1)+'+/-'+str(self.a1_err)
        
def estimate_error_scaling():
    
    params = get_params()
    
    e = create_event(params)
    
    for d in params['datasets']:
        
        d.read_dataset_to_telescope()
    
        e.telescopes.append(d.tel)
        
    (e,params) = create_model(e,params)
    
    residual_lcs = generate_residual_lcs(e,params)
        
    params = estimate_err_scale_factors(residual_lcs,e,params)
    
    params = rescale_lightcurves(e,params)
    
    output_results(params)
    
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
                
            elif 'BINARY_ORIGIN' in line:
                
                entries = line.replace('BINARY_ORIGIN','').replace('\n','').split()
                params['binary_origin'] = entries[0]
                
    for d in params['datasets']:
        print(d.summary())
        
    return params

def create_event(params):
    
    e = event.Event()
    
    e.name = params['name']
    e.ra = params['ra']
    e.dec = params['dec']    

    return e
    
def create_model(e,params,diagnostics=False):
    
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
    
    if 'binary origin' in params.keys():
        model.binary_origin  = params['binary_origin']
        
    model.define_model_parameters()
    
    f.model = model
    
    results = []
    for key in model_params:
        results.append(params[key])
    
    params['fitted_model_params'] = results
    
    results = results + params['phot_params']
    results.append(params['chisq'])
    
    f.fit_results = results
    
    e.fits.append(f)
    
    if diagnostics:
        fig = microloutputs.LM_plot_lightcurves(f)
    
        plt.show()
        
        plt.close()
    
    return e,params

def generate_residual_lcs(e,params):
    
    f = e.fits[-1]
    model = f.model
    
    pylima_params = model.compute_pyLIMA_parameters(params['fitted_model_params'])
    
    residual_lcs = []
    
    for tel in e.telescopes:
        
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

def straightline(x, a0, a1):
    return a0 + a1*x
    
def estimate_err_scale_factors(residual_lcs,e,params):

    print('\nEstimating rescaling factors using slope approach: ')
    
    for i in range(0,len(residual_lcs),1):
        
        d = params['datasets'][i]

        lc = residual_lcs[i]
        
        N = float(len(lc))
        Ndof = N - len(params['model_params'])
        sqrtNdof = np.sqrt(Ndof/N)
        Ndofsq = (Ndof/N)**2
        
        delta = lc[:,1]
        
        idx1 = np.where(abs(delta) <= 0.1)[0]
        idx2 = np.where(lc[:,2] <= 0.05)[0]
        idx = list(set(idx1).intersection(set(idx2)))
        
        sigmas_sq = lc[:,2]**2
        res_sq = (delta/Ndofsq)**2
        
        #(intercept,slope) = optimize.curve_fit(straightline, sigmas[idx], delta[idx])[0]
        (popt,pcov) = optimize.curve_fit(straightline, sigmas_sq[idx], res_sq[idx])
        
        (intercept,slope) = popt
        (intercept_err,slope_err) = np.sqrt(np.diag(pcov))
        
        if intercept > 0:
            d.a0 = np.sqrt(intercept)
            d.a0_err = (intercept_err/intercept)*d.a0
        else:
            median_res = np.median(res_sq[idx])
            mad = np.median(res_sq[idx]-median_res)
            
            d.a0 = np.sqrt(median_res)
            d.a0_err = (mad/median_res)*d.a0
            
        if slope > 0:
            d.a1 = np.sqrt(slope)
            d.a1_err = (slope_err/slope)*d.a1
        else:
            d.a1 = 1.0
            d.a1_err = None
        
        
        res = np.sqrt(res_sq[idx])
        median_res = np.median(np.sqrt(res))
        mad_res = abs(np.median(res-median_res))
        
        print(d.slope_results()+' median residuals ='+\
                str(median_res)+'+/-'+str(mad_res))
        
        x = np.linspace(sigmas_sq[idx].min(),sigmas_sq[idx].max(),20)
        y = straightline(x, intercept, slope)
        
        fig = plt.figure(1,(10,10))
        
        plt.plot(sigmas_sq,res_sq,'b.')
        
        plt.plot(x,y,'k-')
        
        [xmin,xmax,ymin,ymax] = plt.axis()
        plt.axis([xmin,xmax,-0.01,0.25])
        plt.xlabel('$\sigma^{2}$')
        plt.ylabel('$[(Data-model)/\sqrt{N_{dof}/N}]^{2}$')
        plt.title('Residuals for '+d.name)
        
        plt.savefig(path.join(params['output'],
                    'err_factor_fit_'+d.name+'.png'))
    
        plt.close(1)
        
        err_new = np.sqrt(median_res**2 + sigmas_sq[idx])
        err_factor = np.median(err_new/np.sqrt(sigmas_sq[idx]))
        
        params['datasets'][i] = d
    
    return params
    
def objective_function(fit_process_parameters, your_event, your_model, guess):
    """Optimizing function for errorbar rescaling by Etienne Bachelet"""
    
    rescaled_flux = np.copy(fit_process_parameters)
    model_param = np.copy(guess)
    
    pyLIMA_parameters = your_model.compute_pyLIMA_parameters(model_param)
    
    chichi = 0
    for index,telescope in enumerate(your_event.telescopes):
         if (rescaled_flux[2*index]>0) and (rescaled_flux[2*index+1]>=0):
              pass
         else:
              return np.inf
              
    for index,telescope in enumerate(your_event.telescopes):
        # Find the residuals of telescope observation regarding the parameters and model
        model = your_model.compute_the_microlensing_model(telescope, pyLIMA_parameters)
        flux= telescope.lightcurve_flux[:,1]
        
        errmag = np.sqrt(rescaled_flux[2*index]**2*telescope.lightcurve_magnitude[:,2]**2 + \
                            rescaled_flux[2*index+1]**2)
        errflux = errmag*telescope.lightcurve_flux[:,1]*np.log(10)/2.5
        
        residus = (flux - model[0])/errflux 
        chichi += (residus ** 2).sum()+np.sum(2*np.log(errflux))
    
        #print('CHICHI: ',chichi, fit_process_parameters, guess)
    
    #opt = raw_input('Continue?')
    
    return chichi

def rescale_lightcurves(e,params):
    """Function to estimate the re-scaling factors k and emin from the 
    equation:
    
    e' = k sqrt(sigma**2 + emin**2)
    
    Developed by Etienne Bachelet
    """
    
    print('\nEstimating rescaling factors using k,emin approach:')
    rescale_params_guess = []

    for i in e.telescopes:
        
        rescale_params_guess.append(1.0)    # k
        rescale_params_guess.append(0.0)    # emin
        
    fit_model_params = []
    for key in params['model_params']:
        fit_model_params.append(params[key])
    
    result = optimize.minimize(objective_function, rescale_params_guess,
                               args=(e, e.fits[-1].model, fit_model_params),
                               options={'maxiter': 1e5})
    
    covar = result.hess_inv
    
    x_errors = np.sqrt( np.diagonal(covar) )
    
    for i, tel in enumerate(e.telescopes):
        
        d = params['datasets'][i]
        d.k = result.x[i*2]
        d.k_err = x_errors[i*2]
        d.emin = result.x[i*2+1]
        d.emin_err = x_errors[i*2+1]
        params['datasets'][i] = d
        
        print(d.k_emin_results())
        
    return params

def output_results(params):
    
    tel_codes = []
    for d in params['datasets']:

        tel_codes.append(d.name)
    
    tel_codes = np.array(tel_codes)
    idx = np.argsort(tel_codes)
    
    output_file = path.join(params['output'],'fitted_err_rescalings.dat')
    f = open(output_file,'w')
    
    f.write('Results from slope method:\n')
    for i in idx:
        
        d = params['datasets'][i]        
        f.write(d.slope_results()+'\n')
    
    f.write('\n')
    for i in idx:
        
        d = params['datasets'][i]    
        f.write(d.slope_results(latex=True)+'\n')
    
    f.write('\n')
    f.write('Results from k/emin method:\n')
    for i in idx:
        
        d = params['datasets'][i]  
        f.write(d.k_emin_results()+'\n')
    
    f.write('\n')
    for i in idx:
        
        d = params['datasets'][i]  
        f.write(d.k_emin_results(latex=True)+'\n')
        
    f.close()
    
    print('\nOutput summary of results to '+output_file)
    
if __name__ == '__main__':
    estimate_error_scaling()
    