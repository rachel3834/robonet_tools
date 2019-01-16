# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 10:09:36 2019

@author: rstreet
"""

import numpy as np
import copy
from pyLIMA import microltoolbox
import matplotlib.pyplot as plt

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
