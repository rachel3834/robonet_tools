# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 13:16:17 2018

@author: rstreet
"""

from sys import argv
import numpy as np
from scipy import optimize

def transform(x,C,a,b,c,d):
    """General transformation function as defined by 
    Jordi, C et al, (2010), A&A, 523, A48:
    
    C1 = a + bC2 + cC2^2 + dC2^3

    For solution using Halley's method, reformulated to find the roots of the
    equation:
    
    f eqv C1 - a - bC2 + cC2^2 + dC2^3,
    
    where x = C2
    """
    
    f = C - a - b*x - c*x*x - d*x*x*x
    
    return f

def transform_der1(x,C,a,b,c,d):
    """First derivative of the general transformation function with respect
    to C2 == x, i.e.:
    
    df = -b - 2cC2 - 3dC2^2
    """
    
    f = -b -c*x -d*x*x
    
    return f
    
def transform_der2(x,C,a,b,c,d):
    """Second derivative of the general transformation function with respect
    to C2 == x, i.e.:
    
    df^2 = 2c - 3dC2
    """
    
    f = -c -d*x
    
    return f

def Jordi_coefficients(phot_system, A):
    """Function to return the coefficients published by Jordi for a specific
    photometric system.  Options are: SDSS
    """
    
    if phot_system == 'SDSS':
        coefficients = { 0: {'c': [0.3800, 0.8376, -0.0097, -0.0001], 'sig': 0.17},
                         1: {'c': [0.3482, 1.3463, -0.0310, 0.0067], 'sig': 0.35},
                         3: {'c': [0.4649, 2.3969, -0.2976, 0.0207], 'sig': 0.12},
                         5: {'c': [0.4052, 0.6407, -0.0091, 0.0004], 'sig': 0.11}}
    else:
        print('ERROR: Transformation to '+phot_system+' not yet supported')
        
    return coefficients[A]
    
def transform_SDSS_to_Gaia(params):
    """Photometric transformations, based on the tables of coefficients
    published by Jordi, C et al, (2010), A&A, 523, A48.
    
    Relationship of the form:
    C1 = a + bC2 + cC2^2 + dC2^3
    
    where C1 = colour involving at least one Gaia measurement and
    C2 = a colour in the Sloan system
    (Coefficients for the other systems are not yet implemented)
    
    The coefficients are only valid for discrete values of extinction at 
    wavelength=550nm, and interpolations between these values are non-trivial.
    """
    
    # Jordi coefficients are arranged for discrete intervals of extinction
    # at lambda=550nm (used as the dictionary key)"
    data = Jordi_coefficients('SDSS', params['A'])
    
    Gbprp = data['c'][0] + data['c'][1]*params['gi'] + \
                                data['c'][2]*params['gi']**2 + \
                                data['c'][3]*params['gi']**3
    
    sigGbprp = data['sig']
    
    return Gbprp, sigGbprp

def transform_Gaia_to_SDSS(params):
    """Photometric transformation from a Gaia colour measurement (bp-rp)
    to an SDSS colour, based on the tables of coefficients
    published by Jordi, C et al, (2010), A&A, 523, A48, and using
    the Halley method to invert the quoted transformation
    """
    
    coefficients = Jordi_coefficients('SDSS', params['A'])
    
    p = tuple([ params['Gbprp'] ] + coefficients['c'])
    
    gi = optimize.newton(transform, 0.0, fprime=transform_der1, 
                            fprime2=transform_der2, args=p)
    
    siggi = (params['sigGbprp']/params['Gbprp']) * gi
    
    return gi,siggi
    
def estimate_teff_from_Gbprp(Gbprp):
    """Relationship between stellar effective temperature and the Gaia 
    colour measurement CXP = Gbp - Grp from 
    Jordi, C et al, (2010), A&A, 523, A48.
    
    log(Teff) = 3.999−0.654(CXP)+0.709(CXP)^2 −0.316(CXP)^3

    The residual of the fit is quoted to be 0.02 dex, equivalent to a relative 
    error Delta(Teff)/Teff of ∼4.6%.
    """
    
    if Gbprp < 2.0:
        
        log_teff = 3.999 - (0.654*Gbprp) + (0.709*Gbprp*Gbprp) - \
                                (0.316 * Gbprp*Gbprp*Gbprp)
        
        teff = 10**(log_teff)
        
        sig_teff = 0.046 * teff
        
    else:
        
        teff = -99999.9999
        sig_teff = -99999.9999
        
        print('Warning: Gaia colour measurement outside valid range for calculation of Teff')
        
    return teff, sig_teff
    
if __name__ == '__main__':
    
    params = {}
    
    if len(argv) == 1:
        
        print('Useage:')
        print('> python jordi_phot_transforms.py A arg1=val arg2=val arg3=val...')
        print(' where A is the extinction band to use')
        print(' and args are of the form key=value, and optional keys are:')
        print(' gi, siggi  -> (g-i), uncertainty (g-i)')
        print(' Gbprp, sigBbprp  -> G(bp-rp), uncertainty G(bp-rp)')
        exit()
        
    else:
        
        params['A'] = float(argv[1])
        
        for arg in argv[2:]:
            
            if '=' in arg:
                
                (key, value) = arg.split('=')
                
                params[key] = float(value)
    
    if 'gi' in params.keys() and 'siggi' in params.keys():
        (Gbprp,sigGbprp) = transform_SDSS_to_Gaia(params)
        (teff,sig_teff) = estimate_teff_from_Gbprp(Gbprp)
        print('Gaia (bp-rp) = '+str(Gbprp)+'+/-'+str(sigGbprp)+' mag')
        print('Teff = '+str(teff)+'+/-'+str(sig_teff)+' K')
    
    if 'Gbprp' in params.keys():
        (gi, siggi) = transform_Gaia_to_SDSS(params)
        
        print('(g-i) = '+str(gi)+'+/-'+str(siggi)+' mag')
    