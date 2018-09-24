# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 13:16:17 2018

@author: rstreet
"""

from sys import argv
import numpy as np

def transform_Gaia_to_SDSS(gi, siggi, A):
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
    coefficients = { 0: {'c': [0.3800, 0.8376, -0.0097, -0.0001], 'sig': 0.17},
                     1: {'c': [0.3482, 1.3463, -0.0310, 0.0067], 'sig': 0.35},
                     3: {'c': [0.4649, 2.3969, -0.2976, 0.0207], 'sig': 0.12},
                     5: {'c': [0.4052, 0.6407, -0.0091, 0.0004], 'sig': 0.11}}
    
    data = coefficients[A]
    
    Gbprp = data['c'][0] + data['c'][1]*gi + \
                                data['c'][2]*gi**2 + data['c'][3]*gi**3
    
    sigGbprp = data['sig']
    
    return Gbprp, sigGbprp

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
    
    if len(argv) == 1:
        
        gi = float(raw_input('Please enter the (g-i) colour: '))
        siggi = float(raw_input('Please enter the uncertainty on the (g-i) colour: '))
        A = float(raw_input('Please enter the extinction (550nm) band to use [0,1,3 or 5]:'))
        
    else:
        
        gi = float(argv[1])
        siggi = float(argv[2])
        A = float(argv[3])
        
    (Gbprp,sigGbprp) = transform_Gaia_to_SDSS(gi, siggi, A)
    (teff,sig_teff) = estimate_teff_from_Gbprp(Gbprp)
    
    print('Gaia (bp-rp) = '+str(Gbprp)+'+/-'+str(sigGbprp)+' mag')
    print('Teff = '+str(teff)+'+/-'+str(sig_teff)+' K')
    