# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 16:50:57 2018

@author: rstreet
"""

from sys import argv
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

def get_essential_parameters(data=None):
    """Function to provide the published values for the essential parameters
    of Red Clump giant stars.
    Sources:
    Bensby et al. (2017), 2017, A&A, 605A, 89 for V, I bands and
    Ruiz-Dern et al. (2017), A&A, 609A, 116 for Gaia, 2MASS, SDSS and Wise bands.
    """
    
    if data == None:
        data = {}
        
    data['M_I_0'] = -0.12
    data['V-I_0'] = 1.09
    data['M_V_0'] = 0.97
    
    data['M_G_0'] = 0.495
    data['sigMG_0'] = 0.009
    data['M_J_0'] = -0.945
    data['sigMJ_0'] = 0.01
    data['M_H_0'] = -1.450
    data['sigMH_0'] = 0.017
    data['M_Ks_0'] = -1.606
    data['sigMKs_0'] = 0.009
    data['M_g_0'] = 1.331
    data['sigMg_0'] = 0.056
    data['M_r_0'] = 0.552
    data['sigMr_0'] = 0.026
    data['M_i_0'] = 0.262
    data['sigMi_0'] = 0.032
    data['M_W1_0'] = -1.711
    data['sigMW1_0'] = 0.017
    data['M_W2_0'] = -1.585
    data['sigMW2_0'] = 0.016
    data['M_W3_0'] = -1.638
    data['sigMW3_0'] = 0.011
    data['M_W4_0'] = -1.704
    data['sigMW4_0'] = 0.012

    data['J-H_0'] = data['M_J_0'] - data['M_H_0']
    data['sigJH_0'] = np.sqrt( (data['sigMJ_0']*data['sigMJ_0']) + \
                                (data['sigMH_0']*data['sigMH_0']) )
    data['H-K_0'] = data['M_H_0'] - data['M_Ks_0']
    data['sigHK_0'] = np.sqrt( (data['sigMH_0']*data['sigMH_0']) + \
                                (data['sigMKs_0']*data['sigMKs_0']) )
    
    data['g-r_0'] = data['M_g_0'] - data['M_r_0']
    data['siggr_0'] = np.sqrt( (data['sigMg_0']*data['sigMg_0']) + \
                                (data['sigMr_0']*data['sigMr_0']) )
    data['r-i_0'] = data['M_r_0'] - data['M_i_0']
    data['sigri_0'] = np.sqrt( (data['sigMr_0']*data['sigMr_0']) + \
                                (data['sigMi_0']*data['sigMi_0']) )
    
    return data
    
def calc_red_clump_distance(ra,dec):
    """Function to estimate the distance of the Red Clump stars in the
    Galactic Bulge from the observer on Earth, taking into account the 
    bar structure, using the relations from 
    Nataf, D. M., Gould, A., Fouqué, P., et al. 2013, ApJ, 769, 88
    """
    
    c = SkyCoord(ra, dec, unit=(u.hourangle,u.degree), frame='icrs')
    
    print('Galactic coordinates for '+ra+', '+dec+' (l,b) [deg]: '+\
            str(c.galactic.l.deg)+', '+str(c.galactic.b.deg))
    
    R_0 = 8.16 # Kpc
    phi = 40.0 * (np.pi/180.0)
    
    D_RC = R_0 / (np.cos(c.galactic.l.radian) + np.sin(c.galactic.l.radian)*(1.0/np.tan(phi)))
    
    print('Red Clump distance for ('+str(c.galactic.l.deg)+', '+\
                                    str(c.galactic.b.deg)+') = '+\
                                    str(D_RC)+' Kpc')
    
    return D_RC

def calc_I_apparent(D_RC):
    """Function to return the apparent magnitude in I-band of the Red Clump
    stars given their distance from the observer, D_RC, and the measured 
    apparent I magnitude of Red Clump stars in the Galactic Center, 
    published by 
    Nataf, D. M., Gould, A., Fouqué, P., et al. 2013, ApJ, 769, 88
    """

    R_0 = 8.16 # Kpc
    
    delta_I = 5.0 * np.log10(R_0/D_RC)
    
    I_RC_app = 14.443 + delta_I
    
    print('\nOffset in delta_I magnitude of Red Clump apparent magnitude, due to separation from the Galactic centre: '+str(delta_I)+'mag')
    print('Apparent I_RC at this distnace, I_RC_app = '+str(I_RC_app)+'mag')
    
    return I_RC_app

def calc_apparent_magnitudes(data):
    """Function to calculate the apparent magnitudes of Red Clump stars
    at a given distance from the observer, based on their absolute magnitudes
    """
    
    # Absolute magnitudes are for a distance of 10pc by definition:
    R_0 = 10.0 # pc
    
    delta_m = 5.0 * np.log10(data['D_RC']*1000.0) - 5.0
    
    passbands = [ 'G', 'J', 'H', 'Ks', 'g', 'r', 'i', 'W1', 'W2', 'W3', 'W4' ]
    
    for f in passbands:
        
        abs_mag = 'M_'+f+'_0'
        abs_mag_err = 'sigM'+f+'_0'
        app_mag = 'm_'+f+'_0'
        app_mag_err = 'sigm'+f+'_0'
        
        if abs_mag in data.keys():
            
            data[app_mag] = data[abs_mag] + delta_m
            data[app_mag_err] = data[abs_mag_err]
            
    return data
    
def calc_distance(params):
    '''Function to calculate the distance to a star from its apparent and
    absolute magnitudes, using the standard distance modulus expression:
    D = 10^( m - M + 5 )/5
    '''
    
    def calc_D(m,M):
        log_D = (m - M + 5.0 ) / 5.0
        D = 10**( log_D )
        return log_D,D
	
    (log_D,D) = calc_D( params['app_mag'], params['abs_mag'] )
    (log_D_min,D_min) = calc_D( params['app_mag']-params['sig_app_mag'], \
                             params['abs_mag']-params['sig_abs_mag'] )
    (log_D_max,D_max) = calc_D( params['app_mag']+params['sig_app_mag'], \
                             params['abs_mag']+params['sig_abs_mag'] )
    sig_D = ( D_max - D_min ) / 2.0
    
    print 'Distance = '+str(D)+'+/-'+str(sig_D)+'pc'


if __name__ == '__main__':
    
    if len(argv) == 1:
        ra = raw_input('Please enter the RA (sexigesimal, J2000.0): ')
        dec = raw_input('Please enter the Dec (sexigesimal, J2000.0): ')
    else:
        ra = argv[1]
        dec = argv[2]
    
    D_RC = calc_red_clump_distance(ra,dec)
    