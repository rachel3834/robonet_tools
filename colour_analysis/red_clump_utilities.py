# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 16:50:57 2018

@author: rstreet
"""

from sys import argv
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

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
    