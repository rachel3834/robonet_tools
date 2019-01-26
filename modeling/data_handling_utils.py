# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 21:44:01 2019

@author: rstreet
"""

import numpy as np
import os

def apply_error_rescaling(rescaling_coeffs,name,lightcurve,use_error_rescaling):
    
    if len(rescaling) > 0:
        
        lightcurve[:,2] =  (rescaling[0]**2 + rescaling[1]**2*lightcurve[:,2]**2)**0.5
        
        print(' -> Applied error rescaling coefficient for dataset '+name+\
                ' a0='+str(rescaling[0])+', a1='+str(rescaling[1]))
    return lightcurve

def read_data_file(data_file,calibrated=False, gamma=None, rescaling=[]):
    """Function to decide which function to use to read the datafile"""
    
    if 'p.t' in data_file:
        
        tel = read_rbn_lightcurve(data_file, calibrated=calibrated, 
                                  gamma=gamma, rescaling=rescaling)
        
    elif 'cal.t' in data_file:
        
        tel = read_rbn_lightcurve(data_file, calibrated=calibrated, 
                                  gamma=gamma, rescaling=rescaling)
        
    elif 'DK-1.54' in data_file:
        
        tel = read_danish_lightcurve(data_file, gamma=gamma, rescaling=rescaling)
    
    return tel
    
def read_rbn_lightcurve(data_file, calibrated=False, gamma=None, rescaling=[]):
    """Function to read the standard DanDIA-format RoboNet lightcurve file.
    The calibrated option controls whether the (default) raw, instrumental 
    magnitudes are read or the calibrated magnitude columns
    Returns a pyLIMA telescope option.
    """
    
    if calibrated:
        mag_col = 8
        mag_err_col = 9
    else:
        mag_col = 6
        mag_err_col = 7
        
    lightcurve = np.loadtxt(data_file,dtype=str)
    lightcurve = np.c_[lightcurve[:,1],lightcurve[:,mag_col],
                       lightcurve[:,mag_err_col]].astype(float)

    print('\nRead RoboNet lightcurve '+os.path.basename(data_file))
    
    event_telescope = os.path.basename(data_file)
    
    aaa = event_telescope
    name = aaa[14:17]+aaa[21]+'_'+aaa[-4:-2]
    
    lightcurve = apply_error_rescaling(rescaling,name,lightcurve)

    telescope = telescopes.Telescope(name=name, camera_filter=aaa[-4:-2],
                     light_curve_magnitude=lightcurve,
                     light_curve_magnitude_dictionnary={'time': 0, 'mag': 1, 'err_mag': 2})

    if gamma != None:
        gamma = fetch_gamma(os.path.basename(event_telescope))
        telescope.gamma = gamma
        print(' -> Applied gamma='+str(gamma))

    return telescope
    
def read_danish_lightcurve(data_file, gamma=None, rescaling=[]):
    """Function to read the standard MiNDSTEp-format Danish lightcurve file,
    from the DanDIA pipeline adapted for EMCCDs.
    Returns a pyLIMA telescope option.
    """

    lightcurve = np.loadtxt(data_file,dtype=str)
    lightcurve = np.c_[lightcurve[:,1],lightcurve[:,6],lightcurve[:,7]].astype(float)
    
    print('\nRead Danish lightcurve '+os.path.basename(data_file))
    
    name = 'Danish'
    
    lightcurve = apply_error_rescaling(rescaling,name,lightcurve)
    
    telescope = telescopes.Telescope(name=name, camera_filter='Z',
                                        light_curve_magnitude=lightcurve,
                                        light_curve_magnitude_dictionnary={'time': 0, 'mag': 1, 'err_mag': 2})

    if gamma != None:
        gamma = fetch_gamma(os.path.basename(event_telescope))
        telescope.gamma = gamma
        print(' -> Applied gamma='+str(gamma))
    
    return telescope
    