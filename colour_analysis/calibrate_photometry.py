# -*- coding: utf-8 -*-
"""
Created on Wed May  9 16:33:56 2018

@author: rstreet
"""
import os
import sys
from pyDANDIA import pipeline_setup
from pyDANDIA import logs
from pyDANDIA import metadata
import vizier_tools

VERSION = 'calibrate_photometry_0.0.1'

def calibrate_photometry():
    """Function to calculate the photometric transform between the instrumental
    magnitudes produced by the pyDANDIA pipeline and catalog data."""
    
    params = get_args()

    setup = pipeline_setup.pipeline_setup(params)
    
    log = logs.start_stage_log( setup.red_dir, 'phot_calib', version=VERSION )
    
    reduction_metadata = metadata.MetaData()
    reduction_metadata.load_a_layer_from_file( setup.red_dir, 
                                              'pyDANDIA_metadata.fits', 
                                              'reduction_parameters' )
    reduction_metadata.load_a_layer_from_file( setup.red_dir, 
                                              'pyDANDIA_metadata.fits', 
                                              'star_catalog' )

    
def get_args():
    """Function to gather the necessary commandline arguments"""
    
    params = { 'red_dir': '',
               'log_dir': '',
               'catalog': '',
               'pipeline_config_dir': '',
               'software_dir': '', 
               'verbosity': '',
            }

    if len(sys.argv) > 1:
        params['red_dir'] = sys.argv[1]
        params['log_dir'] = sys.argv[2]
        params['field_ra'] = sys.argv[3]
        params['field_dec'] = sys.argv[3]
        params['radius'] = sys.argv[3]
    else:
        params['red_dir'] = raw_input('Please enter the path to the reduction directory: ')
        params['log_dir'] = raw_input('Please enter the path to the log directory: ')
        params['field_ra'] = raw_input('Please the field centre RA, J2000, in sexigesimal format: ')
        params['field_dec'] = raw_input('Please the field centre Dec, J2000, in sexigesimal format: ')
        params['radius'] = raw_input('Please the search radius in arcmin: ')
    
    return params


if __name__ == '__main__':
    
    calibrate_photometry()