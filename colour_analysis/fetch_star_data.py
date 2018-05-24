# -*- coding: utf-8 -*-
"""
Created on Wed May 23 16:59:12 2018

@author: rstreet
"""

import sys
from os import path
from pyDANDIA import metadata
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.coordinates import matching
import astropy.units as u
import numpy as np

def fetch_star_on_position():
    """Function to extract the data for a given star from a star catalogue
    based on its position"""
    
    params = get_args()
    
    star_catalog = select_data_from_star_catalog(params)

    target = find_target_data(params,star_catalog)
    
def find_target_data(params,star_catalog):
    """Function to identify the photometry for a given target, if present
    in the star catalogue"""
    
    target = {}
    
    if params['target_ra'] != None:
        
        target_location = SkyCoord([params['target_ra']], [params['target_dec']], unit=(u.hourangle, u.deg))
                
        stars = SkyCoord(star_catalog['RA'], star_catalog['DEC'], unit="deg")
        
        tolerance = 5.0 * u.arcsec
        
        match_data = matching.search_around_sky(target_location, stars, 
                                                seplimit=tolerance)            
        if len(match_data[0]) > 0:
            
            idx = np.argsort(match_data[2].value)
            
            print str(len(match_data[1]))+' matching objects in the catalog'
            
            target = {'star_index': star_catalog['star_index'][match_data[1][idx[0]]],
                      'x': star_catalog['x_pixel'][match_data[1][idx[0]]],
                      'y': star_catalog['y_pixel'][match_data[1][idx[0]]],
                      'ra': star_catalog['RA'][match_data[1][idx[0]]],
                      'dec': star_catalog['DEC'][match_data[1][idx[0]]],
                      'mag': star_catalog['mag'][match_data[1][idx[0]]],
                      'mag_err': star_catalog['mag_err'][match_data[1][idx[0]]]}
                      
            print 'Target identified as '+repr(target)
            
            print 'List of matching objects: '
            
            for i in range(len(match_data[1])):
                print str(star_catalog['star_index'][match_data[1][idx[i]]])+' '+\
                      str(star_catalog['x_pixel'][match_data[1][idx[i]]])+' '+\
                      str(star_catalog['y_pixel'][match_data[1][idx[i]]])+' '+\
                      str(star_catalog['RA'][match_data[1][idx[i]]])+' '+\
                      str(star_catalog['DEC'][match_data[1][idx[i]]]),match_data[2][idx[i]]
            
    return target

    
def select_data_from_star_catalog(params):
    """Function to read the photometric and star catalog data from a metadata file"""
    
    meta_file = params['metadata']
    
    reduction_metadata = metadata.MetaData()
    reduction_metadata.load_a_layer_from_file( path.dirname(meta_file), path.basename(meta_file), 'star_catalog' )

    star_catalog = Table()
    star_catalog['star_index'] = reduction_metadata.star_catalog[1]['star_index']
    star_catalog['x_pixel'] = reduction_metadata.star_catalog[1]['x_pixel']
    star_catalog['y_pixel'] = reduction_metadata.star_catalog[1]['y_pixel']
    star_catalog['RA'] = reduction_metadata.star_catalog[1]['RA_J2000']
    star_catalog['DEC'] = reduction_metadata.star_catalog[1]['DEC_J2000']
    star_catalog['mag'] = reduction_metadata.star_catalog[1]['ref_mag']
    star_catalog['mag_err'] = reduction_metadata.star_catalog[1]['ref_mag_err']

    return star_catalog

def get_args():
    """Function to gather the required arguments"""
    
    params = {}
    
    if len(sys.argv) == 1:
        
        params['metadata'] = raw_input('Please enter the path to the metadata file: ')
        params['target_ra'] = raw_input('Please the target RA in sexigesimal format: ')
        params['target_dec'] = raw_input('Please enter the target Dec in sexigesimal format: ')
        
    else:

        params['metadata'] = sys.argv[1]
        params['target_ra'] = sys.argv[2]
        params['target_dec'] = sys.argv[3]
    
    params['red_dir'] = path.dirname(params['metadata'])
    
    return params

if __name__ == '__main__':
    
    fetch_star_on_position()