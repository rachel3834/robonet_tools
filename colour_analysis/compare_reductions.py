# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 14:10:45 2019

@author: rstreet
"""

from os import path
from sys import argv
from pyDANDIA import logs, metadata
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import matching
import numpy as np
import matplotlib.pyplot as plt
    
def compare_star_catalogs():
    """Function to compare the photometry from two different star catalogs"""
    
    params = get_args()
    
    log = logs.start_stage_log( '.', 'combine_colour_datasets' )
    
    star_cat1 = extract_star_catalog(params,'cat1',log)
    star_cat2 = extract_star_catalog(params,'cat2',log)

    match_index = cross_match_stars(star_cat1,star_cat2,log)
    
    plot_matched_ref_photometry(star_cat1,star_cat2,match_index,log)
    
    logs.close_log(log)
    
def get_args():
    
    params = {}
    if len(argv) == 1:
        params['cat1'] = raw_input('Please enter the path to the first metadata file: ')
        params['cat2'] = raw_input('Please enter the path to the second metadata file: ')
    else:
        params['cat1'] = argv[1]
        params['cat2'] = argv[2]
    
    return params
    
def extract_star_catalog(params,filter_id,log):
    """Function to read the photometric and star catalog data from a metadata file"""
    
    meta_file = params[filter_id]
    
    reduction_metadata = metadata.MetaData()
    reduction_metadata.load_a_layer_from_file( path.dirname(meta_file), path.basename(meta_file), 'star_catalog' )
    
    star_catalog = Table()
    star_catalog['star_index'] = reduction_metadata.star_catalog[1]['star_index']
    star_catalog['x'] = reduction_metadata.star_catalog[1]['x_pixel']
    star_catalog['y'] = reduction_metadata.star_catalog[1]['y_pixel']
    star_catalog['RA'] = reduction_metadata.star_catalog[1]['RA_J2000']
    star_catalog['DEC'] = reduction_metadata.star_catalog[1]['DEC_J2000']
    star_catalog['mag'] = reduction_metadata.star_catalog[1]['ref_mag']
    star_catalog['mag_err'] = reduction_metadata.star_catalog[1]['ref_mag_err']
    star_catalog['flux'] = reduction_metadata.star_catalog[1]['ref_flux']
    star_catalog['flux_err'] = reduction_metadata.star_catalog[1]['ref_flux_err']

    log.info('Extracted data for '+str(len(star_catalog))+\
            ' stars for dataset '+filter_id)
            
    return star_catalog

def cross_match_stars(star_cat1,star_cat2,log):
    """Function to cross-match stars by position.
    
    Returns:
        :param dict match_index: { Index in vphas_cat: index in star_cat }
    """
    
    stars1 = SkyCoord(star_cat1['RA'], star_cat1['DEC'], unit="deg")
    stars2 = SkyCoord(star_cat2['RA'], star_cat2['DEC'], unit="deg")
    
    tolerance = 1.0 * u.arcsec
    
    match_data = matching.search_around_sky(stars1, stars2, 
                                             seplimit=tolerance)   
                                             
    idx = np.argsort(match_data[2].value)
    
    match_index = np.array(zip(match_data[0][idx],match_data[1][idx]))
    
    log.info('Matched '+str(len(match_index))+' stars between star catalogs')
        
    return match_index

def plot_matched_ref_photometry(star_cat1,star_cat2,match_index,log):
    """Function to compare the reference frame photometry from the two 
    catalogues"""
    
    mag1 = star_cat1['mag']
    magerr1 = star_cat1['mag_err']
    mag2 = star_cat2['mag']
    magerr2 = star_cat2['mag_err']

    idx1 = np.where(magerr1 > 0.0)[0]
    idx2 = np.where(magerr2 > 0.0)[0]
    idx = list(set(idx1).intersection(set(idx2)))
    
    fig = plt.figure(1,(10,10))
    
    plt.plot(star_cat1['mag'][idx], star_cat1['mag_err'][idx], 
                 'm.',markersize=1, label='Catalogue 1')
    plt.plot(star_cat2['mag'][idx], star_cat2['mag_err'][idx], 
                 'g.',markersize=1,label='Catalogue 2')
    
    (xmin,xmax,ymin,ymax) = plt.axis()
    plt.axis([xmin,xmax,0.0,ymax])
    
    plt.xlabel('Magnitude')
    plt.ylabel('Magnitude error')
    plt.legend()
    
    plt.savefig('compare_ref_photometry.png')
    
    plt.close(1)
    
    fig = plt.figure(2,(10,10))
    
    plt.plot(star_cat1['mag'][idx], star_cat2['mag'][idx], 'm.', markersize=1)
    
    plt.xlabel('Magnitude catalog 1')
    plt.ylabel('Magnitude catalog 2')
    
    plt.savefig('compare_ref_magnitudes.png')
    
    plt.close(2)
    
if __name__ == '__main__':

    compare_star_catalogs()