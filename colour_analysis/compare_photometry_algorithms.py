# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 10:58:53 2019

@author: rstreet
"""

from os import path
from sys import argv
from pyDANDIA import logs, metadata
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt

def compare_pydandia_dandia_ref_analysis():
    """Function to compare the reference frame photometry from DanDIA and
    pyDANDIA, assuming both pipelines have been used to reduce the same 
    reference image.  The DanDIA data may be cropped but the pyDANDIA analysis
    is assumed to take place over the fullframe image."""
    
    params = get_args()
    
    star_catalog = load_pydandia_star_catalog(params)
    
    starlist = load_dandia_starlist(params)
    
    starlist = offset_starlist_coords(starlist, params)

    (pydandia_positions, dandia_positions) = extract_pixel_positions(star_catalog, starlist)
    
    match_index = match_on_pixel_position(pydandia_positions, dandia_positions, 
                                          2.0, verbose=False)
    
    compare_photometry(star_catalog, starlist, match_index)
    
def get_args():
    
    params = {}
    if len(argv) == 1:
        params['metadata'] = raw_input('Please enter the path to the pyDANDIA metadata file: ')
        params['starlist'] = raw_input('Please enter the path to the DanDIA stackref.starlist file: ')
        params['cropx'] = float(raw_input('Please enter the minimum x-coordinate of the DanDIA crop frame in fullframe pixel coordinates: '))
        params['cropy'] = float(raw_input('Please enter the minimum y-coordinate of the DanDIA crop frame in fullframe pixel coordinates: '))
    else:
        params['metadata'] = argv[1]
        params['starlist'] = argv[2]
        params['cropx'] = float(argv[3])
        params['cropy'] = float(argv[4])
    
    return params
    
def load_pydandia_star_catalog(params):
    """Function to read the photometric and star catalog data from a metadata file"""
    
    meta_file = params['metadata']
    
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

    print('Extracted data for '+str(len(star_catalog))+\
            ' stars from the pyDANDIA analysis')
            
    return star_catalog

def load_dandia_starlist(params):
    
    starlist = np.loadtxt(params['starlist'],dtype=str)
    starlist = starlist.astype(float)
    
    print('Extracted data for '+str(len(starlist))+\
            ' stars from the DanDIA analysis')
            
    return starlist

def offset_starlist_coords(starlist, params):
    
    starlist[:,1] = starlist[:,1] + params['cropx']
    starlist[:,3] = starlist[:,3] + params['cropy']
    
    print('Offset DanDIA catalog coordinates by '+str(params['cropx'])+\
                                             ', '+str(params['cropy']))
    
    return starlist

def extract_pixel_positions(star_catalog, starlist):
    
    pydandia_positions = np.zeros( (len(star_catalog),2) )
    pydandia_positions[:,0] = star_catalog['x'].data
    pydandia_positions[:,1] = star_catalog['y'].data
    
    dandia_positions = np.zeros( (len(starlist),2) )
    dandia_positions[:,0] = starlist[:,1]
    dandia_positions[:,1] = starlist[:,3]
    
    return pydandia_positions, dandia_positions
    
def match_on_pixel_position(catalog1, catalog2, tolerance, verbose=False):
    """Function to find matches between two arrays of (x,y) pixel positions.
    Catalog 2 is assumed to cover a subset of catalog1.
    Returned:
        match_index int array [ index_catalog1, index_catalog2 ]
    """
    
    print('Matching stars by pixel positions')
    
    match_index = []
    
    for j in range(0,len(catalog2),1):
        
        dx = catalog1[:,0] - catalog2[j,0]
        dy = catalog1[:,1] - catalog2[j,1]
        
        sep = np.sqrt(dx*dx + dy*dy)
        
        idx = sep.argsort()
        
        if sep[idx[0]] <= tolerance:
            
            match_index.append( [idx[0], j] )
            
            if verbose:
                print(match_index[-1], catalog1[idx[0],0], catalog1[idx[0],1], \
                    catalog2[j,0], catalog2[j,1], sep[idx[0]])
    
    print('-> Matched '+str(len(match_index))+' stars')
    
    return np.array(match_index)

def compare_photometry(star_catalog, starlist, match_index):
    
    (dandia_mag, dandia_mag_err) = convert_flux_to_mag(starlist[match_index[:,1],5],
                                                        starlist[match_index[:,1],6])

    pydandia_mag = star_catalog['mag'].data[match_index[:,0]]
    pydandia_mag_err = star_catalog['mag_err'].data[match_index[:,0]]
    
    idx1 = np.where(dandia_mag > 0.0)
    idx2 = np.where(pydandia_mag > 0.0)
    idx = list(set(idx1[0]).intersection(set(idx2[0])))
    
    fig = plt.figure(1,(10,10))
    
    plt.plot(dandia_mag[idx], pydandia_mag[idx], 'm.', markersize=1)
    
    plt.xlabel('DanDIA magnitude')
    plt.ylabel('pyDANDIA magnitude')
    
    plt.savefig('compare_ref_magnitudes_algo.png')
    
    plt.close(1)
    
def convert_flux_to_mag(flux, flux_err):
    """Function to convert the flux of a star from its fitted PSF model 
    and its uncertainty onto the magnitude scale.
    
    :param float flux: Total star flux
    :param float flux_err: Uncertainty in star flux
    
    Returns:
    
    :param float mag: Measured star magnitude
    :param float mag_err: Uncertainty in measured magnitude
    """
    
    ZP = 25.0

    mag = np.zeros(len(flux))
    mag_err = np.zeros(len(flux))
    
    idx1 = np.where(flux > 0.0)
    idx2 = np.where(flux_err > 0.0)
    idx = list(set(idx1[0]).intersection(set(idx2[0])))
    
    mag[idx] = ZP - 2.5 * np.log10(flux[idx])

    mag_err[idx] = (2.5 / np.log(10.0)) * flux_err[idx] / flux[idx]

    return mag, mag_err
    
if __name__ == '__main__':

    compare_pydandia_dandia_ref_analysis()