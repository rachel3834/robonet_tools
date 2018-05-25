# -*- coding: utf-8 -*-
"""
Created on Tue May 15 20:48:12 2018

@author: rstreet
"""

from os import path
from sys import argv
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import matching
from astropy.table import Table
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

def plot_colour_diagrams():
    """Function to plot colour magnitude and colour-colour plots"""
        
    params = get_args()
    
    (star_catalog,catalog_header) = read_combined_star_catalog(params)
    
    target = find_target_data(params,star_catalog)

    print target
    
    deltas = calibrate_instrumental_colour_colour_diagram(params,star_catalog,catalog_header,target)
    
    plot_colour_mag_diagram(params,star_catalog,deltas,catalog_header,target)
    
    plot_colour_colour_diagram(params,star_catalog,deltas,catalog_header,target)
    
def get_args():
    """Function to gather the necessary commandline arguments"""

    params = {}
    
    if len(argv) == 1:
        
        params['catalog_file'] = raw_input('Please enter the path to combined star catalog file: ')
        params['red_dir'] = raw_input('Please enter the path to the output directory: ')
        params['target_ra'] = raw_input('Please enter the RA of target to highlight, if any [sexigesimal format or none]:')
        params['target_dec'] = raw_input('Please enter the Dec of target to highlight, if any [sexigesimal format or none]:')
        
    else:

        params['catalog_file'] = argv[1]
        params['red_dir'] = argv[2]
        params['target_ra'] = argv[3]
        params['target_dec'] = argv[4]
    
    for key, value in params.items():
        
        if 'none' in str(value).lower():
            
            value = None
        
            params[key] = value
            
    return params
    
def read_combined_star_catalog(params):
    """Function to read the photometric and star catalog data from a metadata file"""
    
    if path.isfile(params['catalog_file']) == False:
        
        return np.zeros(1)
    
    hdulist = fits.open(params['catalog_file'])
    
    data = hdulist[1].data
    
    header = hdulist[0].header
    
    star_catalog = Table(data)
    
    return star_catalog, header

def find_target_data(params,star_catalog):
    """Function to identify the photometry for a given target, if present
    in the star catalogue"""
    
    target = {}
    
    if params['target_ra'] != None:
        
        target_location = SkyCoord([params['target_ra']], [params['target_dec']], unit=(u.hourangle, u.deg))
                
        stars = SkyCoord(star_catalog['RA'], star_catalog['DEC'], unit="deg")
        
        tolerance = 2.0 * u.arcsec
        
        match_data = matching.search_around_sky(target_location, stars, 
                                                seplimit=tolerance)    
                                                
        idx = np.argsort(match_data[2].value)
    
        if len(match_data[0]) > 0:
            target['star_index'] = star_catalog['star_index'][match_data[1][idx[0]]]
            target['ra'] = star_catalog['RA'][match_data[1][idx[0]]]
            target['dec'] = star_catalog['DEC'][match_data[1][idx[0]]]
            target['ref_mag_ip'] = star_catalog['ref_mag_ip'][match_data[1][idx[0]]]
            target['ref_mag_rp'] = star_catalog['ref_mag_rp'][match_data[1][idx[0]]]
            target['separation'] = match_data[2][idx[0]].to_string(unit=u.arcsec)
            try:
                target['ref_mag_gp'] = star_catalog['ref_mag_gp'][match_data[1][idx[0]]]
            except AttributeError:
                pass
            print 'Target identified as '+repr(target)
    
    return target

def index_valid_star_entries(star_catalog,valid_cat=False):
    """Function to return an index of all stars with both full instrumental and
    catalogue entries"""
    
    idx1 = np.where(star_catalog['ref_mag_ip'] > 0.0)[0]
    idx2 = np.where(star_catalog['ref_mag_rp'] > 0.0)[0]
    idx3 = np.where(star_catalog['ref_mag_gp'] > 0.0)[0]
    idx4 = np.where(star_catalog['imag'] > 0.0)[0]
    idx5 = np.where(star_catalog['rmag'] > 0.0)[0]
    idx6 = np.where(star_catalog['gmag'] > 0.0)[0]
    idx = set(idx1).intersection(set(idx2))
    idx = idx.intersection(set(idx3))
    
    if valid_cat == False:
        return list(idx)
        
    idx = idx.intersection(set(idx4))
    idx = idx.intersection(set(idx5))
    idx = list(idx.intersection(set(idx6)))
    
    return idx
    
def plot_colour_mag_diagram(params,star_catalog,deltas,catalog_header,target):
    """Function to plot a colour-magnitude diagram"""
    
    filters = { 'ip': 'SDSS-i', 'rp': 'SDSS-r', 'gp': 'SDSS-g' }
    
    staridx = index_valid_star_entries(star_catalog,valid_cat=True)
    
    inst_i = star_catalog['ref_mag_ip'][staridx]
    inst_r = star_catalog['ref_mag_rp'][staridx]
    cal_i = star_catalog['imag'][staridx]
    cal_r = star_catalog['rmag'][staridx]
    inst_ri = inst_r - inst_i    # Catalogue column order is red -> blue
    cal_ri = cal_r - cal_i
    
    if len(target) > 0 and target['ref_mag_ip'] > 0.0 and \
            target['ref_mag_rp'] > 0.0:
        target_ri = target['ref_mag_rp'] - target['ref_mag_ip']
    
    fig = plt.figure(1)

    plt.scatter(inst_ri+deltas['dri'], 
                inst_r+deltas['dr'], 
                 c='#2b8c85', 
                 marker='.', s=1)
    
#    plt.scatter(cal_ri, cal_r, c='#8c6931', 
#             marker='.', s=1, 
#             label='Catalogue')
    
    if len(target) > 0 and target['ref_mag_ip'] > 0.0 and \
            target['ref_mag_rp'] > 0.0:
        plt.plot(target_ri+deltas['dri'], 
                 target['ref_mag_rp']+deltas['dr'],
                 'md',markersize=6)
        
    plt.xlabel('SDSS (r-i) [mag]')

    plt.ylabel('SDSS-r [mag]')
    
    [xmin,xmax,ymin,ymax] = plt.axis()
    
    plt.axis([xmin,xmax,ymax,ymin])
    
    plot_file = path.join(params['red_dir'],'colour_magnitude_diagram.png')

    plt.grid()
    
#    plt.legend()
    
    plt.axis([0.0,2.0,20.0,14.0])
    
    plt.savefig(plot_file)

    plt.close(1)
    
    print 'Colour-magnitude diagram output to '+plot_file
    
def plot_colour_colour_diagram(params,star_catalog,deltas,catalog_header,target):
    """Function to plot a colour-colour diagram, if sufficient data are
    available within the given star catalog"""
    
    filters = { 'ip': 'SDSS-i', 'rp': 'SDSS-r', 'gp': 'SDSS-g' }
    
    try:
    
        staridx = index_valid_star_entries(star_catalog,valid_cat=True)
        
        inst_i = star_catalog['ref_mag_ip'][staridx]
        inst_r = star_catalog['ref_mag_rp'][staridx]
        inst_g = star_catalog['ref_mag_gp'][staridx]
        inst_gr = inst_g - inst_r    # Catalogue column order is red -> blue
        inst_ri = inst_r - inst_i   
        
        if len(target) > 0 and target['ref_mag_ip'] > 0.0 and \
            target['ref_mag_rp'] > 0.0 and target['ref_mag_gp']:
            target_gr = target['ref_mag_gp'] - target['ref_mag_rp']
            target_ri = target['ref_mag_rp'] - target['ref_mag_ip']
    
        fig = plt.figure(1)
    
        plt.scatter(inst_gr+deltas['dgr'], 
                    inst_ri+deltas['dri'], 
                     c='#2b8c85', 
                     marker='.', s=1)
        
        if len(target) > 0 and target['ref_mag_ip'] > 0.0 and \
            target['ref_mag_rp'] > 0.0 and target['ref_mag_gp']:
            plt.plot(target_gr+deltas['dgr'], 
                     target_ri+deltas['dri'],
                     'md',markersize=6)
            
        plt.xlabel('SDSS (g-r) [mag]')
    
        plt.ylabel('SDSS (r-i) [mag]')
        
        plot_file = path.join(params['red_dir'],'colour_colour_diagram.png')
        
        plt.axis([0.5,3.0,0.0,2.0])
    
        plt.grid()
        
        plt.savefig(plot_file)
    
        plt.close(1)
        
        print 'Colour-colour diagram output to '+plot_file
        
    except AttributeError:
            
        print 'Warning: Insufficient data for colour-colour diagram'
        
def calibrate_instrumental_colour_colour_diagram(params,star_catalog,catalog_header,target):
    """Function to plot a colour-colour diagram, if sufficient data are
    available within the given star catalog"""
    
    filters = { 'ip': 'SDSS-i', 'rp': 'SDSS-r', 'gp': 'SDSS-g' }
    
    try:
    
        staridx = index_valid_star_entries(star_catalog,valid_cat=True)
        
        inst_i = star_catalog['ref_mag_ip'][staridx]
        inst_r = star_catalog['ref_mag_rp'][staridx]
        inst_g = star_catalog['ref_mag_gp'][staridx]
        inst_gr = inst_g - inst_r    # Catalogue column order is red -> blue
        inst_ri = inst_r - inst_i   
        
        cal_i = star_catalog['imag'][staridx]
        cal_r = star_catalog['rmag'][staridx]
        cal_g = star_catalog['gmag'][staridx]
        cal_gr = cal_g - cal_r    # Catalogue column order is red -> blue
        cal_ri = cal_r - cal_i    
        
        deltas = {}
        deltas['di'] = np.median(cal_i - inst_i)
        deltas['dr'] = np.median(cal_r - inst_r)
        deltas['dg'] = np.median(cal_g - inst_g)
        deltas['dgr'] = np.median(cal_gr - inst_gr)
        deltas['dri'] = np.median(cal_ri - inst_ri)
        
        if len(target) > 0 and target['ref_mag_ip'] > 0.0 and \
            target['ref_mag_rp'] > 0.0 and target['ref_mag_gp']:
            target_gr = target['ref_mag_gp'] - target['ref_mag_rp']
            target_ri = target['ref_mag_rp'] - target['ref_mag_ip']
    
        fig = plt.figure(1)
    
        plt.plot(inst_gr, inst_ri,'.',
                 color='#2b8c85',markersize=1,alpha=0.4,
                 label='Instrumental')
        plt.plot(cal_gr, cal_ri,'.',
                 color='#8c6931',markersize=1,alpha=0.4,
                 label='Catalogue')
        
        plt.plot(inst_gr+deltas['dgr'], inst_ri+deltas['dri'],'.',
                 color='m',markersize=1,alpha=0.4,
                 label='Calibrated')
                 
#        if len(target) > 0 and target['mag1'] > 0.0 and \
#            target['mag2'] > 0.0 and target['mag2']:
#            plt.plot(target_colour2, target_colour1,'md',markersize=6)
            
        plt.xlabel('SDSS (g-r) [mag]')
    
        plt.ylabel('SDSS (r-i) [mag]')
        
        plot_file = path.join(params['red_dir'],'colour_colour_diagram.png')
        
        plt.grid()
        
        plt.legend()
        
        plt.savefig(plot_file)
    
        plt.close(1)
        
        print 'Colour-colour diagram output to '+plot_file
        
    except AttributeError:
        
        deltas = {}
        print 'Warning: Insufficient data for colour-colour diagram'
    
    return deltas
    
if __name__ == '__main__':
    
    plot_colour_diagrams()
    