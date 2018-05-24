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
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

def plot_colour_diagrams():
    """Function to plot colour magnitude and colour-colour plots"""
        
    params = get_args()
    
    (star_catalog,catalog_header) = read_combined_star_catalog(params)
    
    target = find_target_data(params,star_catalog)
    
    print target
    
    #plot_colour_mag_diagram(params,star_catalog,catalog_header,target)
    
    #plot_colour_colour_diagram(params,star_catalog,catalog_header,target)

    plot_instrumental_colour_colour_diagram(params,star_catalog,catalog_header,target)
    
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
    
    star_catalog = []

    for i in range(0,len(data),1):
        
        star_catalog.append( list( data[i] ) )
    
    star_catalog = np.array( star_catalog )
        
    return star_catalog, header

def find_target_data(params,star_catalog):
    """Function to identify the photometry for a given target, if present
    in the star catalogue"""
    
    target = {}
    
    if params['target_ra'] != None:
        
        target_location = SkyCoord([params['target_ra']], [params['target_dec']], unit=(u.hourangle, u.deg))
                
        stars = SkyCoord(star_catalog[:,1], star_catalog[:,2], unit="deg")
        
        tolerance = 2.0 * u.arcsec
        
        match_data = matching.search_around_sky(target_location, stars, 
                                                seplimit=tolerance)    
        print match_data
        
        if len(match_data[0]) > 0:
            target['star_index'] = star_catalog[match_data[1][0],0]
            target['ra'] = star_catalog[match_data[1][0],1]
            target['dec'] = star_catalog[match_data[1][0],2]
            target['mag1'] = star_catalog[match_data[1][0],3]
            target['cal_mag1'] = star_catalog[match_data[1][0],5]
            target['mag2'] = star_catalog[match_data[1][0],7]
            target['cal_mag2'] = star_catalog[match_data[1][0],9]
            target['separation'] = match_data[2][0].to_string(unit=u.arcsec)
            if star_catalog.shape[1] == 15:
                target['mag3'] = star_catalog[match_data[1][0],11]
                target['cal_mag3'] = star_catalog[match_data[1][0],13]
            print 'Target identified as '+repr(target)
    
    return target
    
def plot_colour_mag_diagram(params,star_catalog,catalog_header,target):
    """Function to plot a colour-magnitude diagram"""
    
    filters = { 'ip': 'SDSS-i', 'rp': 'SDSS-r', 'gp': 'SDSS-g' }
    
    idx1 = np.where(star_catalog[:,3] > 10.0)[0]
    idx2 = np.where(star_catalog[:,5] > 10.0)[0]
    idx = list(set(idx1).intersection(set(idx2)))
    
    mag1 = star_catalog[idx,3]
    mag2 = star_catalog[idx,5]
    colour = mag2 - mag1    # Catalogue column order is red -> blue
    
    if len(target) > 0 and target['mag1'] > 0.0 and \
            target['mag2'] > 0.0:
        target_colour = target['mag2'] - target['mag1']
    
    fig = plt.figure(1)

    plt.plot(colour, mag2,'k.',markersize=1)
    
    if len(target) > 0 and target['mag1'] > 0.0 and \
            target['mag2'] > 0.0:
        plt.plot(target_colour, target['mag2'],'md',markersize=6)
        
    plt.xlabel(filters[catalog_header['FILTER2']]+'-'+\
                filters[catalog_header['FILTER1']]+' [mag]')

    plt.ylabel(filters[catalog_header['FILTER2']]+' [mag]')
    
    [xmin,xmax,ymin,ymax] = plt.axis()
    
    plt.axis([xmin,xmax,ymax,ymin])
    
    plot_file = path.join(params['red_dir'],'colour_magnitude_diagram.png')
    
    plt.savefig(plot_file)

    plt.close(1)
    
    print 'Colour-magnitude diagram output to '+plot_file
    
def plot_colour_colour_diagram(params,star_catalog,catalog_header,target):
    """Function to plot a colour-colour diagram, if sufficient data are
    available within the given star catalog"""
    
    filters = { 'ip': 'SDSS-i', 'rp': 'SDSS-r', 'gp': 'SDSS-g' }
    
    if star_catalog.shape[1] == 9:
    
        idx1 = np.where(star_catalog[:,3] > 0.0)[0]
        idx2 = np.where(star_catalog[:,5] > 0.0)[0]
        idx3 = np.where(star_catalog[:,7] > 0.0)[0]
        idx = set(idx1).intersection(set(idx2))
        idx = list(idx.intersection(set(idx3)))
        
        mag1 = star_catalog[idx,3]
        mag2 = star_catalog[idx,5]
        mag3 = star_catalog[idx,7]
        colour1 = mag2 - mag1    # Catalogue column order is red -> blue
        colour2 = mag3 - mag2   
        
        if len(target) > 0 and target['mag1'] > 0.0 and \
            target['mag2'] > 0.0 and target['mag2']:
            target_colour1 = target['mag2'] - target['mag1']
            target_colour2 = target['mag3'] - target['mag2']
    
        fig = plt.figure(1)
    
        plt.plot(colour2, colour1,'k.',markersize=1)
        
        if len(target) > 0 and target['mag1'] > 0.0 and \
            target['mag2'] > 0.0 and target['mag2']:
            plt.plot(target_colour2, target_colour1,'md',markersize=6)
            
        plt.xlabel(filters[catalog_header['FILTER3']]+'-'+\
                    filters[catalog_header['FILTER2']]+' [mag]')
    
        plt.ylabel(filters[catalog_header['FILTER2']]+'-'+\
                    filters[catalog_header['FILTER1']]+' [mag]')
        
        plot_file = path.join(params['red_dir'],'colour_colour_diagram.png')
        
        plt.savefig(plot_file)
    
        plt.close(1)
        
        print 'Colour-colour diagram output to '+plot_file
        
    else:
            
        print 'Warning: Insufficient data for colour-colour diagram'
        
def plot_instrumental_colour_colour_diagram(params,star_catalog,catalog_header,target):
    """Function to plot a colour-colour diagram, if sufficient data are
    available within the given star catalog"""
    
    filters = { 'ip': 'SDSS-i', 'rp': 'SDSS-r', 'gp': 'SDSS-g' }
    
    if star_catalog.shape[1] == 15:
    
        idx1 = np.where(star_catalog[:,5] > 0.0)[0]
        idx2 = np.where(star_catalog[:,9] > 0.0)[0]
        idx3 = np.where(star_catalog[:,13] > 0.0)[0]
        idx = set(idx1).intersection(set(idx2))
        idx = list(idx.intersection(set(idx3)))
        
        inst_mag1 = star_catalog[idx,3]
        inst_mag2 = star_catalog[idx,7]
        inst_mag3 = star_catalog[idx,11]
        inst_colour1 = inst_mag3 - inst_mag2    # Catalogue column order is red -> blue
        inst_colour2 = inst_mag2 - inst_mag1    
        
        jdx = np.where(inst_colour1 > 10.0)   
        print jdx
        print star_catalog[idx[jdx[0][0]],:]
        
        cal_mag1 = star_catalog[idx,5]
        cal_mag2 = star_catalog[idx,9]
        cal_mag3 = star_catalog[idx,13]
        cal_colour1 = cal_mag3 - cal_mag2    # Catalogue column order is red -> blue
        cal_colour2 = cal_mag2 - cal_mag1    
        
        if len(target) > 0 and target['mag1'] > 0.0 and \
            target['mag2'] > 0.0 and target['mag2']:
            target_colour1 = target['mag3'] - target['mag2']
            target_colour2 = target['mag2'] - target['mag1']
    
        fig = plt.figure(1)
    
        plt.plot(inst_colour1, inst_colour2,'.',
                 color='#8c6931',markersize=1,alpha=0.4,
                 label='Instrumental')
#        plt.plot(cal_colour1, cal_colour2,'.',
#                 color='#2b8c85',markersize=1,alpha=0.4,
#                 label='Calibrated')
        
#        if len(target) > 0 and target['mag1'] > 0.0 and \
#            target['mag2'] > 0.0 and target['mag2']:
#            plt.plot(target_colour2, target_colour1,'md',markersize=6)
            
        plt.xlabel(filters[catalog_header['FILTER3']]+'-'+\
                    filters[catalog_header['FILTER1']]+' [mag]')
    
        plt.ylabel(filters[catalog_header['FILTER2']]+'-'+\
                    filters[catalog_header['FILTER1']]+' [mag]')
        
        plot_file = path.join(params['red_dir'],'colour_colour_diagram.png')
        
        plt.grid()
        
        plt.legend()
        
        plt.savefig(plot_file)
    
        plt.close(1)
        
        print 'Colour-colour diagram output to '+plot_file
        
    else:
            
        print 'Warning: Insufficient data for colour-colour diagram'

if __name__ == '__main__':
    
    plot_colour_diagrams()
    