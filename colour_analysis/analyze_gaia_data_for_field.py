# -*- coding: utf-8 -*-
"""
Created on Sun Jul 29 14:27:23 2018

@author: rstreet
"""

from sys import argv
from os import path, remove
import glob
import json
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import wget
import analyze_gaia_results
import select_colour_sample

def analyze_gaia_data_for_field():
    """Function to analyze the data returned from a search of the Gaia archive
    around the location of a target of known location, firstly to check
    for a match and if one is found, calculate the measured distance, if
    possible"""

    params = get_args()

    photometry = select_colour_sample.read_colour_photometry_file(params)

    photometry = create_gaia_sources(params, photometry)

    photometry = query_bailer_jones_catalog_field(params,photometry)

    output_annotated_photometry(params, photometry)

def create_gaia_sources(params, photometry):
    """Function to perform a single object search of the Gaia catalogue,
    using its API service for Table Access Protocol (TAP)
    """

    photometry['GaiaSources'] = [ None ]*len(photometry)

    for j,star in enumerate(photometry):
        if 'none' not in str(star['Gaia_ID'].data):
            star_data = {'source_id': star['Gaia_ID'].data,
                     'ra': star['ra_deg'].data,
                     'dec': star['dec_deg'].data}
            source = analyze_gaia_results.GaiaSource(star_data)
            photometry['GaiaSources'][j] = source

    return photometry

def query_bailer_jones_catalog_field(params,photometry):
    """Function to query the API for the catalogue of Gaia source positions
    published in Bailer-Jones et al. (2018), AJ, 156, 58"""

    for j,source in enumerate(photometry['GaiaSources'].data):
        source = analyze_gaia_results.query_bailer_jones_catalog(params,source)
        photometry['GaiaSources'][j] = source

    return photometry


def ensure_matching_star(params,source):
    """Function to compare the target location with those of the returned
    search results to look for a match"""

    target = SkyCoord(params['target_ra'],params['target_dec'],unit=(u.hourangle,u.deg))

    gaia_star = SkyCoord(source.ra,source.dec,unit=(u.deg,u.deg))

    separation = target.separation(gaia_star)

    setattr( source, 'separation', (separation.value * 3600.0) )

    if source.separation > params['max_separation']:
        source = None

    return source

def get_args():
    """Function to read in commandline arguments if available or otherwise
    ask the user for the required parameters"""

    params = {}

    if len(argv) > 1:

        params['phot_file'] = argv[1]
        params['red_dir'] = argv[2]
        params['max_separation'] = argv[3]

    else:
        params['phot_file'] = input('Please enter the path to the colour photometry file: ')
        params['red_dir'] = input('Please enter the path to the data directory: ')
        params['max_separation'] = input('Please enter the maximum acceptable separation to match star sky positions [arcsec]: ')

    return params

def output_annotated_photometry(params, photometry):

    if str(config['photometry_data_file']).lower() != 'none':

        fname = params['phot_file'].split('.')+'_gaia.txt'
        f = open(path.dirpath(params['phot_file'],fname), 'w')
        f.write('# All measured floating point quantities in units of magnitude\n')
        f.write('# Selected indicates whether a star lies within the selection radius of a given location, if any.  1=true, 0=false\n')
        f.write('# Star   x_pix    y_pix   ra_deg   dec_deg   g  sigma_g    r  sigma_r    i  sigma_i   (g-i)  sigma(g-i) (g-r)  sigma(g-r)  (r-i) sigma(r-i)  Selected  Gaia_ID Distance Distance_lo Distance_hi\n')

        for j in range(0,len(photometry['i']),1):
            source = photometry['GaiaSources'][j]
            f.write( str(photometry['Star'][j])+' '+\
                        str(photometry['x_pix'][j])+' '+str(photometry['y_pix'][j])+' '+\
                        str(photometry['ra'][j])+' '+str(photometry['dec'][j])+' '+\
                        str(photometry['g'][j])+' '+str(photometry['sigma_g'][j])+' '+\
                        str(photometry['r'][j])+' '+str(photometry['sigma_r'][j])+' '+\
                        str(photometry['i'][j])+' '+str(photometry['sigma_i'][j])+' '+\
                        str(photometry['(g-i)'][j])+' '+str(photometry['sigma(g-i)'][j])+' '+\
                        str(photometry['(g-r)'][j])+' '+str(photometry['sigma(g-r)'][j])+' '+\
                        str(photometry['(r-i)'][j])+' '+str(photometry['sigma(r-i)'][j])+' '+\
                        str(photometry['selected'][j])+' '+str(photometry['gaia_source_id'])+' '+\
                        str(source.r_est)+' '+str(source.r_lo)+' '+str(source.r_hi)+'\n' )

        f.close()


if __name__ == '__main__':

    analyze_gaia_data_for_field()
