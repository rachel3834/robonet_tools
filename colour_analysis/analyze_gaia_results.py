# -*- coding: utf-8 -*-
"""
Created on Sun Jul 29 14:27:23 2018

@author: rstreet
"""

from sys import argv
from os import path
import json
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u

def analyze_gaia_data():
    """Function to analyze the data returned from a search of the Gaia archive
    around the location of a target of known location, firstly to check 
    for a match and if one is found, calculate the measured distance, if 
    possible"""
    
    params = get_args()
    
    gaia_data = load_gaia_archive_json(params['json_file'])
        
    separations = check_for_matching_star(params,gaia_data)
    
    gaia_data = calc_distance_from_parallax(gaia_data)
    
    print gaia_data
        
def calc_distance_from_parallax(gaia_data):
    """Function to calculate the distance to an object based on its parallax"""
    
    d = 1.0 / abs(gaia_data['parallax']/1000.0)
    
    sig_d = abs(gaia_data['parallax_error']/gaia_data['parallax']) * d
    
    gaia_data['distance'] = d
    gaia_data['distance_error'] = sig_d
    
    return gaia_data
    
def check_for_matching_star(params,gaia_data):
    """Function to compare the target location with those of the returned
    search results to look for a match"""
    
    target = SkyCoord(params['target_ra'],params['target_dec'],unit=(u.hourangle,u.deg))
    
    gaia_stars = SkyCoord(gaia_data['ra'],gaia_data['dec'],unit=(u.deg,u.deg))
    
    separations = target.separation(gaia_stars)
    
    for j,s in enumerate(separations):
    
        gaia_data['separation'][j] = s.value * 3600.0
    
    return gaia_data
    
def get_args():
    """Function to read in commandline arguments if available or otherwise
    ask the user for the required parameters"""
    
    params = {}
    
    if len(argv) > 1:
        
        params['json_file'] = argv[1]
        params['target_ra'] = argv[2]
        params['target_dec'] = argv[3]
        
    else:
        
        params['json_file'] = raw_input('Please enter the path to the JSON file from the Gaia archive:')
        params['target_ra'] = raw_input('Please enter the target RA in sexigesimal format: ')
        params['target_dec'] = raw_input('Please enter the target Dec in sexigesimal format: ')
    
    return params
    
def load_gaia_archive_json(json_file):
    """Function to read in search results from the Gaia archive around
    a given location in JSON format"""
    
    if path.isfile(json_file) == False:
        print('ERROR: Cannot find JSON-format search results from the Gaia archive, file '+json_file)
        exit()
    
    with open(json_file) as input_file:  
        data = json.load(input_file)
    
    cols = [ 'source_id', 'ra', 'dec', 'parallax', 'parallax_error', 'phot_g_mean_mag', 'distance', 'distance_error', 'separation' ]
    formats = [ 'S25', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8' ]
    idx = []
    
    for c in cols:
        for i,item in enumerate(data['metadata']):
            if item['name'] == c:
                idx.append(i)
    
    table_data = []
    
    for i,row in enumerate(data['data']):
        row_data = []
        for k,j in enumerate(idx):
            if 'S' in formats[k]:
                row_data.append(data['data'][i][j])
            elif 'f' in formats[k]:
                try:
                    row_data.append(float(data['data'][i][j]))
                except ValueError:
                    row_data.append(None)
                    
        row_data.append(0.0)
        row_data.append(0.0)
        row_data.append(0.0)
        
        table_data.append( tuple(row_data) )
    
    t = Table(rows=table_data, names=tuple(cols), 
              dtype=tuple(formats))
              
    return t
    
if __name__ == '__main__':
    
    analyze_gaia_data()