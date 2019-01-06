# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 16:42:56 2019

@author: rstreet
"""

from os import path, makedirs
from sys import argv
from astropy.io import fits
from datetime import datetime
import glob
from shutil import copy

def search_lco_local_archive():
    """Function to search the LCO local data archive, primarily for data
    taken prior to 2014 that is not available through the online Cloud archive
    service
    
    Data structure is:
    /net/archive/archive/engineering/site-code/camera/day-date/
    -> cat/ preview/ processed/ raw/ wcs/
    -> -e91 frames are in processed, fits.fz compressed
    """
    
    search_params = get_args()
    
    active_sites = ['coj', 'cpt', 'elp', 'lsc', 'ogg']
    active_cameras_prefix = ['fl', 'fa', 'kb']
    header_params = ['DATE-OBS', 'SITEID', 'ENCID', 'TELESCOP', 'INSTRUME', 
                     'GROUPID', 'OBJECT']
    
    search_dirs = locate_search_subdirs(search_params, active_sites, 
                                        active_cameras_prefix)
    
    scan_for_data(search_params, search_dirs, header_params)
    
def locate_search_subdirs(search_params, active_sites, active_cameras_prefix):
    
    search_dirs = []
                
    for site in active_sites:
        
        for camera_prefix in active_cameras_prefix:

            cameras = glob.glob(path.join(search_params['archive_path'],\
                                            site,camera_prefix+'*'))
                                            
            for camera_dir in cameras:
                
                date_list = glob.glob(path.join(camera_dir,'20??????'))
                
                for dpath in date_list:
                    
                    ddate = path.basename(dpath)
                    ddate = datetime.strptime(ddate, '%Y%m%d')
                    
                    if ddate >= search_params['start_date'] and \
                           ddate <= search_params['end_date']:
                           
                           search_dirs.append(dpath)
    
    print('Identified '+str(len(search_dirs))+' data directories to search')
    
    return search_dirs
    
def get_args():
    
    search_params = {'copy': False}
    
    if len(argv) == 1:
        search_params['archive_path'] = raw_input('Please enter the path to the top-level archive directory [e.g. /net/archive/archive/engineering/]: ')
        search_params['start_date'] = raw_input('Please enter the start date for a search: ')
        search_params['end_date'] = raw_input('Please enter the end date for a search: ')
        search_params['log_dir'] = raw_input('Please enter the directory path for output: ')
        search_params['search_key'] = raw_input('Please enter header parameter to be used to identify data [OBJECT or GROUPID]: ')
        search_params['search_substr'] = raw_input('Please enter the sub-string which must be included in the search_key value for data to be selected [e.g. OGLE, no wildcards required]: ')
        opt = raw_input('Copy selected data?  Y or n: ')
        
        if 'Y' in str(opt).upper():
            search_params['copy'] = True
            
    else:
        search_params['archive_path'] = argv[1]
        search_params['start_date'] = argv[2]
        search_params['end_date'] = argv[3]
        search_params['log_dir'] = argv[4]
        search_params['search_key'] = argv[5]
        search_params['search_substr'] = argv[6]
        
        if '-copy' in argv:
            search_params['copy'] = True
            
    for d in ['start_date', 'end_date']:
        t = datetime.strptime(search_params[d], '%Y-%m-%d')
        search_params[d] = t
    
    print('Searching data archive with the following parameters:')
    print(search_params)
    
    return search_params
    
def get_image_header_info(image_path, params):
    """Function to extract the values of a set of parameter keywords from a
    compressed FITS image header
    """
    
    with fits.open(image_path) as hdul: 
        header = hdul[1].header

    header_info = {}
    
    for key in params:
        header_info[key] = header[key]

    return header_info

def scan_for_data(search_params, search_dirs, header_params):
    
    log_path = path.join(search_params['log_dir'], 
                         'archive_walker_search_results.log')
    
    data_log = open(log_path,'w')
    data_log.write('# ' + ' '.join(header_params) + ' Path\n')
    
    for d in search_dirs:
        
        #dpath = path.join(d,'processed')
        dpath = path.join(d,'raw')
        
        file_list = glob.glob(path.join(dpath, '*.fits.fz'))
        
        print(' -> Found '+str(len(file_list))+' frames in '+dpath)
        
        for f in file_list:
            
            header_info = get_image_header_info(f, header_params)
            
            if search_params['search_substr'] in header_info[search_params['search_key']]:
                output = ''
                for key in header_params:
                    output += ' ' + str(header_info[key])
                
                output += ' ' + f + '\n'
            
            data_log.write(output)
    
    data_log.close()

def copy_data_from_archive(search_params,file_path):
    
    if search_params['copy'] == True:
        
        (src,file_name) = str(path.split(file_path))
        dest = src.replace(search_params['archive_path'],search_params['log_dir'])
        
        if path.isdir(dest) == False:
            makedirs(dest)
        
        copy(file_path, path.join(dest,file_name))
        
        print('-> Transferred '+file_name+' to '+dest)
        
if __name__ == '__main__':
        
    search_lco_local_archive()