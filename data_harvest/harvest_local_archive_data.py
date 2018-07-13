# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 13:40:54 2018

@author: rstreet
"""

from sys import argv
from os import path
from datetime import datetime
from commands import getstatusoutput
from shutil import copy

def archive_spider():
    """Walker function to explore the local LCO science data archive, find
    ROME/REA data and copy them to local datastore."""
    
    params = get_args()
    
    if path.isdir(params['archive_path']) == False:
        
        print('Error: Cannot find archive data at '+params['archive_path'])
        exit()
    
        for site in params['site_list']:
            
            camera_list = glob.glob(path.join(params['archive_path'],site,'fl*'))
            
            for camera_dir in camera_list:
                
                camera = path.basename(camera_dir)
                
                date_list = glob.glob(path.join(camera_dir,'20??????'))
                
                for date_dir in date_list:
                    
                    date = datetime.strptime(date_dir,"%Y%m%d")
                    
                    if date >= params['rome_start_date']:
                        
                        file_list = glob.glob(path.join(date_dir,'????????-'+camera+'-????????-????-e91.fits.fz'))
                        
                        for f in file_list:
                                                        
                            (rome_image,image_target) = check_matching_image(params,f)
                            
                            copy_image_to_local_archive(params,f,image_target)
                            
def get_args():

    params = {}
    
    if len(argv) == 1:
        
        params['archive_path'] = raw_input('Please enter the path to the top-level archive directory [e.g. /archive/engineering/]: ')
        params['field'] = raw_input('Please enter which field you want to search for (or ALL): ')
    
    else:
    
        params['archive_path'] = argv[1]
        params['field'] = argv[2]
    
    params['field'] = str(params['field']).upper()
    params['site_list'] = [ 'coj', 'cpt', 'lsc' ]
    params['rome_start_date'] = datetime.strptime('2017-04-01','%Y-%m-%d')
    
    params['local_archive'] = {
                            'ROME-FIELD-01': '/data09/ROME_Archive/ROME-FIELD-01',
                            'ROME-FIELD-02': '/data09/ROME_Archive/ROME-FIELD-02',
                            'ROME-FIELD-03': '/data08/ROME_Archive/ROME-FIELD-03',
                            'ROME-FIELD-04': '/data09/ROME_Archive/ROME-FIELD-04',
                            'ROME-FIELD-05': '/data09/ROME_Archive/ROME-FIELD-05',
                            'ROME-FIELD-06': '/data09/ROME_Archive/ROME-FIELD-06',
                            'ROME-FIELD-07': '/data10/ROME_Archive/ROME-FIELD-07',
                            'ROME-FIELD-08': '/data10/ROME_Archive/ROME-FIELD-08',
                            'ROME-FIELD-09': '/data10/ROME_Archive/ROME-FIELD-09',
                            'ROME-FIELD-10': '/data10/ROME_Archive/ROME-FIELD-10',
                            }
    
    if params['field'] not in params['local_archive'].keys():
        
        print('Warning: No local archive location yet established for '+params['field'])
        exit()
        
    return params
    
def check_matching_image(params,f):
    """Function to check whether the image given matches the fields selected
    for transfer to the local archive or not"""
    
    image_target = getstatusoutput('gethead OBJECT '+f)
    
    if params['field'] != 'ALL':
        
        template = params['field']
        
    else:
        template = 'ROME-FIELD'
    
    if template in coutput:
        
        return True, image_target
        
    else:
        
        return False, image_target
        
def copy_image_to_local_archive(params,f,image_target):
    """Function to copy the selected frame to the appropriate location in the
    local data repository, distributing the data across disks as necessary"""
    
    if image_target in params['local_archive'].keys():
        
        copy(f,params['local_archive'][image_target])
        
        print('Copied '+image_target+' '+path.basename(f)+' -> '+params['local_archive'][image_target])
        
    
if __name__ == '__main__':
    
    archive_spider()