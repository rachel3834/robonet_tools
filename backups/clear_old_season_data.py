# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 10:09:58 2018

@author: rstreet
"""

import os
import sys

def clear_old_season_data_dirs():
    """WARNING: DELETES DATA
    
    This function is designed to clear the raw data that accumulates in the
    operational observing system directions.  
    """
    
    print('WARNING: THIS FUNCTION WILL REMOVE DATA!')
    dir_path = raw_input('Please enter the path to the top level data directory to be cleared: ')
    
    if os.sys.path(dir_path):
        
        opt = raw_input('Are you sure you want to remove all image data from this location? Y or N: ')
        
        if str(opt).upper() == 'Y':
            
            remove_raw_image_data(dir_path)

        print('Completed clean up')
        
        
def remove_raw_image_data(dir_path):
    """Function to walk through the directory tree of the ROMEREA project
    and remove all raw images."""
    
    for root, dirs, files in os.walk(dir_path):

        for name in files:
            
            local_file_path = os.path.join(root,name)
            
            if '.fits' in local_file_path:
            
                print('Removing '+local_file_path)
                
                #os.remove(local_file_path)


if __name__ == '__main__':
    
    clear_old_season_data_dirs()