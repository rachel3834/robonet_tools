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
    
    if os.path.isdir(dir_path):
        
        opt = raw_input('Are you sure you want to remove all image data from this location? Y or N: ')
        
        if str(opt).upper() == 'Y':
            
            remove_raw_image_data(dir_path)

        print('Completed clean up')
        
        
def remove_raw_image_data(dir_path):
    """Function to walk through the directory tree of the ROMEREA project
    and remove all raw images."""
    
    top_dirs = ['bad', 'good']
    
    for d in top_dirs:
        
        for s in [ 'rea', 'rome' ]:
            
            for t in [ 'coj', 'cpt', 'lsc' ]:
                
                cameras = glob.glob(os.path.join(dir_path,d,s,t,'*'))
                
                for c in cameras:
                    
                    filters = glob.glob(os.path.join(dir_path,d,s,t,c,'*'))
                    
                    for f in filters:
                        
                        fields = glob.glob(os.path.join(dir_path,d,s,t,c,f,'*'))
                        
                        for field in fields:
                            
                            files = glob.glob(os.path.join(dir_path,d,s,t,c,f,field,'*.fits'))
                            
                            for i in files:
                                
                                print('Removing '+i)
                                #os.remove(i)


if __name__ == '__main__':
    
    clear_old_season_data_dirs()