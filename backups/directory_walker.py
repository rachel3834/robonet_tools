# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 12:31:18 2018

@author: rstreet
"""

from os import path
from sys import argv
import glob
import compress_data_products

def walk_ml2016_structure(top_dir):
    """Function to walk the data reduction structure used for the ML2016 data, 
    identify directories to be compressed and compress them."""
    
    if path.isdir(top_dir) == False:
        print('Error: Cannot find top level directory '+top_dir)
        exit()
    
    event_dirs = glob.glob(path.join(top_dir,'*'))
    print(event_dirs)
    
    for event in event_dirs:
        
        for binning in ['11', '22']:
            
            if path.isdir(path.join(event,binning)):
                print(path.join(event,binning))
                camera_dirs = glob.glob(path.join(event,binning,'*'))
                print(camera_dirs)
                for camera in camera_dirs:
                    
                    compress_dandia_reduced_data_products(camera)


if __name__ == '__main__':
    
    if len(argv) == 1:
        
        top_dir = raw_input('Please enter the top level directory path: ')
        
    else:
        
        top_dir = argv[1]
    
    walk_ml2016_structure(top_dir)
    