# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 11:03:41 2019

@author: rstreet
"""

import os
import sys
import glob
import subprocess
import compress_data_products

def seek_and_compress_fits(top_level_dir):
    """Function to walk over every subdirectory of the directory given to 
    identify FITS files and compress them using fpack"""
    
    dir_listing = os.walk(top_level_dir)
    
    for sub_dir in os.walk(top_level_dir):
        
        compress_data_products.bzip2_image_dir(sub_dir[0])
        

if __name__ == '__main__':
    
    if len(sys.argv) == 1:
        top_level_dir = input('Please enter the path to the top level directory: ')
    else:
        top_level_dir = sys.argv[1]
        
    seek_and_compress_fits(top_level_dir)