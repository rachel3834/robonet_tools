# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 09:07:46 2018

@author: rstreet
"""

import os
import sys
import pickle
import glob
from astropy.io import fits

def process_image_set():
    """Function to update the image headers of a set of images"""
    
    (dir_path,wcs_path) = get_args()
    
    wcs_data = load_header_data(wcs_path)
    
    file_list = glob.glob(os.path.join(dir_path,'*.fits'))
    
    for image_file in file_list:
        
        update_image_wcs(image_file,wcs_data)
        
def get_args():
    """Function to obtain or prompt for commandline arguments"""
    
    if len(sys.argv) == 1:
        
        dir_path = raw_input('Please enter the path to the image directory: ')
    
        wcs_path = raw_input('Please enter the path to the directory containing the WCS header pickle files: ')
        
    else:
        
        dir_path = sys.argv[1]
        
        wcs_path = sys.argv[2]

    return dir_path, wcs_path
    
def update_image_wcs(image_file,wcs_data):
    """Function to update the WCS header keywords of a given image file"""
    
    old_wcs_keys = ['CTYPE1','CTYPE2','CRPIX1','CRPIX2','CRVAL1','CRVAL2',
                'CUNIT1','CUNIT2','CD1_1','CD1_2','CD2_1','CD2_2']    
    new_wcs_keys = ['WCSIMCAT', 'WCSMATCH', 'WCSNREF', 'WCSTOL', 'RA', 'DEC',
                    'WEQUINOX', 'WEPOCH', 'RADECSYS', 'CDELT1', 'CDELT2', 
                    'CROTA1', 'CROTA2', 'SECPIX1', 'SECPIX2', 'WCSSEP', 
                    'EQUINOX', 'CD1_1','CD1_2','CD2_1','CD2_2','EPOCH',
                    'IMWCS']
                    
    if os.path.isfile(image_file):
        
        if os.path.basename(image_file) in wcs_data.keys():
            
            wcs_header = wcs_data[os.path.basename(image_file)]
        
            header = fits.getheader(image_file)
            data = fits.getdata(image_file)
            
            for key in old_wcs_keys:
                
                if key in header:
                    
                    header.pop(key,None)
            
            hdr_ok = True
            
            for key in new_wcs_keys:
                
                try:
                    
                    header[key] = wcs_header[key]
                
                except KeyError:
                
                    print('ERROR: Missing header keys from WCS solution')
                    
                    hdr_ok = False
                    
            header['WCSERR'] = 0
            
            if hdr_ok:
                
                new_image_file = image_file.replace('.fits','_wcs.fits')
            
                fits.writeto(new_image_file, data, header)
            
        else:
        
            print('ERROR: No entry for '+os.path.basename(image_file)+' in WCS data')
            
            
def load_header_data(wcs_path):
    """Function to load image header data from pickled dictionary files"""

    file_list = [ 'coj_header_imwcs_update.pickle',
                 'cpt_header_imwcs_update.pickle',
                 'lsc_header_imwcs_update.pickle' ]
    
    wcs_data = {}
    
    for f in file_list:         
        
        file_path = os.path.join(wcs_path,f)
        
        with open(file_path, "rb") as input_file:
    
            data = pickle.load(input_file)
            
            wcs_data.update(data)
    
    return wcs_data
    

if __name__ == '__main__':
    
    process_image_set()