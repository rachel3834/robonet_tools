# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 15:13:18 2018

@author: rstreet
"""

import os
import sys
from astropy import coordinates as coords
from astropy import units as u
from astropy import wcs
from astropy.io import fits
import glob
import numpy as np

def process_image_set():
    """Function to update the image headers of a set of images"""
    
    params = get_args()
    
    file_list = glob.glob(os.path.join(params['dir_path'],'*.fits'))
    
    target = coords.SkyCoord(params['target_ra']+' '+params['target_dec'],
                             unit=(u.hourangle, u.deg),frame='icrs')
    
    for image in file_list:
        
        crop_image(image,target,params)


def crop_image(image_file,target,params):
    """Function to crop a single image around the target location given"""

    new_wcs_keys = ['WCSIMCAT', 'WCSMATCH', 'WCSNREF', 'WCSTOL', 'RA', 'DEC',
                    'WEQUINOX', 'WEPOCH', 'RADECSYS', 'CDELT1', 'CDELT2', 
                    'CROTA1', 'CROTA2', 'SECPIX1', 'SECPIX2', 'WCSSEP', 
                    'EQUINOX', 'CD1_1','CD1_2','CD2_1','CD2_2','EPOCH',
                    'IMWCS']
    
    header = fits.getheader(image_file)
    data = fits.getdata(image_file)
    
    image_centre = coords.SkyCoord(header['RA']+' '+header['DEC'],
                             unit=(u.hourangle, u.deg),frame='icrs')
    
    crop_half_width_pix = (params['sub_image_width']*60.0)/float(header['SECPIX1'])
    
    w = wcs.WCS(naxis=2)
    w.wcs.cdelt = [header['CDELT1'],header['CDELT2']]
    w.wcs.crval = [image_centre.ra.deg,image_centre.dec.deg]
    w.wcs.crpix = [header['NAXIS1']/2, header['NAXIS2']/2]
    w.wcs.cd = np.array([[header['CD1_1'],header['CD1_2']],[header['CD2_1'],header['CD2_2']]])
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.radesys = 'ICRS'
    w.wcs.cunit = [u.deg,u.deg]
    w.wcs.equinox = header['EQUINOX']
    
    target_pixel = wcs.utils.skycoord_to_pixel(target,w,origin=1)
    
    xlimits = [ int(target_pixel[0] - crop_half_width_pix), 
               int(target_pixel[0] + crop_half_width_pix) ]
    ylimits = [ int(target_pixel[1] - crop_half_width_pix), 
               int(target_pixel[1] + crop_half_width_pix) ]
    
    if xlimits[0] > 0.0 and ylimits[0] > 0.0 and \
        xlimits[1] < header['NAXIS2'] and ylimits[1] < header['NAXIS1']:
            
            new_data = data[ylimits[0]:ylimits[1],xlimits[0]:xlimits[1]]
            
            header['NAXIS1'] = data.shape[0]
            header['NAXIS2'] = data.shape[1]
            header['RA'] = params['target_ra']
            header['DEC'] = params['target_dec']
            header['CRPIX1'] = header['NAXIS1']/2.0
            header['CRPIX2'] = header['NAXIS2']/2.0
            
            new_image_file = image_file.replace('.fits','_crop.fits')
            
            fits.writeto(new_image_file, new_data, header, overwrite=True)
                
def get_args():
    """Function to obtain or prompt for commandline arguments"""
    
    params = {}
    
    if len(sys.argv) == 1:
        
        params['dir_path'] = raw_input('Please enter the path to the image directory: ')
        params['target_ra'] = raw_input('Please enter the target RA: ')    
        params['target_dec'] = raw_input('Please enter the target Dec: ')    
        params['sub_image_width'] = float(raw_input('Please enter the width of the cropped image [arcmin]: '))
        
    else:
        
        params['dir_path'] = sys.argv[1]
        params['target_ra'] = sys.argv[2]
        params['target_dec'] = sys.argv[3]
        params['sub_image_width'] = float(sys.argv[4])
        
    return params


if __name__ == '__main__':
    
    process_image_set()