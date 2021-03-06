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
    
    header = fits.getheader(file_list[0])
    
    pixscale = get_pixel_scale(header)
    
    if params['option'] == 'arcmin':
        crop_half_width_pix = int((params['sub_image_width_arcmin']*60.0)/pixscale)/2
    else:
        crop_half_width_pix = int(params['sub_image_width_pixels']/2.0)
        
    print('Half width of sub-image will be '+str(crop_half_width_pix)+' pix')
    
    if params['xoffset'] != 0.0 or params['yoffset'] != 0.0:
        
        print('Applying offset to cropped section: dx='+str(params['xoffset'])+\
                                                ', dy='+str(params['yoffset']))
                            
    for image in file_list:
        
        crop_image(image,target,params,crop_half_width_pix)

def get_pixel_scale(header):
    """Function to return the pixel scale from the image header, 
    accommodating different keyword conventions"""
    
    keywords = [ 'SECPIX1', 'PIXSCALE' ]
    
    pixscale = None
    
    for key in keywords:
        
        if key in header.keys():
            
            pixscale = header[key]
    
    return float(pixscale)
    
def crop_image(image_file,target,params,crop_half_width_pix):
    """Function to crop a single image around the target location given"""

    new_wcs_keys = ['WCSIMCAT', 'WCSMATCH', 'WCSNREF', 'WCSTOL', 'RA', 'DEC',
                    'WEQUINOX', 'WEPOCH', 'RADECSYS', 'CDELT1', 'CDELT2', 
                    'CROTA1', 'CROTA2', 'SECPIX1', 'SECPIX2', 'WCSSEP', 
                    'EQUINOX', 'CD1_1','CD1_2','CD2_1','CD2_2','EPOCH',
                    'IMWCS']
    
    hdu = fits.open(image_file)
    header = hdu[0].header
    data = hdu[0].data
    cat = None
    bpm = None
    if len(hdu) == 2:
        bpm = hdu[1]
    elif len(hdu) == 3:
        cat = hdu[1]
        bpm = hdu[2]
    
    image_centre = coords.SkyCoord(header['RA']+' '+header['DEC'],
                             unit=(u.hourangle, u.deg),frame='icrs')
    
    
    w = wcs.WCS(naxis=2)
    if 'CDELT1' in header:
        w.wcs.cdelt = [header['CDELT1'],header['CDELT2']]
    w.wcs.crval = [image_centre.ra.deg,image_centre.dec.deg]
    w.wcs.crpix = [header['NAXIS1']/2, header['NAXIS2']/2]
    w.wcs.cd = np.array([[header['CD1_1'],header['CD1_2']],[header['CD2_1'],header['CD2_2']]])
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.radesys = 'ICRS'
    w.wcs.cunit = [u.deg,u.deg]
    if 'EQUINOX' in header:
        w.wcs.equinox = header['EQUINOX']
    
    target_pixel = wcs.utils.skycoord_to_pixel(target,w,origin=1)
    
    xlimits = [ int(target_pixel[0] - crop_half_width_pix + params['xoffset']), 
               int(target_pixel[0] + crop_half_width_pix + params['xoffset']) ]
    ylimits = [ int(target_pixel[1] - crop_half_width_pix + params['yoffset']), 
               int(target_pixel[1] + crop_half_width_pix + params['yoffset']) ]
    
    print(os.path.basename(image_file)+' cropsize= '+str(xlimits[1] - xlimits[0])+
            ' x '+str(ylimits[1] - ylimits[0])+' pix')
    
    if xlimits[0] > 0.0 and ylimits[0] > 0.0 and \
        xlimits[1] < header['NAXIS2'] and ylimits[1] < header['NAXIS1']:
            
            new_data = data[ylimits[0]:ylimits[1],xlimits[0]:xlimits[1]]
            
            if bpm != None:
                bpm.data = bpm.data[ylimits[0]:ylimits[1],xlimits[0]:xlimits[1]]
                
            header['NAXIS1'] = new_data.shape[0]
            header['NAXIS2'] = new_data.shape[1]
            header['RA'] = params['target_ra']
            header['DEC'] = params['target_dec']
            header['CRPIX1'] = header['NAXIS1']/2.0
            header['CRPIX2'] = header['NAXIS2']/2.0
            header['CRVAL1'] = image_centre.ra.deg
            header['CRVAL2'] = image_centre.dec.deg
            
            new_image_file = image_file.replace('.fits','_crop.fits')
            
            new_hdu = fits.HDUList()
            new_hdu.append(fits.PrimaryHDU(data=new_data, header=header))
            
            if cat != None:
                try:
                    new_hdu.append(fits.BinTableHDU(data=cat.data, header=cat.header))
                except TypeError:
                    pass

            if bpm != None:
                new_hdu.append(fits.ImageHDU(data=bpm.data, header=bpm.header))

            new_hdu.writeto(new_image_file, overwrite=True)
            
    else:
        
        print('--> Warning: crop limits exceed image boundaries ('+\
                repr(xlimits)+', '+repr(ylimits)+')')
                
def get_args():
    """Function to obtain or prompt for commandline arguments"""
    
    params = {}
    
    if len(sys.argv) == 1:
        
        params['dir_path'] = input('Please enter the path to the image directory: ')
        params['target_ra'] = input('Please enter the target RA: ')
        params['target_dec'] = input('Please enter the target Dec: ')
        params['option'] = input('Specify width of cropped image in arcmin or pixels? [arcmin,pixels]: ')
        params['option'] = str(params['option']).lower()
        if params['option'] == 'arcmin':
            params['sub_image_width_arcmin'] = float(input('Please enter the width of the cropped image [arcmin]: '))
        else:
            params['sub_image_width_pixels'] = float(input('Please enter the width of the cropped image [pixels]: '))
        params['xoffset'] = float(input('Please enter the box offset in the X-axis: '))
        params['yoffset'] = float(input('Please enter the box offset in the Y-axis: '))
        
    else:
        
        params['dir_path'] = sys.argv[1]
        params['target_ra'] = sys.argv[2]
        params['target_dec'] = sys.argv[3]
        params['option'] = sys.argv[4]
        if params['option'] == 'arcmin':
            params['sub_image_width_arcmin'] = float(sys.argv[5])
        else:
            params['sub_image_width_pixels'] = float(sys.argv[5])
        params['xoffset'] = float(sys.argv[6])
        params['yoffset'] = float(sys.argv[7])
        
    return params


if __name__ == '__main__':
    
    process_image_set()