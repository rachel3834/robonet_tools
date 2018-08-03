# -*- coding: utf-8 -*-
"""
Created on Sun Jul 29 15:42:40 2018

@author: rstreet
"""

from os import path
from sys import argv
import numpy as np
from pyDANDIA import pipeline_setup
from pyDANDIA import logs
from pyDANDIA import metadata
from pyDANDIA import catalog_utils
from astropy.table import Table
import logging

def calc_calibrated_lightcurve():
    """Function to apply a measured photometric offset to calibrate the 
    instrumental magnitudes provided in a lightcurve"""
    
    params = get_args()
    
    setup = pipeline_setup.pipeline_setup(params)
    
    log = logs.start_stage_log( setup.red_dir, 'lc_calib' )

    log.info('Calibrating lightcurve.  Initial parameters:')    
    for key, value in params.items():
        log.info(key+': '+str(value))
        
    (reduction_metadata, params, star_catalog) = fetch_metadata(setup,params,log)
    
    lc_data = read_rbn_lc(params['lc_file'],log)
    
    params = calc_mag_offset(params,lc_data,reduction_metadata,star_catalog,log)
    
    apply_mag_offset_to_lc(params,lc_data,log)
    
    logging.shutdown()
    
def apply_mag_offset_to_lc(params,lc_data,log):
    """Function to apply the calculated magnitude offset to the lightcurve data
    and output"""
    
    var_delta_m = params['sig_delta_m']*params['sig_delta_m']
    
    f = open(params['cal_lc_file'],'w')
    
    f.write('# Image : HJD (days) : df (ADU) : df err (ADU) : rf (ADU) : rf err (ADU) : mag : mag err : cal mag : cal mag err : Exposure Time (s) : FWHM (pixel) : Sky (ADU) : Airmass : Photometric Scale Factor : Chi Squared Per Pixel :     dx (pix) : dy (pix) : a : b : c : d\n')
    
    for i in range(0,len(lc_data),1):
        
        row = ''
        for j in range(0,8,1):
            row += lc_data[i,j]+' '
        
        cal_mag = float(lc_data[i,6]) + params['delta_m']
        
        cal_mag_err = np.sqrt( float(lc_data[i,7])*float(lc_data[i,7]) + var_delta_m )
        
        row += str(round(cal_mag,5))+' '+str(round(cal_mag_err,9))+' '
        
        for j in range(9,20,1):
            row += lc_data[i,j]+' '
        
        f.write(row+'\n')
    
    f.close()
    
    log.info('Applied magnitude correction to lightcurve data, output to '+params['cal_lc_file'])
    
def calc_mag_offset(params,lc_data,reduction_metadata,star_catalog,log):
    """When comparing the magnitudes measured from DanDIA in the lightcurve
    with those from pyDANDIA in the calibration data, the total correction 
    to the magnitude consists of:
    
    corr_mag = lc_mag + delta_m
    
    where delta_m = lc_mag[ref_frame] - cal_ref_mag[star_catalog]
    
    and:
    lc_mag[ref_frame] = magnitude of star in lightcurve entry for reference frame
    cal_ref_mag[star_catalog] = magnitude of star from metadatas star_catalog
    
    Note: The reference frame used to produce both the lightcurve and the
    metadata must be the same.  
    """
    
    # First find the lightcurve entry for the reference image (which must be 
    # common to both lightcurve and reference frame reductions).  
    
    images = lc_data[:,0].tolist()
    
    refidx = images.index(params['refimage'].replace('.fits','_crop'))
    
    log.info('Reference image at index '+str(refidx)+' in lightcurve')
    
    lc_mag = float(lc_data[refidx,6])
    lc_mag_err = float(lc_data[refidx,7])
    
    log.info('Stars measured lightcurve mag in ref frame: '+\
            str(lc_mag)+' +/- '+str(lc_mag_err))
            
    sidx = int(params['staridx'])

    cal_ref_mag = star_catalog['cal_ref_mag'][sidx-1]
    cal_ref_mag_err = star_catalog['cal_ref_mag_err'][sidx-1]
    
    log.info('Stars calibrated pyDANDIA mag in ref frame: '+\
            str(cal_ref_mag)+' +/- '+str(cal_ref_mag_err))
    
    params['delta_m'] = cal_ref_mag - lc_mag
    params['sig_delta_m'] = np.sqrt( (cal_ref_mag_err*cal_ref_mag_err) + \
                            (lc_mag_err*lc_mag_err))
    
    log.info('Photometric offset, delta_m: '+str(params['delta_m'])+' +/- '+str(params['sig_delta_m']))
    log.info('Uncertainties in offset calculation:')
    log.info('Error on calibrated reference magnitude: '+str(cal_ref_mag_err))
    log.info('Error on lightcurve magnitude: '+str(lc_mag_err))
    
    return params
    
def fetch_metadata(setup,params,log):
    """Function to extract the information necessary for the photometric
    calibration from a metadata file, adding information to the params 
    dictionary"""
    
    reduction_metadata = metadata.MetaData()
    
    reduction_metadata.load_a_layer_from_file( setup.red_dir, 
                                              params['metadata'],
                                              'data_architecture' )
    reduction_metadata.load_a_layer_from_file( setup.red_dir, 
                                              params['metadata'],
                                              'star_catalog' )
    reduction_metadata.load_a_layer_from_file( setup.red_dir, 
                                              params['metadata'],
                                              'phot_calib' )
                                              
    params['refimage'] = reduction_metadata.data_architecture[1]['REF_IMAGE'][0]
    
    log.info('Using data for reference image: '+params['refimage'])
    
    star_catalog = Table()
    star_catalog['star_index'] = reduction_metadata.star_catalog[1]['star_index']
    star_catalog['mag'] = reduction_metadata.star_catalog[1]['ref_mag']
    star_catalog['mag_err'] = reduction_metadata.star_catalog[1]['ref_mag_err']
    star_catalog['cal_ref_mag'] = reduction_metadata.phot_calib[1]['cal_ref_mag']
    star_catalog['cal_ref_mag_err'] = reduction_metadata.phot_calib[1]['cal_ref_mag_err']
    
    log.info('Extracted star catalog')
    
    return reduction_metadata, params, star_catalog
    
def read_rbn_lc(lc_file,log):
    """Function to read in the lightcurve data in RoboNet standard format"""
    
    if path.isfile(lc_file) == False:
        
        print('ERROR: Cannot find lightcurve file '+lc_file)
        log.info('ERROR: Cannot find lightcurve file '+lc_file)
        
        logging.shutdown()

        exit()
        
    f = open(lc_file,'r').readlines()
    
    header = f[1]
    data = []
    
    for line in f[1:]:
        entries = line.replace('\n','').split()
        data.append(entries)
    
    data = np.array(data)
    
    log.info('Read in lightcurve file')
    
    return data
    
def get_args():
    """Function to acquired the required input parameters"""
    
    params = {}
    
    if len(argv) > 1:
        
        params['lc_file'] = argv[1]
        params['red_dir'] = argv[2]
        params['log_dir'] = argv[3]
        params['metadata'] = argv[4]
        params['staridx'] = argv[5]
        params['cal_lc_file'] = argv[8]
        
    else:
        
        params['lc_file'] = raw_input('Please enter the path to the lightcurve file: ')
        params['red_dir'] = raw_input('Please enter the path to the reduction directory: ')
        params['log_dir'] = raw_input('Please enter the path to the log directory: ')
        params['metadata'] = raw_input('Please enter the path to the metadata file: ')
        params['staridx'] = raw_input('Please enter the index number of the star in the metadatas starcatalog table: ')
        params['cal_lc_file'] = raw_input('Please enter the output file for the calibrated lightcurve: ')
    
    return params


if __name__ == '__main__':
    
    calc_calibrated_lightcurve()
    