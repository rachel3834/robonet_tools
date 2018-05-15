# -*- coding: utf-8 -*-
"""
Created on Wed May  9 16:33:56 2018

@author: rstreet
"""
import os
import sys
from pyDANDIA import pipeline_setup
from pyDANDIA import logs
from pyDANDIA import metadata
from astropy.coordinates import SkyCoord
import vizier_tools
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

VERSION = 'calibrate_photometry_0.0.1'

def calibrate_photometry():
    """Function to calculate the photometric transform between the instrumental
    magnitudes produced by the pyDANDIA pipeline and catalog data."""
    
    params = get_args()

    setup = pipeline_setup.pipeline_setup(params)
    
    log = logs.start_stage_log( setup.red_dir, 'phot_calib', version=VERSION )
    
    (reduction_metadata, params) = fetch_metadata(setup,params,log)
        
    vphas_cat = fetch_catalog_sources_within_image(params,log)
    
    match_index = match_stars_by_position(reduction_metadata,vphas_cat,log)
    
    calc_phot_calib(params,reduction_metadata,vphas_cat,match_index,log)
    
    logs.close_log(log)
    
def get_args():
    """Function to gather the necessary commandline arguments"""
    
    params = { 'red_dir': '',
               'metadata': '',
               'log_dir': '',
               'catalog': '',
               'pipeline_config_dir': '',
               'software_dir': '', 
               'verbosity': '',
            }

    if len(sys.argv) > 1:
        params['red_dir'] = sys.argv[1]
        params['metadata'] = sys.argv[2]
        params['log_dir'] = sys.argv[3]
    else:
        params['red_dir'] = raw_input('Please enter the path to the reduction directory: ')
        params['metadata'] = raw_input('Please enter the name of the metadata file: ')
        params['log_dir'] = raw_input('Please enter the path to the log directory: ')
    
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
                                              'reduction_parameters' )
    reduction_metadata.load_a_layer_from_file( setup.red_dir, 
                                              params['metadata'], 
                                              'headers_summary' )
    reduction_metadata.load_a_layer_from_file( setup.red_dir, 
                                              params['metadata'],
                                              'star_catalog' )
    
    params['fov'] = reduction_metadata.reduction_parameters[1]['FOV'][0]
    params['refimage'] = reduction_metadata.data_architecture[1]['REF_IMAGE'][0]
    iref = reduction_metadata.headers_summary[1]['IMAGES'].tolist().index(params['refimage'])
    params['ra'] = reduction_metadata.headers_summary[1]['RAKEY'][iref]
    params['dec'] = reduction_metadata.headers_summary[1]['DECKEY'][iref]
    params['filter'] = reduction_metadata.headers_summary[1]['FILTKEY'][iref]
    
    log.info('Gathered information from metadata file '+params['metadata']+':')
    log.info('Image field of view: '+str(params['fov'])+'sq. deg')
    log.info('Reference image: '+params['refimage']+\
                ', index '+str(iref)+' in dataset')
    log.info('Filter used for dataset: '+params['filter'])
    log.info('Pointing center coordinates: '+params['ra']+' '+params['dec'])
    
    return reduction_metadata, params
    
def fetch_catalog_sources_within_image(params,log):
    """Function to extract the objects from the VPHAS+ catalogue within the
    field of view of the reference image, based on the metadata information."""
    
    params['radius'] = (np.sqrt(params['fov'])/2.0)*60.0
    
    log.info('Search radius: '+str(params['radius'])+' arcmin')
    
    vphas_cat = vizier_tools.search_vizier_for_sources(params['ra'], 
                                                       params['dec'], 
                                                        params['radius'], 
                                                        'VPHAS+')
        
    log.info('VPHAS+ search returned '+str(len(vphas_cat))+' entries')
    
    return vphas_cat    

def match_stars_by_position(reduction_metadata,vphas_cat,log):
    """Function to cross-match stars by position.
    
    Returns:
        :param dict match_index: { Index in vphas_cat: index in star_cat }
    """
    
    ra_stars = reduction_metadata.star_catalog[1]['RA_J2000']
    dec_stars = reduction_metadata.star_catalog[1]['DEC_J2000']
    
    stars = SkyCoord(ra_stars, dec_stars, unit="deg")
    
    # Position matching tolerance in deg
    tolerance = 1.0 / 3600.0
    
    match_index = {}
    
    for i in range(0,len(vphas_cat),1):
        
        vstar = SkyCoord(vphas_cat['_RAJ2000'][i], vphas_cat['_DEJ2000'][i], unit="deg")
        
        sep = stars.separation(vstar)
        
        idx = sep.argsort()
        
        if sep[idx[0]].value < tolerance:
            
            match_index[i] = idx[0]
            
            log.info('VPHAS star '+str(i)+' at ('+vstar.to_string()+') '+
                ' matches detected star '+str(idx[0])+' at ('+
                stars[idx[0]].to_string()+') with separation '+
                sep[idx[0]].to_string())
    
    log.info('Matched '+str(len(match_index)))
    
    return match_index

def calc_phot_calib(params,reduction_metadata,vphas_cat,match_index,log):
    """Function to plot the photometric calibration"""
    
    cat_col = params['filter'].replace('p','') + 'mag'
    
    fig = plt.figure(1)
    
    vidx = match_index.keys()
    didx = match_index.values()
    
    cat_mags = []
    instr_mags = []
    instr_mag_errs = []
    
    for vidx in match_index.keys():
        didx = match_index[vidx]
        if reduction_metadata.star_catalog[1]['ref_mag'][didx] > 0.0 and\
            reduction_metadata.star_catalog[1]['ref_mag_err'][didx] > 0.0:
            cat_mags.append(vphas_cat[cat_col][vidx])
            instr_mags.append(reduction_metadata.star_catalog[1]['ref_mag'][didx])
            instr_mag_errs.append(reduction_metadata.star_catalog[1]['ref_mag_err'][didx])
    cat_mags = np.array(cat_mags)
    instr_mags = np.array(instr_mags)
    instr_mag_errs = np.array(instr_mag_errs)
    
    fit = calc_transform(cat_mags, instr_mags)
    
    log.info('Initial fit: '+repr(fit))
    
    idx = exclude_outliers(cat_mags, instr_mags, instr_mag_errs, fit)
    
    fit = calc_transform(cat_mags[idx], instr_mags[idx])
    
    log.info('Second fit: '+repr(fit))
    
    xplot = np.linspace(cat_mags[idx].min(),cat_mags[idx].max(),50)
    yplot = phot_func(fit,xplot)
    
    plt.plot(cat_mags, instr_mags,'m.')
    
    plt.plot(xplot, yplot,'k-')
    
    plt.xlabel('VPHAS+ catalog magnitude')

    plt.ylabel('Instrumental magnitude')
    
    [xmin,xmax,ymin,ymax] = plt.axis()
    
    plt.axis([xmax,xmin,ymax,ymin])
    
    plt.savefig(os.path.join(params['red_dir'],'phot_calib.png'))

    plt.close(1)

def phot_func(p,mags):
    """Photometric function"""
        
    return p[0] + p[1]*mags
    
def calc_transform(cat_mags, instr_mags):
    """Function to calculate the photometric transformation between a set
    of catalogue magnitudes and the instrumental magnitudes for the same stars
    """
    
    errfunc = lambda p, x, y: phot_func(p,x) - y

    pinit = [ 0.0, 1.0 ]
    
    (pfit,iexec) = optimize.leastsq(errfunc,pinit,args=(cat_mags,instr_mags))

    return pfit

def exclude_outliers(cat_mags, instr_mags, instr_mag_errs, fit):
    
    pred_mags = phot_func(fit,cat_mags)
    
    idx = np.isnan(pred_mags)
        
    dif = ((instr_mags[~idx] - pred_mags[~idx])**2).sum()
    unc = (1.0/(instr_mag_errs[~idx]*instr_mag_errs[~idx])).sum()
    rms = np.sqrt( dif/unc )
    
    delta = instr_mags - pred_mags
    delta[np.isnan(delta)] = 10.0
    
    idx = np.where(abs(delta) <= 10.0*rms)
        
    return idx

if __name__ == '__main__':
    
    calibrate_photometry()