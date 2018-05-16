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
from astropy.coordinates import matching
import astropy.units as u
from astropy.table import Table
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
    
    (reduction_metadata, params, star_catalog) = fetch_metadata(setup,params,log)
    
    star_catalog = select_good_detected_stars(star_catalog,log)
    
    vphas_cat = fetch_catalog_sources_within_image(params,log)
    
    vphas_cat = select_calibration_stars(vphas_cat,params,log)
    
    match_index = match_stars_by_position(star_catalog,vphas_cat,log)
    
    fit = calc_phot_calib(params,star_catalog,vphas_cat,match_index,log)
    
    star_catalog = apply_phot_calib(star_catalog,fit,log)
    
    phot_catalog = np.zeros([len(star_catalog),2])
    phot_catalog[:,0] = star_catalog['star_index'][:]
    phot_catalog[:,1] = star_catalog['cal_ref_mag'][:]
    
    reduction_metadata.create_phot_calibration_layer(phot_catalog,log=log)
    
    reduction_metadata.save_a_layer_to_file(setup.red_dir, 
                                            params['metadata'],
                                            'phot_calib', log=log)
                                                
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
    params['cat_mag_col'] = params['filter'].replace('p','') + 'mag'
    params['cat_err_col'] = 'e_'+params['filter'].replace('p','') + 'mag'
    
    log.info('Gathered information from metadata file '+params['metadata']+':')
    log.info('Image field of view: '+str(params['fov'])+'sq. deg')
    log.info('Reference image: '+params['refimage']+\
                ', index '+str(iref)+' in dataset')
    log.info('Filter used for dataset: '+params['filter'])
    log.info('Pointing center coordinates: '+params['ra']+' '+params['dec'])
    
    star_catalog = Table()
    star_catalog['star_index'] = reduction_metadata.star_catalog[1]['star_index']
    star_catalog['RA'] = reduction_metadata.star_catalog[1]['RA_J2000']
    star_catalog['DEC'] = reduction_metadata.star_catalog[1]['DEC_J2000']
    star_catalog['mag'] = reduction_metadata.star_catalog[1]['ref_mag']
    star_catalog['mag_err'] = reduction_metadata.star_catalog[1]['ref_mag_err']
    star_catalog['clean'] = np.zeros(len(reduction_metadata.star_catalog[1]['ref_mag']))
    star_catalog['cal_ref_mag'] = np.zeros(len(reduction_metadata.star_catalog[1]['ref_mag']))
    
    log.info('Extracted star catalog')
    
    return reduction_metadata, params, star_catalog
    
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

def select_good_detected_stars(star_catalog,log):
    """Function to identify and flag stars detected in the reference image
    for which the photometry is poor."""
    
    idx1 = np.where(star_catalog['mag'] > 0.0)[0].tolist()
    idx2 = np.where(star_catalog['mag_err'] > 0.0)[0].tolist()
    
    idx = set(idx1).intersection(set(idx2))
    
    star_catalog['clean'][list(idx)] = 1.0

    log.info('Identified '+str(len(idx))+' detected stars with good photometry')
    
    return star_catalog
    
def select_calibration_stars(vphas_cat,params,log):
    """Function to remove the entries from the VPHAS catalog with poor 
    photometry. 
    Based on code by Y. Tsapras.
    
    Returns the same vphas_cat Table, but uses the clean Column to indicate
    selected stars with 1 and deselected stars with 0. 
    """
    
    cerr = params['cat_err_col']
    cmag = params['cat_mag_col']
    
    vphas_cat['clean'] = 0
    
    med = np.median(vphas_cat[cerr][np.where(vphas_cat[cerr]>0)])
    max_err = 2.0 * med
    
    log.info('Median photometric uncertainty of catalog stars: '+str(med))
    log.info('Excluding catalog stars with uncertainty > '+str(max_err))
    
    idx1 = np.where(vphas_cat[cerr] <= max_err)
    idx2 = np.where(vphas_cat[cerr] > 0)
    idx = set(idx1[0]).intersection(set(idx2[0]))
    
    limit_mag = 17.5
    if params['filter'] == 'gp':
        limit_mag = 22.0
        
    log.info('Using limiting mag '+str(limit_mag)+\
            ' for catalog selection for filter '+params['filter'])
                
    idx2 = np.where(vphas_cat[cmag] < limit_mag)
    idx = idx.intersection(set(idx2[0]))
    
    vphas_cat['clean'][list(idx)] = 1
    
    log.info('Selected '+str(len(idx))+' stars suitable for use in photometric calibration')
    
    return vphas_cat

def match_stars_by_position(star_catalog,vphas_cat,log):
    """Function to cross-match stars by position.
    
    Returns:
        :param dict match_index: { Index in vphas_cat: index in star_cat }
    """
    
    ddx = np.where(star_catalog['clean'] == 1.0)[0]
    
    det_stars = SkyCoord(star_catalog['RA'][ddx], star_catalog['DEC'][ddx], unit="deg")
    
    vdx = np.where(vphas_cat['clean'] == 1.0)[0]

    cat_stars = SkyCoord(vphas_cat['_RAJ2000'][vdx], vphas_cat['_DEJ2000'][vdx], unit="deg")
    
    tolerance = 1.0 * u.arcsec
    
    match_data = matching.search_around_sky(det_stars, cat_stars, 
                                             seplimit=tolerance)    
    
    match_index = np.array(zip(ddx[match_data[0]],vdx[match_data[1]]))
    
    log.info('Matched '+str(len(match_index)))
        
    return match_index

def calc_phot_calib(params,star_catalog,vphas_cat,match_index,log):
    """Function to plot the photometric calibration"""

    cmag = params['cat_mag_col']
    cerr = params['cat_err_col']
    
    for i in range(0,4,1):

        fit = initial_guess(params,star_catalog,vphas_cat,match_index,log)
        
        cat_mags = vphas_cat[cmag][match_index[:,1]]
        vpivot = np.median(cat_mags)
        cat_mags = cat_mags - vpivot
        
        det_mags = star_catalog['mag'][match_index[:,0]]
        dpivot = np.median(det_mags)
        det_mags = det_mags - dpivot
        
        fit = calc_transform(fit, cat_mags, det_mags)
        
        log.info('Fit result ['+str(i)+']: '+repr(fit))
        
        match_index = exclude_outliers(vphas_cat,star_catalog,params,
                                        match_index,fit,vpivot,dpivot,log)
    
    fit = fit.tolist() +  [vpivot,dpivot]
    
    fig = plt.figure(1)
    xplot = np.linspace(vphas_cat[params['cat_mag_col']][match_index[:,1]].min(),
                        vphas_cat[params['cat_mag_col']][match_index[:,1]].max(),50)
    yplot = phot_func(fit,xplot-vpivot) + dpivot
    
    plt.errorbar(vphas_cat[cmag][match_index[:,1]], 
             star_catalog['mag'][match_index[:,0]],
             yerr=star_catalog['mag_err'][match_index[:,0]],
             xerr=vphas_cat[cerr][match_index[:,1]],
             color='m', fmt='none')
    
    plt.plot(xplot, yplot,'k-')
    
    plt.xlabel('VPHAS+ catalog magnitude')

    plt.ylabel('Instrumental magnitude')
    
    [xmin,xmax,ymin,ymax] = plt.axis()
    
    plt.axis([xmax,xmin,ymax,ymin])
    
    plt.savefig(os.path.join(params['red_dir'],'phot_calib.png'))

    plt.close(1)

    log.info('Final fitted photometric calibration: '+repr(fit))

    return fit
    
def initial_guess(params,star_catalog,vphas_cat,match_index,log):
    """Function to make an initial guess at the fit parameters"""
    
    cmag = params['cat_mag_col']
    
    cat_mags = vphas_cat[cmag][match_index[:,1]]
    det_mags = star_catalog['mag'][match_index[:,0]]
    
    grad = (cat_mags[:-1] - cat_mags[1:]) / (det_mags[:-1] - det_mags[1:])
    
#    if grad.mean() < 1.0: 
    fit = [ 0.0, grad.mean() ]
#    else:
#        fit = [ 0.0, 0.5 ]
    
    log.info('Fit initial guess: '+repr(fit))
    
    return fit    
    
def phot_func(p,mags):
    """Photometric transform function"""
    
    return p[0] + p[1]*mags

def errfunc(p,x,y):
    """Function to calculate the residuals on the photometric transform"""
    
    return y - phot_func(p,x)
    
def calc_transform(pinit, x, y):
    """Function to calculate the photometric transformation between a set
    of catalogue magnitudes and the instrumental magnitudes for the same stars
    """
        
    (pfit,iexec) = optimize.leastsq(errfunc,pinit,args=(x,y))
    
    return pfit

def exclude_outliers(vphas_cat,star_catalog,params,match_index,fit,vpivot,dpivot,log):
    
    cmag = params['cat_mag_col']
    
    pred_mags = phot_func(fit, (vphas_cat[cmag][match_index[:,1]]-vpivot) )
    
    residuals = star_catalog['mag'][match_index[:,0]] - dpivot - pred_mags
    
    (median,MAD) = calc_MAD(residuals)
    
    log.info('Median, MAD of photometric residuals: '+str(median)+' '+str(MAD))
    
    jdx = np.where(residuals >= (median - 3.0*MAD))[0]
    kdx = np.where(residuals <= (median + 3.0*MAD))[0]
    dx = list(set(jdx).intersection(set(kdx)))

    log.info('Excluded '+str(len(match_index)-len(dx))+', '+\
            str(len(match_index))+' stars remaining')
    
    match_index = np.array(zip(match_index[dx,0],match_index[dx,1]))
    
    return match_index

def calc_MAD(x):
    """Function to calculate the Median Absolute Deviation from a single-column
    array of floats"""
    
    median = np.median(x)
    MAD = np.median(abs(x - np.median(x)))
    
    return median, MAD

def apply_phot_calib(star_catalog,fit_params,log):
    """Function to apply the computed photometric calibration to calculate 
    calibrated magnitudes for all detected stars"""
    
    mags = star_catalog['mag'] - fit_params[3]
    
    cal_mags = phot_func(fit_params,mags) + fit_params[2]
    
    idx = np.where(star_catalog['mag'] == 0.0)
    cal_mags[idx] = 0.0
        
    star_catalog['cal_ref_mag'] = cal_mags

    log.info('Calculated calibrated reference magnitudes for all detected stars')
    
    return star_catalog

                
if __name__ == '__main__':
    
    calibrate_photometry()