# -*- coding: utf-8 -*-
"""
Created on Tue May 15 20:48:12 2018

@author: rstreet
"""

from os import path, remove
from sys import argv
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import matching
from astropy.table import Table
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import star_colour_data
import spectral_type_data
import jester_phot_transforms
import bilir_phot_transforms
import red_clump_utilities
import interp_Bessell_Brett
import isochrone_utilities
import photometry_classes
import logging
import stellar_radius_relations
from astropy import constants

def perform_colour_analysis():
    """Function to plot colour magnitude and colour-colour plots"""
    
    tol = 2.0       # Arcmin
    calib_on_colours = False
    
    params = get_args()
    
    log = start_log(params)
    
    (star_catalog,catalog_header) = read_combined_star_catalog(params,log)
    
    target = find_target_data(params,star_catalog,log)

    (source, blend) = calc_source_blend_params(params,log)
    
    (det_idx, cat_idx, close_cat_idx) = index_valid_star_entries(star_catalog,
                                                                target,tol,log,
                                                                valid_cat=True)
    
    deltas = calibrate_instrumental_colour_colour_diagram(params,star_catalog,
                                                 catalog_header,target,
                                                 det_idx,cat_idx,close_cat_idx,
                                                 log,
                                                 calib=calib_on_colours)

    RC = localize_red_clump(star_catalog,close_cat_idx,log)
    
    analyse_colour_mag_diagrams(params,star_catalog,catalog_header,
                                target, source,blend,RC,
                                det_idx,cat_idx,close_cat_idx,log)
                                                 
    RC = measure_RC_offset(params,RC,target,log)
    
    (target,source,blend) = calc_phot_properties(target, source, blend, RC, log)
    
    
    plot_colour_colour_diagram(params,star_catalog,catalog_header,
                               target, source, blend, RC,
                               det_idx,cat_idx,close_cat_idx, log)
    
    star_data = isochrone_utilities.analyze_isochrones(source.gr_0,source.ri_0, 
                                                       params['isochrone_file'],
                                                       log=log)
    source.mass = star_data[0]
    source.sig_mass = star_data[1]
    source.teff = star_data[2]
    source.sig_teff = star_data[3]
    source.logg = star_data[4]
    source.sig_logg = star_data[5]
    
    star_data = isochrone_utilities.analyze_isochrones(blend.gr_0,blend.ri_0, 
                                                       params['isochrone_file'],
                                                       log=log)
    blend.mass = star_data[0]
    blend.sig_mass = star_data[1]
    blend.teff = star_data[2]
    blend.sig_teff = star_data[3]
    blend.logg = star_data[4]
    blend.sig_logg = star_data[5]
    
    (source, blend) = calc_source_blend_ang_radii(source, blend, log)
    
    (source, blend) = calc_source_blend_physical_radii(source, blend, log)
    
    (source,blend) = calc_source_blend_distance(source, blend, log)
    
    output_red_clump_data_latex(params,RC,log)
    
    output_source_blend_data_latex(params,source,blend,log)
    
def start_log(params, console=False):
    """Function to initialise a log file"""
    
    log_file = path.join(params['red_dir'],'colour_analysis.log')
    
    if path.isfile(log_file):
        remove(log_file)
        
    log = logging.getLogger( 'colour_analysis' )
    
    if len(log.handlers) == 0:
        log.setLevel( logging.INFO )
        file_handler = logging.FileHandler( log_file )
        file_handler.setLevel( logging.INFO )
        
        if console == True:
            console_handler = logging.StreamHandler()
            console_handler.setLevel( logging.INFO )
    
        formatter = logging.Formatter( fmt='%(asctime)s %(message)s', \
                                    datefmt='%Y-%m-%dT%H:%M:%S' )
        file_handler.setFormatter( formatter )

        if console == True:        
            console_handler.setFormatter( formatter )
    
        log.addHandler( file_handler )
        if console == True:            
            log.addHandler( console_handler )
    
    log.info('Analyzing event colour data')
    log.info('Initial parameters:')
    for key, value in params.items():
        log.info(key+': '+str(value))
        
    return log
    
def get_args():
    """Function to gather the necessary commandline arguments"""

    params = {}
    
    if len(argv) == 1:
        
        params['catalog_file'] = raw_input('Please enter the path to combined star catalog file: ')
        params['red_dir'] = raw_input('Please enter the path to the output directory: ')
        params['target_ra'] = raw_input('Please enter the RA of target [sexigesimal]:')
        params['target_dec'] = raw_input('Please enter the Dec of target [sexigesimal]:')
        params['star_class'] = raw_input('Please input the likely luminosity class of the source: ')
        params['isochrone_file'] = raw_input('Please input the path to an isochrone data file or None: ')
        params['f_s_g'] = float(raw_input('Please input the source flux (SDSS-g) or None: '))
        params['sig_f_s_g'] = float(raw_input('Please input the uncertainty for the source flux (SDSS-g) or None: '))
        params['f_b_g'] = float(raw_input('Please input the blend flux (SDSS-g) or None: '))
        params['sig_f_b_g'] = float(raw_input('Please input the uncertainty for the blend flux (SDSS-g) or None: '))
        params['f_s_r'] = float(raw_input('Please input the source flux (SDSS-r) or None: '))
        params['sig_f_s_r'] = float(raw_input('Please input the uncertainty for the source flux (SDSS-r) or None: '))
        params['f_b_r'] = float(raw_input('Please input the blend flux (SDSS-r) or None: '))
        params['sig_f_b_r'] = float(raw_input('Please input the uncertainty for the blend flux (SDSS-r) or None: '))
        params['f_s_i'] = float(raw_input('Please input the source flux (SDSS-i) or None: '))
        params['sig_f_s_i'] = float(raw_input('Please input the uncertainty for the source flux (SDSS-i) or None: '))
        params['f_b_i'] = float(raw_input('Please input the blend flux (SDSS-i) or None: '))
        params['sig_f_b_i'] = float(raw_input('Please input the uncertainty for the blend flux (SDSS-i) or None: '))
        
    else:

        params['catalog_file'] = argv[1]
        params['red_dir'] = argv[2]
        params['target_ra'] = argv[3]
        params['target_dec'] = argv[4]
        params['star_class'] = argv[5]
        params['isochrone_file'] = argv[6]
        params['f_s_g'] = float(argv[7])
        params['sig_f_s_g'] = float(argv[8])
        params['f_b_g'] = float(argv[9])
        params['sig_f_b_g'] = float(argv[10])
        params['f_s_r'] = float(argv[11])
        params['sig_f_s_r'] = float(argv[12])
        params['f_b_r'] = float(argv[13])
        params['sig_f_b_r'] = float(argv[14])
        params['f_s_i'] = float(argv[15])
        params['sig_f_s_i'] = float(argv[16])
        params['f_b_i'] = float(argv[17])
        params['sig_f_b_i'] = float(argv[18])
    
    for key, value in params.items():
        
        if 'none' in str(value).lower():
            
            value = None
        
            params[key] = value
            
    return params
    
def read_combined_star_catalog(params,log):
    """Function to read the photometric and star catalog data from a metadata file"""
    
    if path.isfile(params['catalog_file']) == False:
        
        return np.zeros(1)
    
    hdulist = fits.open(params['catalog_file'])
    
    data = hdulist[1].data
    
    header = hdulist[0].header
    
    star_catalog = Table(data)
    
    log.info('Read data from combined colour star catalog')
    
    return star_catalog, header

def calc_source_blend_params(params,log):
    """Function to construct a dictionary of needed parameters for the 
    source and blend"""
    
    source = photometry_classes.Star()
    
    source.fs_g = params['f_s_g']
    source.sig_fs_g = params['sig_f_s_g']
    (source.g, source.sig_g) = flux_to_mag_pylima(source.fs_g,source.sig_fs_g)
    
    source.fs_r = params['f_s_r']
    source.sig_fs_r = params['sig_f_s_r']
    (source.r, source.sig_r) = flux_to_mag_pylima(source.fs_r,source.sig_fs_r)
    
    source.fs_i = params['f_s_i']
    source.sig_fs_i = params['sig_f_s_i']
    (source.i, source.sig_i) = flux_to_mag_pylima(source.fs_i,source.sig_fs_i)
    
    source.compute_colours(use_inst=True)
    source.transform_to_JohnsonCousins()
    
    log.info('\n')
    log.info('Source measured photometry:')
    log.info(source.summary(show_mags=True))
    log.info(source.summary(show_mags=False,show_colours=True))
    log.info(source.summary(show_mags=False,johnsons=True))
    
    blend = photometry_classes.Star()
    
    blend.fs_g = params['f_b_g']
    blend.sig_fs_g = params['sig_f_b_g']
    (blend.g, blend.sig_g) = flux_to_mag_pylima(blend.fs_g,blend.sig_fs_g)
    
    blend.fs_r = params['f_b_r']
    blend.sig_fs_r = params['sig_f_b_r']
    (blend.r, blend.sig_r) = flux_to_mag_pylima(blend.fs_r,blend.sig_fs_r)
    
    blend.fs_i = params['f_b_i']
    blend.sig_fs_i = params['sig_f_b_i']
    (blend.i, blend.sig_i) = flux_to_mag_pylima(blend.fs_i,blend.sig_fs_i)
    
    blend.compute_colours(use_inst=True)
    blend.transform_to_JohnsonCousins()
    
    log.info('\n')
    log.info('Blend measured photometry:')
    log.info(blend.summary(show_mags=True))
    log.info(blend.summary(show_mags=False,show_colours=True))
    log.info(blend.summary(show_mags=False,johnsons=True))
    
    return source, blend
    
def flux_to_mag_pylima(flux,flux_err):
    """Function to convert the flux and flux uncertainty measured by 
    modeling in pyLIMA to magnitudes

    Uses default pyLIMA zeropoint = 27.4 mag    
    """
    
    def flux2mag(ZP, flux):
        
        return ZP - 2.5 * np.log10(flux)
    
    ZP = 27.40
    
    if flux < 0.0 or flux_err < 0.0:
        
        mag = 0.0
        mag_err = 0.0

    else:        

        mag = flux2mag(ZP, flux)
        
        mag_err = (2.5/np.log(10.0))*flux_err/flux
        
    return mag, mag_err

def find_target_data(params,star_catalog,log):
    """Function to identify the photometry for a given target, if present
    in the star catalogue"""
    
    target = photometry_classes.Star()
    
    if params['target_ra'] != None:
        
        target_location = SkyCoord([params['target_ra']], [params['target_dec']], unit=(u.hourangle, u.deg))
                
        stars = SkyCoord(star_catalog['RA'], star_catalog['DEC'], unit="deg")
        
        tolerance = 2.0 * u.arcsec
        
        match_data = matching.search_around_sky(target_location, stars, 
                                                seplimit=tolerance)    
                                                
        idx = np.argsort(match_data[2].value)
    
        if len(match_data[0]) > 0:
            target.star_index = star_catalog['star_index'][match_data[1][idx[0]]]
            target.ra = star_catalog['RA'][match_data[1][idx[0]]]
            target.dec = star_catalog['DEC'][match_data[1][idx[0]]]
            target.i = star_catalog['cal_ref_mag_ip'][match_data[1][idx[0]]]
            target.sig_i = star_catalog['cal_ref_mag_err_ip'][match_data[1][idx[0]]]
            target.r = star_catalog['cal_ref_mag_rp'][match_data[1][idx[0]]]
            target.sig_r = star_catalog['cal_ref_mag_err_rp'][match_data[1][idx[0]]]
            target.i_inst = star_catalog['ref_mag_ip'][match_data[1][idx[0]]]
            target.sig_i_inst = star_catalog['ref_mag_err_ip'][match_data[1][idx[0]]]
            target.r_inst = star_catalog['ref_mag_rp'][match_data[1][idx[0]]]
            target.sig_r_inst = star_catalog['ref_mag_err_rp'][match_data[1][idx[0]]]
            target.separation = match_data[2][idx[0]].to_string(unit=u.arcsec)
            try:
                target.g = star_catalog['cal_ref_mag_gp'][match_data[1][idx[0]]]
                target.sig_g = star_catalog['cal_ref_mag_err_gp'][match_data[1][idx[0]]]
                target.g_inst = star_catalog['ref_mag_gp'][match_data[1][idx[0]]]
                target.sig_g_inst = star_catalog['ref_mag_err_gp'][match_data[1][idx[0]]]
            except AttributeError:
                pass
            
            log.info('\n')
            log.info('Target identified as star '+str(target.star_index)+\
                        ' in the combined ROME catalog, with parameters:')
            log.info('RA = '+str(target.ra)+' Dec = '+str(target.dec))
            log.info('Measured ROME photometry, calibrated to the VPHAS+ scale:')
            log.info(target.summary(show_mags=True))
            
        if target.i != None and target.r != None:
            
            target.compute_colours(use_inst=True)
            
            log.info(target.summary(show_mags=False,show_colours=True))
            
        target.transform_to_JohnsonCousins()
        
        log.info(target.summary(show_mags=False,johnsons=True))
    
    return target


def index_valid_star_entries(star_catalog,target,tol,log,valid_cat=False):
    """Function to return an index of all stars with both full instrumental and
    catalogue entries"""
    
    idx1 = np.where(star_catalog['cal_ref_mag_ip'] > 0.0)[0]
    idx2 = np.where(star_catalog['cal_ref_mag_ip'] <= 22.0)[0]
    idx3 = np.where(star_catalog['cal_ref_mag_rp'] > 0.0)[0]
    idx4 = np.where(star_catalog['cal_ref_mag_rp'] <= 22.0)[0]
    idx5 = np.where(star_catalog['cal_ref_mag_gp'] > 0.0)[0]
    idx6 = np.where(star_catalog['cal_ref_mag_gp'] <= 22.0)[0]
    
    det_idx = set(idx1).intersection(set(idx2))
    det_idx = det_idx.intersection(set(idx3))
    det_idx = det_idx.intersection(set(idx4))
    det_idx = det_idx.intersection(set(idx5))
    det_idx = det_idx.intersection(set(idx6))
    
    log.info('Identified '+str(len(det_idx))+\
            ' detected stars with valid measurements in gri')
    
    if valid_cat == False:
        return list(det_idx), None, None
        
    idx4 = np.where(star_catalog['imag'] > 0.0)[0]
    idx5 = np.where(star_catalog['rmag'] > 0.0)[0]
    idx6 = np.where(star_catalog['gmag'] > 0.0)[0]
    
    cat_idx = det_idx.intersection(set(idx4))
    cat_idx = cat_idx.intersection(set(idx5))
    cat_idx = list(cat_idx.intersection(set(idx6)))
    det_idx = list(det_idx)
    
    log.info('Identified '+str(len(cat_idx))+\
            ' detected stars with valid catalogue entries in gri')
    
    close_idx = find_stars_close_to_target(star_catalog, target, tol, log)
    
    close_cat_idx = list(set(cat_idx).intersection(set(close_idx)))
    
    log.info('Identified '+str(len(close_cat_idx))+\
            ' stars close to the target with valid catalogue entries in gri')
            
    return det_idx, cat_idx, close_cat_idx

def find_stars_close_to_target(star_catalog, target, tol, log):
    """Function to identify stars close to the target"""
    
    tol = tol / 60.0        # Select stars within 2 arcmin of target
    det_stars = SkyCoord(star_catalog['RA'], star_catalog['DEC'], unit="deg")
    
    t = SkyCoord(target.ra, target.dec, unit="deg")
    
    seps = det_stars.separation(t)
    
    jdx = np.where(seps.deg < tol)[0]
    
    log.info('Identified '+str(len(jdx))+' stars within '+str(round(tol*60.0,1))+\
            'arcmin of the target')
    
    return jdx
    
def analyse_colour_mag_diagrams(params,star_catalog,catalog_header,
                                target,source,blend,RC,
                                det_idx,cat_idx,close_cat_idx,log):
    """Function to plot a colour-magnitude diagram"""
    
    tol = 2.0
    
    filters = { 'ip': 'SDSS-i', 'rp': 'SDSS-r', 'gp': 'SDSS-g' }
    
    inst_i = star_catalog['cal_ref_mag_ip'][det_idx]
    inst_r = star_catalog['cal_ref_mag_rp'][det_idx]
    inst_g = star_catalog['cal_ref_mag_gp'][det_idx]
    cal_i = star_catalog['imag'][cat_idx]
    cal_r = star_catalog['rmag'][cat_idx]
    cal_g = star_catalog['gmag'][cat_idx]
    inst_ri = inst_r - inst_i    # Catalogue column order is red -> blue
    inst_gr = inst_g - inst_r    
    inst_gi = inst_g - inst_i    
    cal_ri = cal_r - cal_i 
    cal_gr = cal_g - cal_r
    cal_gi = cal_g - cal_i
    
    linst_i = star_catalog['cal_ref_mag_ip'][close_cat_idx]
    linst_r = star_catalog['cal_ref_mag_rp'][close_cat_idx]
    linst_g = star_catalog['cal_ref_mag_gp'][close_cat_idx]
    lcal_i = star_catalog['imag'][close_cat_idx]
    lcal_r = star_catalog['rmag'][close_cat_idx]
    lcal_g = star_catalog['gmag'][close_cat_idx]
    linst_ri = linst_r - linst_i    # Catalogue column order is red -> blue
    linst_gr = linst_g - linst_r
    linst_gi = linst_g - linst_i
    lcal_ri = lcal_r - lcal_i
    lcal_gr = lcal_g - lcal_r
    lcal_gi = lcal_g - lcal_i
    
    plot_colour_mag_diagram(params,inst_i, inst_ri, linst_i, linst_ri, target, 
                            source, blend, RC, 'r', 'i', 'i', tol, log)
                            
    plot_colour_mag_diagram(params,inst_r, inst_ri, linst_r, linst_ri, target, 
                            source, blend, RC, 'r', 'i', 'r', tol, log)
                            
    plot_colour_mag_diagram(params,inst_g, inst_gr, linst_g, linst_gr, target, 
                            source, blend, RC, 'g', 'r', 'g', tol, log)
                            
    plot_colour_mag_diagram(params,inst_g, inst_gi, linst_g, linst_gi, target, 
                            source, blend, RC, 'g', 'i', 'g', tol, log)
    
    
def plot_colour_mag_diagram(params,mags, colours, local_mags, local_colours, 
                            target, source, blend, RC, blue_filter, red_filter, 
                            yaxis_filter, tol, log):
    """Function to plot a colour-magnitude diagram, highlighting the data for 
    local stars close to the target in a different colour from the rest, 
    and indicating the position of both the target and the Red Clump centroid.
    """
    
    fig = plt.figure(1,(10,10))
    
    plt.rcParams.update({'font.size': 18})
        
    plt.scatter(colours,mags,
                 c='#E1AE13', marker='.', s=1, label=None)
    
    plt.scatter(local_colours,local_mags,
                 c='#8c6931', marker='*', s=4, 
                 label='Stars < '+str(round(tol,1))+'arcmin of target')
    
    col_key = blue_filter+red_filter
    
    if getattr(source,blue_filter) != None and getattr(source,red_filter) != None:
        
        plt.errorbar(getattr(source,col_key), getattr(source,yaxis_filter), 
                 yerr = getattr(source,'sig_'+yaxis_filter),
                 xerr = getattr(source,'sig_'+col_key), color='m',
                 marker='d',markersize=6, label='Source')
    
    if getattr(blend,blue_filter) != None and getattr(blend,red_filter) != None:
        
        plt.errorbar(getattr(blend,col_key), getattr(blend,yaxis_filter), 
                 yerr = getattr(blend,'sig_'+yaxis_filter),
                 xerr = getattr(blend,'sig_'+col_key), color='b',
                 marker='+',markersize=6, label='Blend')
                
    if getattr(target,blue_filter) != None and getattr(target,red_filter) != None:
        
        plt.errorbar(getattr(target,col_key), getattr(target,yaxis_filter), 
                 yerr = getattr(target,'sig_'+yaxis_filter),
                 xerr = getattr(target,'sig_'+col_key), color='k',
                 marker='x',markersize=6, label='Target at baseline')
                 
    plt.errorbar(getattr(RC,col_key), getattr(RC,yaxis_filter), 
                 yerr=getattr(RC,'sig_'+yaxis_filter), 
                 xerr=getattr(RC,'sig_'+col_key),
                 color='g', marker='s',markersize=6, label='Red Clump centroid')
                 
    plt.xlabel('SDSS ('+blue_filter+'-'+red_filter+') [mag]')

    plt.ylabel('SDSS-'+yaxis_filter+' [mag]')
    
    [xmin,xmax,ymin,ymax] = plt.axis()
    
    plt.axis([xmin,xmax,ymax,ymin])
    
    plot_file = path.join(params['red_dir'],'colour_magnitude_diagram_'+\
                                            yaxis_filter+'_vs_'+blue_filter+red_filter\
                                            +'.png')

    plt.grid()
        
    if red_filter == 'i' and blue_filter == 'r' and yaxis_filter == 'i':
        plt.axis([-1.0,2.0,20.2,13.5])
        plt.legend(loc=4)  # Lower right
    
    if red_filter == 'i' and blue_filter == 'r' and yaxis_filter == 'r':
        plt.axis([-1.0,2.0,21.0,13.5])
        plt.legend(loc=4)  # Lower right
        
    if red_filter == 'r' and blue_filter == 'g':
        plt.axis([-1.0,3.0,21.0,13.5])
        plt.legend(loc=2)  # Upper left
        
    plt.savefig(plot_file)

    plt.close(1)
    
    log.info('Colour-magnitude diagram output to '+plot_file)
    
def plot_colour_colour_diagram(params,star_catalog,catalog_header,
                               target, source, blend, RC,
                               det_idx,cat_idx,close_cat_idx,log):
    """Function to plot a colour-colour diagram, if sufficient data are
    available within the given star catalog"""
    
    tol = 2.0
    
    filters = { 'ip': 'SDSS-i', 'rp': 'SDSS-r', 'gp': 'SDSS-g' }
    
    try:
    
        inst_i = star_catalog['cal_ref_mag_ip'][det_idx]
        inst_r = star_catalog['cal_ref_mag_rp'][det_idx]
        inst_g = star_catalog['cal_ref_mag_gp'][det_idx]
        inst_gr = inst_g - inst_r - RC.Egr
        inst_ri = inst_r - inst_i - RC.Eri
        
        linst_i = star_catalog['cal_ref_mag_ip'][close_cat_idx]
        linst_r = star_catalog['cal_ref_mag_rp'][close_cat_idx]
        linst_g = star_catalog['cal_ref_mag_gp'][close_cat_idx]
        lcal_i = star_catalog['imag'][close_cat_idx]
        lcal_r = star_catalog['rmag'][close_cat_idx]
        lcal_g = star_catalog['gmag'][close_cat_idx]
        linst_gr = linst_g - linst_r - RC.Egr
        linst_ri = linst_r - linst_i - RC.Eri
        
        fig = plt.figure(1,(10,10))
        
        ax = plt.axes()
        
        ax.scatter(inst_gr, inst_ri, 
                   c='#E1AE13', marker='.', s=1, label=None)
        
        ax.scatter(linst_gr, linst_ri, marker='*', s=4, c='#8c6931',
                   label='Stars < '+str(round(tol,1))+'arcmin of target')
                     
        if source.gr_0 != None and source.ri_0 != None:
            plt.plot(source.gr_0, source.ri_0,'md',markersize=6, label='Source')
            
        if blend.gr_0 != None and blend.ri_0 != None:
            plt.plot(blend.gr_0, blend.ri_0,'bp',markersize=6, label='Blend')
        
        (spectral_type, luminosity_class, gr_colour, ri_colour) = spectral_type_data.get_spectral_class_data()
        
        plot_dwarfs = True
        plot_giants = False
        for i in range(0,len(spectral_type),1):
            
            spt = spectral_type[i]+luminosity_class[i]
            
            if luminosity_class[i] == 'V':
                c = 'g'
            else:
                c = 'k'
                        
            if luminosity_class[i] == 'III' and plot_giants:
                
                plt.plot(gr_colour[i], ri_colour[i], marker='s', color=c, alpha=0.5)

                plt.annotate(spt, (gr_colour[i], ri_colour[i]-0.1), 
                                color=c, size=10, rotation=-30.0, alpha=1.0)

            if luminosity_class[i] == 'V' and plot_dwarfs:
                
                plt.plot(gr_colour[i], ri_colour[i], marker='v', color=c, alpha=0.5)

                plt.annotate(spt, (gr_colour[i], 
                               ri_colour[i]+0.1), 
                                 color=c, size=10, 
                                 rotation=-30.0, alpha=1.0)

        plt.xlabel('SDSS (g-r) [mag]')
    
        plt.ylabel('SDSS (r-i) [mag]')
        
        plot_file = path.join(params['red_dir'],'colour_colour_diagram.png')
        
        plt.axis([-1.2,2.0,-1.5,1.0])
    
        plt.grid()
        
        plt.legend(loc=4) # Lower right
    
        plt.savefig(plot_file)
    
        plt.close(1)
        
        log.info('Colour-colour diagram output to '+plot_file)
        
    except AttributeError:
            
        log.info('Warning: Insufficient data for colour-colour diagram')
        
def calibrate_instrumental_colour_colour_diagram(params,star_catalog,
                                                 catalog_header,target,
                                                 det_idx,cat_idx,close_cat_idx,
                                                 log,
                                                 calib=True):
    """Function to plot a colour-colour diagram, if sufficient data are
    available within the given star catalog"""
    
    tol = 2.0
    
    filters = { 'ip': 'SDSS-i', 'rp': 'SDSS-r', 'gp': 'SDSS-g' }
    
    deltas = {}
    deltas['di'] = 0.0
    deltas['dr'] = 0.0
    deltas['dg'] = 0.0
    deltas['dgr'] = 0.0
    deltas['dri'] = 0.0
    
    if calib:
        
        try:
        
            inst_i = star_catalog['cal_ref_mag_ip'][close_cat_idx]
            inst_r = star_catalog['cal_ref_mag_rp'][close_cat_idx]
            inst_g = star_catalog['cal_ref_mag_gp'][close_cat_idx]
            inst_gr = inst_g - inst_r    # Catalogue column order is red -> blue
            inst_ri = inst_r - inst_i   
            
            cal_i = star_catalog['imag'][close_cat_idx]
            cal_r = star_catalog['rmag'][close_cat_idx]
            cal_g = star_catalog['gmag'][close_cat_idx]
            cal_gr = cal_g - cal_r    # Catalogue column order is red -> blue
            cal_ri = cal_r - cal_i    
            
            plot_phot_transform(params, inst_i, cal_i, 'SDSS-i')
            plot_phot_transform(params, inst_r, cal_r, 'SDSS-r')
            plot_phot_transform(params, inst_g, cal_g, 'SDSS-g')
            
            deltas['di'] = np.median(cal_i - inst_i)
            deltas['dr'] = np.median(cal_r - inst_r)
            deltas['dg'] = np.median(cal_g - inst_g)
            deltas['dgr'] = np.median(cal_gr - inst_gr)
            deltas['dri'] = np.median(cal_ri - inst_ri)
            
            if len(target) > 0 and target['cal_ref_mag_ip'] > 0.0 and \
                target['cal_ref_mag_rp'] > 0.0 and target['cal_ref_mag_gp']:
                target_gr = target['cal_ref_mag_gp'] - target['cal_ref_mag_rp']
                target_ri = target['cal_ref_mag_rp'] - target['cal_ref_mag_ip']
        
            fig = plt.figure(1, (10,10))
            
            plt.rcParams.update({'font.size': 18})
            
            plt.plot(inst_gr, inst_ri,'+',
                     color='#2b8c85',markersize=4,alpha=0.4,
                     label='Instrumental')
            plt.plot(cal_gr, cal_ri,'.',
                     color='#8c6931',markersize=4,alpha=0.4,
                     label='Catalogue')
            
    #        plt.plot(inst_gr+deltas['dgr'], inst_ri+deltas['dri'],'.',
    #                 color='m',markersize=2,alpha=0.4,
    #                 label='Calibrated')
                     
    #        if len(target) > 0 and target['mag1'] > 0.0 and \
    #            target['mag2'] > 0.0 and target['mag2']:
    #            plt.plot(target_colour2, target_colour1,'md',markersize=6)
                
            plt.xlabel('SDSS (g-r) [mag]')
        
            plt.ylabel('SDSS (r-i) [mag]')
            
            plot_file = path.join(params['red_dir'],'calib_colour_colour_diagram.png')
            
            plt.grid()
            
            plt.axis([-1.0,5.0,0.0,2.0])
            plt.legend()
            
            plt.savefig(plot_file)
        
            plt.close(1)
            
            log.info('Calibration colour-colour diagram output to '+plot_file)
            
            log.info('Measured offsets in photometry:')
            
            for key, value in deltas.items():
                log.info(key+' = '+str(value))
            
        except AttributeError:
            
            deltas = {}
            log.info('Warning: Insufficient data for colour-colour diagram')
        
    return deltas

def plot_phot_transform(params, inst_mag, cal_mag, bandpass):
    """Function to plot the relationship between the catalogue and instrumental
    magnitudes of a single passband"""
    
    fig = plt.figure(2)

    plt.plot(cal_mag, inst_mag,'k.')

    plt.xlabel('Catalog magnitude')

    plt.ylabel('Instrumental magnitude')
    
    plt.title('Relation between instrumental and catalogue magnitudes in '+\
                bandpass)
                
    [xmin,xmax,ymin,ymax] = plt.axis()
    
    plt.axis([xmax,xmin,ymax,ymin])
    
    plt.savefig(path.join(params['red_dir'],
                'phot_transform_'+bandpass+'.png'))

    plt.close(2)

def localize_red_clump(star_catalog,close_cat_idx,log):
    """Function to calculate the centroid of the Red Clump stars in a 
    colour-magnitude diagram"""
    
    def select_within_range(mags, colours, mag_min, mag_max, col_min, col_max):
        """Function to identify the set of array indices with values
        between the range indicated"""
        
        idx1 = np.where(colours >= col_min)[0]
        idx2 = np.where(colours <= col_max)[0]
        idx3 = np.where(mags >= mag_min)[0]
        idx4 = np.where(mags <= mag_max)[0]
        idx = set(idx1).intersection(set(idx2))
        idx = idx.intersection(set(idx3))
        idx = list(idx.intersection(set(idx4)))
        
        return idx
    
    RC = photometry_classes.Star()
    
    inst_i = star_catalog['cal_ref_mag_ip'][close_cat_idx]
    inst_r = star_catalog['cal_ref_mag_rp'][close_cat_idx]
    inst_g = star_catalog['cal_ref_mag_gp'][close_cat_idx]
    cal_i = star_catalog['imag'][close_cat_idx]
    cal_r = star_catalog['rmag'][close_cat_idx]
    cal_g = star_catalog['gmag'][close_cat_idx]
    inst_ri = inst_r - inst_i    # Catalogue column order is red -> blue
    inst_gi = inst_g - inst_i 
    inst_gr = inst_g - inst_r  
    cal_ri = cal_r - cal_i 
    cal_gi = cal_g - cal_i
    cal_gr = cal_g - cal_r
    
    log.info('\n')
    log.info('Localizing the Red Clump')
    log.info('Median (r-i), i: '+str(np.median(inst_ri))+', '+str(np.median(inst_i)))
    log.info('Median (g-i), i: '+str(np.median(inst_gi))+', '+str(np.median(inst_i)))
    log.info('Median (g-r), g: '+str(np.median(inst_gr))+', '+str(np.median(inst_g)))
    
    ri_min = 0.8 
    ri_max = 1.1 
    i_min = 15.5
    i_max = 16.5
    
    r_min = 16.2
    r_max = 17.5
    
    gi_min = 2.0 
    gi_max = 4.0
    
    gr_min = 1.5 
    gr_max = 2.2 
    g_min = 18.5
    g_max = 19.5
    
    log.info('Selected Red Clump giants between:')
    log.info('i = '+str(i_min)+' to '+str(i_max))
    log.info('r = '+str(r_min)+' to '+str(r_max))
    log.info('(r-i) = '+str(ri_min)+' to '+str(ri_max))
    log.info('g = '+str(g_min)+' to '+str(g_max))
    log.info('(g-r) = '+str(gr_min)+' to '+str(gr_max))
    log.info('(g-i) = '+str(gi_min)+' to '+str(gi_max))
    
    idx = select_within_range(inst_i, inst_ri, i_min, i_max, ri_min, ri_max)
    
    (RC.ri, RC.sig_ri, RC.i, RC.sig_i) = calc_distribution_centroid_and_spread_2d(inst_ri[idx], inst_i[idx], use_iqr=True)
    
    idx = select_within_range(inst_r, inst_ri, r_min, r_max, ri_min, ri_max)
    
    (RC.r, RC.sig_r) = calc_distribution_centre_and_spread(inst_r[idx], use_iqr=True)
    
    idx = select_within_range(inst_g, inst_gr, g_min, g_max, gr_min, gr_max)
    
    (RC.gr, RC.sig_gr, RC.g, RC.sig_g) = calc_distribution_centroid_and_spread_2d(inst_gr[idx], inst_g[idx], use_iqr=True)
    
    idx = select_within_range(inst_g, inst_gi, g_min, g_max, gi_min, gi_max)
    
    (RC.gi, RC.sig_gi, RC.g, RC.sig_g) = calc_distribution_centroid_and_spread_2d(inst_gi[idx], inst_g[idx], use_iqr=True)
    
    log.info('\n')
    log.info('Centroid of Red Clump Stars at:')
    log.info(RC.summary(show_mags=True))
    log.info(RC.summary(show_mags=False,show_colours=True))
    
    RC.transform_to_JohnsonCousins()
    
    log.info(RC.summary(show_mags=False,johnsons=True))
    
    return RC

def calc_distribution_centroid_and_spread_2d(xdata, ydata, use_iqr=False):
    """Function to calculate the centroid of a 2D distribution and
    estimate the uncertainty on those values by different statistics"""
    
    (xcentre, sig_x) = calc_distribution_centre_and_spread(xdata, use_iqr=use_iqr)
    (ycentre, sig_y) = calc_distribution_centre_and_spread(ydata, use_iqr=use_iqr)
    
    return xcentre, sig_x, ycentre, sig_y

def calc_distribution_centre_and_spread(xdata, use_iqr=False):
    """Function to calculate the centroid of a 1D distribution and
    estimate the uncertainty on those values by different statistics"""
    
    xcentre = np.median(xdata)
    
    xmad = np.median(abs(xdata - xcentre))
    
    xiq_min = np.percentile(xdata,25.0)
    xiq_max = np.percentile(xdata,75.0)
    xiqr = (xiq_max - xiq_min)/2.0
    
    if use_iqr:
        sig_x = xiqr
    else:
        sig_x = xmad
    
    return xcentre, sig_x
    
def measure_RC_offset(params,RC,target,log):
    """Function to calculate the offset of the Red Clump from its expected 
    values, taken from Bensby et al. (2017), 2017, A&A, 605A, 89 for V, I bands and
    Hawkins et al. (2017) MNRAS, 471, 722 for 2MASS J,H,Ks.
    """
    
    in_use = False
    use_2mass = False
    
    RC = red_clump_utilities.get_essential_parameters(RC=RC)
    
    log.info('\n Red Clump colours and absolute SDSS magnitudes:')
    log.info('Mg_RC,0 = '+str(RC.M_g_0)+' +/- '+str(RC.sig_Mg_0)+'mag')
    log.info('Mr_RC,0 = '+str(RC.M_r_0)+' +/- '+str(RC.sig_Mr_0)+'mag')
    log.info('Mi_RC,0 = '+str(RC.M_i_0)+' +/- '+str(RC.sig_Mi_0)+'mag')
    log.info('MI_RC,0 = '+str(RC.M_I_0)+' +/- '+str(RC.sig_MI_0)+'mag')
    log.info('MV_RC,0 = '+str(RC.M_V_0)+' +/- '+str(RC.sig_MV_0)+'mag')
    log.info('(g-r)_RC,0 = '+str(RC.gr_0)+' +/- '+str(RC.sig_gr_0)+'mag')
    log.info('(g-i)_RC,0 = '+str(RC.gi_0)+' +/- '+str(RC.sig_gi_0)+'mag')
    log.info('(r-i)_RC,0 = '+str(RC.ri_0)+' +/- '+str(RC.sig_ri_0)+'mag')
    log.info('(V-I)_RC,0 = '+str(RC.VI_0)+' +/- '+str(RC.sig_VI_0)+'mag')
    
    if use_2mass:
        RC.transform_2MASS_to_SDSS()
        
        log.info('\n Red Clump NIR colours and magnitudes:')
        log.info('J_RC,0 = '+str(RC.M_J_0)+' +/- '+str(RC.sig_MJ_0)+'mag')
        log.info('H_RC,0 = '+str(RC.M_H_0)+' +/- '+str(RC.sig_MH_0)+'mag')
        log.info('Ks_RC,0 = '+str(RC.M_Ks_0)+' +/- '+str(RC.sig_MKs_0)+'mag')
        log.info('(J-H)_RC,0 = '+str(RC.JH_0)+' +/- '+str(RC.sig_JH_0)+'mag')
        log.info('(H-Ks)_RC,0 = '+str(RC.HK_0)+' +/- '+str(RC.sig_HK_0)+'mag')
    
    RC.D = red_clump_utilities.calc_red_clump_distance(params['target_ra'],params['target_dec'],log=log)
    RC = red_clump_utilities.calc_apparent_magnitudes(RC)
    
    log.info('\n Red Clump apparent SDSS magnitudes at distance '+str(RC.D)+'Kpc')
    log.info('g_RC,app = '+str(RC.m_g_0)+' +/- '+str(RC.sig_mg_0)+'mag')
    log.info('r_RC,app = '+str(RC.m_r_0)+' +/- '+str(RC.sig_mr_0)+'mag')
    log.info('i_RC,app = '+str(RC.m_i_0)+' +/- '+str(RC.sig_mi_0)+'mag')
    log.info('V_RC,app = '+str(RC.m_V_0)+' +/- '+str(RC.sig_mV_0)+'mag')
    log.info('I_RC,app = '+str(RC.m_I_0)+' +/- '+str(RC.sig_mI_0)+'mag')
    
    if in_use:
        RC.transform_to_JohnsonCousins()
        
        log.info('\n Derived Red Clump instrumental colours and magnitudes:')
        log.info(RC.summary(show_mags=False,johnsons=True))
    
    RC.A_g = RC.g - RC.m_g_0
    RC.sig_A_g = np.sqrt(RC.sig_mg_0*RC.sig_mg_0)
    RC.A_r = RC.r - RC.m_r_0
    RC.sig_A_r = np.sqrt(RC.sig_mr_0*RC.sig_mr_0)
    RC.A_i = RC.i - RC.m_i_0
    RC.sig_A_i = np.sqrt(RC.sig_mi_0*RC.sig_mi_0)
    
    RC.A_I = RC.I - RC.m_I_0
    RC.sig_A_I = RC.sig_mI_0
    RC.A_V = RC.V - RC.m_V_0
    RC.sig_A_V = RC.sig_mV_0
        
    RC.Egr = RC.gr - RC.gr_0
    RC.sig_Egr = np.sqrt( (RC.sig_gr_0*RC.sig_gr_0) )
    RC.Egi = RC.gi - RC.gi_0
    RC.sig_Egi = np.sqrt( (RC.sig_gi_0*RC.sig_gi_0) )
    RC.Eri = RC.ri - RC.ri_0
    RC.sig_Eri = np.sqrt( (RC.sig_ri_0*RC.sig_ri_0) )

    RC.EVI = RC.VI - RC.VI_0
    RC.sig_EVI = RC.sig_VI_0
    
    log.info('\n')
    log.info('Extinction, d(g) = '+str(RC.A_g)+' +/- '+str(RC.sig_A_g)+'mag')
    log.info('Extinction, d(r) = '+str(RC.A_r)+' +/- '+str(RC.sig_A_r)+'mag')
    log.info('Extinction, d(i) = '+str(RC.A_i)+' +/- '+str(RC.sig_A_i)+'mag')
    log.info('Reddening, E(g-r) = '+str(RC.Egr)+' +/- '+str(RC.sig_Egr)+'mag')
    log.info('Reddening, E(g-i) = '+str(RC.Egi)+' +/- '+str(RC.sig_Egi)+'mag')
    log.info('Reddening, E(r-i) = '+str(RC.Eri)+' +/- '+str(RC.sig_Eri)+'mag')
    
    log.info('\n')
    log.info('Extinction, d(V) = '+str(RC.A_V)+' +/- '+str(RC.sig_A_V)+'mag')
    log.info('Extinction, d(I) = '+str(RC.A_I)+' +/- '+str(RC.sig_A_I)+'mag')
    log.info('Reddening, E(V-I) = '+str(RC.EVI)+' +/- '+str(RC.sig_EVI)+'mag')
    
    return RC
    
def calc_phot_properties(target, source, blend, RC, log):
    """Function to calculate the de-reddened and extinction-corrected 
    photometric properties of the target
    """
    in_use = False
    
    target.calibrate_phot_properties(RC,log=log)
    source.calibrate_phot_properties(RC,log=log)
    blend.calibrate_phot_properties(RC,log=log)

    log.info('\nSource star extinction-corrected magnitudes and de-reddened colours:\n')
    log.info(source.summary(show_mags=False,show_cal=True))
    log.info(source.summary(show_mags=False,show_cal=True,show_colours=True))
    log.info(source.summary(show_mags=False,johnsons=True,show_cal=True))
    
    log.info('\nBlend extinction-corrected magnitudes and de-reddened colours:\n')
    log.info(blend.summary(show_mags=False,show_cal=True))
    log.info(blend.summary(show_mags=False,show_cal=True,show_colours=True))
    log.info(blend.summary(show_mags=False,johnsons=True,show_cal=True))
    
    return target,source,blend

def calc_source_blend_ang_radii(source, blend, log):
    """Function to calculate the angular radius of the source star"""
    
    log.info('\n')
    log.info('Calculating the angular radius of the source star:')
    source.calc_stellar_ang_radius(log)
    log.info('Source angular radius (from SDSS (g-i), Boyajian+ 2014 relations) = '+str(source.ang_radius)+' '+str(source.sig_ang_radius))
    
    log.info('\n')
    log.info('Calculating the angular radius of the blend:')
    blend.calc_stellar_ang_radius(log)
    log.info('Blend angular radius (from SDSS (g-i), Boyajian+ 2014 relations) = '+str(blend.ang_radius)+' '+str(blend.sig_ang_radius))
    
    return source, blend
    
def calc_source_blend_physical_radii(source, blend, log):
    """Function to infer the physical radius of the source star from the 
    Torres mass-radius relation based on Teff, logg, and Fe/H
    
    Assumes a solar metallicity of Zsol = 0.0152.    
    """
    
    source.calc_physical_radius(log)
                    
    blend.calc_physical_radius(log)
    
    log.info('\n')
    log.info('Source radius from Torres relation: '+\
                    str(round(source.radius,2))+' +/- '+str(round(source.sig_radius,2))+' Rsol')
    log.info('Blend radius from Torres relation: '+\
                    str(round(blend.radius,2))+' +/- '+str(round(blend.sig_radius,2))+' Rsol')
                    
    return source, blend
    
def convert_ndp(value,ndp):
    """Function to convert a given floating point value to a string, 
    rounded to the given number of decimal places, and suffix with zero
    if the value rounds to fewer decimal places than expected"""
    
    value = str(round(value,ndp))
    
    dp = value.split('.')[-1]
    
    while len(dp) < ndp:
        
        value = value + '0'
        
        dp = value.split('.')[-1]
    
    return value

def calc_source_blend_distance(source,blend,log):
    """Function to calculate the distance to the source star, given the
    angular and physical radius estimates"""
    
    
    log.info('\n')
    
    source.calc_distance(log)
    
    log.info('Inferred source distance: '+str(source.D)+' +/- '+str(source.sig_D)+' pc')
    
    
    blend.calc_distance(log)
    
    log.info('Inferred blend distance: '+str(blend.D)+' +/- '+str(blend.sig_D)+' pc')
    
    return source, blend

def output_red_clump_data_latex(params,RC,log):
    """Function to output a LaTeX format table with the data for the Red Clump"""
    
    file_path = path.join(params['red_dir'],'red_clump_data_table.tex')
    
    t = open(file_path, 'w')
    
    t.write('\\begin{table}[h!]\n')
    t.write('\\centering\n')
    t.write('\\caption{Photometric properties of the Red Clump, with absolute magnitudes (M) taken from \cite{Ruiz-Dern2018}, and the measured properties from ROME data.} \label{tab:RCproperties}\n')
    t.write('\\begin{tabular}{ll}\n')
    t.write('\\hline\n')
    t.write('\\hline\n')
    t.write('$M_{g,RC,0}$ & '+convert_ndp(RC.M_g_0,3)+' $\pm$ '+convert_ndp(RC.sig_Mg_0,3)+'\,mag\\\\\n')
    t.write('$M_{r,RC,0}$ & '+convert_ndp(RC.M_r_0,3)+' $\pm$ '+convert_ndp(RC.sig_Mr_0,3)+'\,mag\\\\\n')
    t.write('$M_{i,RC,0}$ & '+convert_ndp(RC.M_i_0,3)+' $\pm$ '+convert_ndp(RC.sig_Mi_0,3)+'\,mag\\\\\n')
    t.write('$(g-r)_{RC,0}$ & '+convert_ndp(RC.gr_0,3)+' $\pm$ '+convert_ndp(RC.sig_gr_0,3)+'\,mag\\\\\n')
    t.write('$(g-i)_{RC,0}$ & '+convert_ndp(RC.gi_0,3)+' $\pm$ '+convert_ndp(RC.sig_gi_0,3)+'\,mag\\\\\n')
    t.write('$(r-i)_{RC,0}$ & '+convert_ndp(RC.ri_0,3)+' $\pm$ '+convert_ndp(RC.sig_ri_0,3)+'\,mag\\\\\n')
    t.write('$m_{g,RC,0}$ & '+convert_ndp(RC.m_g_0,3)+' $\pm$ '+convert_ndp(RC.sig_mg_0,3)+'\,mag\\\\\n')
    t.write('$m_{r,RC,0}$ & '+convert_ndp(RC.m_r_0,3)+' $\pm$ '+convert_ndp(RC.sig_mr_0,3)+'\,mag\\\\\n')
    t.write('$m_{i,RC,0}$ & '+convert_ndp(RC.m_i_0,3)+' $\pm$ '+convert_ndp(RC.sig_mi_0,3)+'\,mag\\\\\n')
    t.write('$m_{g,RC,\\rm{centroid}}$  & '+convert_ndp(RC.g,2)+' $\pm$ '+convert_ndp(RC.sig_g,2)+'\,mag\\\\\n')
    t.write('$m_{r,RC,\\rm{centroid}}$  & '+convert_ndp(RC.r,2)+' $\pm$ '+convert_ndp(RC.sig_r,2)+'\,mag\\\\\n')
    t.write('$m_{i,RC,\\rm{centroid}}$  & '+convert_ndp(RC.i,2)+' $\pm$ '+convert_ndp(RC.sig_i,2)+'\,mag\\\\\n')
    t.write('$(g-r)_{RC,\\rm{centroid}}$ & '+convert_ndp(RC.gr,2)+'  $\pm$ '+convert_ndp(RC.sig_gr,2)+'\,mag\\\\\n')
    t.write('$(r-i)_{RC,\\rm{centroid}}$ & '+convert_ndp(RC.ri,2)+' $\pm$ '+convert_ndp(RC.sig_ri,2)+'\,mag\\\\\n')
    t.write('$A_{g}$ & '+convert_ndp(RC.A_g,3)+' $\pm$ '+convert_ndp(RC.sig_A_g,3)+'\,mag\\\\\n')
    t.write('$A_{r}$ & '+convert_ndp(RC.A_r,3)+' $\pm$ '+convert_ndp(RC.sig_A_r,3)+'\,mag\\\\\n')
    t.write('$A_{i}$ & '+convert_ndp(RC.A_i,3)+' $\pm$ '+convert_ndp(RC.sig_A_i,3)+'\,mag\\\\\n')
    t.write('$E(g-r)$ & '+convert_ndp(RC.Egr,3)+' $\pm$ '+convert_ndp(RC.sig_Egr,3)+'\,mag\\\\\n')
    t.write('$E(r-i)$ & '+convert_ndp(RC.Eri,3)+' $\pm$ '+convert_ndp(RC.sig_Eri,3)+'\,mag\\\\\n')
    t.write('\\hline\n')
    t.write('\\end{tabular}\n')
    t.write('\\end{table}\n')

    t.close()
    
    log.info('\n')
    log.info('Output red clump data in laTex table to '+file_path)
    
def output_source_blend_data_latex(params,source,blend,log):
    """Function to output a LaTex format table with the source and blend data"""

    file_path = path.join(params['red_dir'],'source_blend_data_table.tex')
    
    t = open(file_path, 'w')

    t.write('\\begin{table}[h!]\n')
    t.write('\\centering\n')
    t.write('\\caption{Photometric properties of the source star (s) and blend (b).} \label{tab:targetphot}\n')
    t.write('\\begin{tabular}{llll}\n')
    t.write('\\hline\n')
    t.write('\\hline\n')
    t.write('$m_{g,s}$ & '+convert_ndp(source.g,3)+' $\pm$ '+convert_ndp(source.sig_g,3)+'\,mag & $m_{g,b}$ & '+convert_ndp(blend.g,3)+' $\pm$ '+convert_ndp(blend.sig_g,3)+'\,mag\\\\\n')
    t.write('$m_{r,s}$ & '+convert_ndp(source.r,3)+' $\pm$ '+convert_ndp(source.sig_r,3)+'\,mag & $m_{r,b}$ & '+convert_ndp(blend.r,3)+' $\pm$ '+convert_ndp(blend.sig_r,3)+'\,mag\\\\\n')
    t.write('$m_{i,s}$ & '+convert_ndp(source.i,3)+' $\pm$ '+convert_ndp(source.sig_i,3)+'\,mag & $m_{i,b}$ & '+convert_ndp(blend.i,3)+' $\pm$ '+convert_ndp(blend.sig_i,3)+'\,mag\\\\\n')
    t.write('$(g-r)_{s}$ & '+convert_ndp(source.gr,3)+' $\pm$ '+convert_ndp(source.sig_gr,3)+'\,mag & $(g-r)_{b}$ & '+convert_ndp(blend.gr,3)+' $\pm$ '+convert_ndp(blend.sig_gr,3)+'\,mag\\\\\n')
    t.write('$(r-i)_{s}$ & '+convert_ndp(source.ri,3)+' $\pm$ '+convert_ndp(source.sig_ri,3)+'\,mag & $(r-i)_{b}$ & '+convert_ndp(blend.ri,3)+' $\pm$ '+convert_ndp(blend.sig_ri,3)+'\,mag\\\\\n')
    t.write('$m_{g,s,0}$ & '+convert_ndp(source.g_0,3)+' $\pm$ '+convert_ndp(source.sig_g_0,3)+'\,mag & $m_{g,b,0}$ & '+convert_ndp(blend.g_0,3)+' $\pm$ '+convert_ndp(blend.sig_g_0,3)+'\,mag\\\\\n')
    t.write('$m_{r,s,0}$ & '+convert_ndp(source.r_0,3)+' $\pm$ '+convert_ndp(source.sig_r_0,3)+'\,mag & $m_{r,b,0}$ & '+convert_ndp(blend.r_0,3)+' $\pm$ '+convert_ndp(blend.sig_r_0,3)+'\,mag\\\\\n')
    t.write('$m_{i,s,0}$ & '+convert_ndp(source.i_0,3)+' $\pm$ '+convert_ndp(source.sig_i_0,3)+'\,mag & $m_{i,b,0}$ & '+convert_ndp(blend.i_0,3)+' $\pm$ '+convert_ndp(blend.sig_i_0,3)+'\,mag\\\\\n')
    t.write('$(g-r)_{s,0}$ & '+convert_ndp(source.gr_0,3)+' $\pm$ '+convert_ndp(source.sig_gr_0,3)+'\,mag & $(g-r)_{b,0}$ & '+convert_ndp(blend.gr_0,3)+' $\pm$ '+convert_ndp(blend.sig_gr_0,3)+'\,mag\\\\\n')
    t.write('$(r-i)_{s,0}$ & '+convert_ndp(source.ri_0,3)+' $\pm$ '+convert_ndp(source.sig_ri_0,3)+'\,mag & $(r-i)_{b,0}$ & '+convert_ndp(blend.ri_0,3)+' $\pm$ '+convert_ndp(blend.sig_ri_0,3)+'\,mag\\\\\n')
    t.write('\\hline\n')
    t.write('\\end{tabular}\n')
    t.write('\\end{table}\n')

    t.close()

    log.info('Output source and blend data in laTex table to '+file_path)
    
if __name__ == '__main__':
    
    perform_colour_analysis()
    