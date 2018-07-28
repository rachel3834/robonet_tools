# -*- coding: utf-8 -*-
"""
Created on Tue May 15 20:48:12 2018

@author: rstreet
"""

from os import path
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

def perform_colour_analysis():
    """Function to plot colour magnitude and colour-colour plots"""
    
    tol = 2.0       # Arcmin
    calib_on_colours = False
    
    params = get_args()
    
    (star_catalog,catalog_header) = read_combined_star_catalog(params)
    
    target = find_target_data(params,star_catalog)

    (source, blend) = calc_source_blend_params(params)
    
    (det_idx, cat_idx, close_cat_idx) = index_valid_star_entries(star_catalog,
                                                                target,tol,
                                                                valid_cat=True)
    
    deltas = calibrate_instrumental_colour_colour_diagram(params,star_catalog,
                                                 catalog_header,target,
                                                 det_idx,cat_idx,close_cat_idx,
                                                 calib=calib_on_colours)

    RC = localize_red_clump(star_catalog,close_cat_idx)
    
    analyse_colour_mag_diagrams(params,star_catalog,catalog_header,
                                source,blend,RC,
                                det_idx,cat_idx,close_cat_idx)
                                                 
    RC = measure_RC_offset(params,RC,target)
    
    (target,source,blend) = calc_phot_properties(target, source, blend, RC)
    
    
    plot_colour_colour_diagram(params,star_catalog,catalog_header,target,RC,
                                                 det_idx,cat_idx,close_cat_idx)

    
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
        params['f_b_g'] = float(raw_input('Please input the blend flux (SDSS-g) or None: '))
        params['f_s_r'] = float(raw_input('Please input the source flux (SDSS-r) or None: '))
        params['f_b_r'] = float(raw_input('Please input the blend flux (SDSS-r) or None: '))
        params['f_s_i'] = float(raw_input('Please input the source flux (SDSS-i) or None: '))
        params['f_b_i'] = float(raw_input('Please input the blend flux (SDSS-i) or None: '))
        
    else:

        params['catalog_file'] = argv[1]
        params['red_dir'] = argv[2]
        params['target_ra'] = argv[3]
        params['target_dec'] = argv[4]
        params['star_class'] = argv[5]
        params['isochrone_file'] = argv[6]
        params['f_s_g'] = float(argv[7])
        params['f_b_g'] = float(argv[8])
        params['f_s_r'] = float(argv[9])
        params['f_b_r'] = float(argv[10])
        params['f_s_i'] = float(argv[11])
        params['f_b_i'] = float(argv[12])
    
    for key, value in params.items():
        
        if 'none' in str(value).lower():
            
            value = None
        
            params[key] = value
            
    return params
    
def read_combined_star_catalog(params):
    """Function to read the photometric and star catalog data from a metadata file"""
    
    if path.isfile(params['catalog_file']) == False:
        
        return np.zeros(1)
    
    hdulist = fits.open(params['catalog_file'])
    
    data = hdulist[1].data
    
    header = hdulist[0].header
    
    star_catalog = Table(data)
    
    return star_catalog, header

def calc_source_blend_params(params):
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
    
    print('Source measured photometry:')
    print(target.summary(show_mags=True))
    print(target.summary(show_colours=True))
    print(target.summary(johnsons=True))
    
    blend = photometry_classes.Star()
    
    blend.fs_g = params['f_s_g']
    blend.sig_fs_g = params['sig_f_s_g']
    (blend.g, blend.sig_g) = flux_to_mag_pylima(blend.fs_g,blend.sig_fs_g)
    
    blend.fs_r = params['f_s_r']
    blend.sig_fs_r = params['sig_f_s_r']
    (blend.r, blend.sig_r) = flux_to_mag_pylima(blend.fs_r,blend.sig_fs_r)
    
    blend.fs_i = params['f_s_i']
    blend.sig_fs_i = params['sig_f_s_i']
    (blend.i, blend.sig_i) = flux_to_mag_pylima(blend.fs_i,blend.sig_fs_i)
    
    blend.compute_colours(use_inst=True)
    blend.transform_to_JohnsonCousins()
    
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

def find_target_data(params,star_catalog):
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
            
            print('Target measured photometry:')
            print(target.summary(show_mags=True))
            
        if target.i != None and target.r != None:
            
            target.compute_colours(use_inst=True)
            
            print(target.summary(show_colours=True))
            
        target.transform_to_JohnsonCousins()
        
        print(target.summary(johnsons=True))
    
    return target


def index_valid_star_entries(star_catalog,target,tol,valid_cat=False):
    """Function to return an index of all stars with both full instrumental and
    catalogue entries"""
    
    idx1 = np.where(star_catalog['cal_ref_mag_ip'] > 0.0)[0]
    idx2 = np.where(star_catalog['cal_ref_mag_rp'] > 0.0)[0]
    idx3 = np.where(star_catalog['cal_ref_mag_gp'] > 0.0)[0]
    
    det_idx = set(idx1).intersection(set(idx2))
    det_idx = det_idx.intersection(set(idx3))
    
    print 'Identified '+str(len(det_idx))+' detected stars with valid measurements in gri'
    
    if valid_cat == False:
        return list(det_idx), None, None
        
    idx4 = np.where(star_catalog['imag'] > 0.0)[0]
    idx5 = np.where(star_catalog['rmag'] > 0.0)[0]
    idx6 = np.where(star_catalog['gmag'] > 0.0)[0]
    
    cat_idx = det_idx.intersection(set(idx4))
    cat_idx = cat_idx.intersection(set(idx5))
    cat_idx = list(cat_idx.intersection(set(idx6)))
    det_idx = list(det_idx)
    
    print 'Identified '+str(len(cat_idx))+' detected stars with valid catalogue entries in gri'
    
    close_idx = find_stars_close_to_target(star_catalog, target, tol)
    
    close_cat_idx = list(set(cat_idx).intersection(set(close_idx)))
    
    print 'Identified '+str(len(close_cat_idx))+\
            ' stars close to the target with valid catalogue entries in gri'
            
    return det_idx, cat_idx, close_cat_idx

def find_stars_close_to_target(star_catalog, target, tol):
    """Function to identify stars close to the target"""
    
    tol = tol / 60.0        # Select stars within 2 arcmin of target
    det_stars = SkyCoord(star_catalog['RA'], star_catalog['DEC'], unit="deg")
    
    t = SkyCoord(target['ra'], target['dec'], unit="deg")
    
    seps = det_stars.separation(t)
    
    jdx = np.where(seps.deg < tol)[0]
    
    print 'Identified '+str(len(jdx))+' stars within '+str(round(tol*60.0,1))+\
            'arcmin of the target'
    
    return jdx
    
def analyse_colour_mag_diagrams(params,star_catalog,catalog_header,
                                source,blend,RC,
                                det_idx,cat_idx,close_cat_idx):
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
    cal_ri = cal_r - cal_i 
    cal_gr = cal_g - cal_r
    
    linst_i = star_catalog['cal_ref_mag_ip'][close_cat_idx]
    linst_r = star_catalog['cal_ref_mag_rp'][close_cat_idx]
    linst_g = star_catalog['cal_ref_mag_gp'][close_cat_idx]
    lcal_i = star_catalog['imag'][close_cat_idx]
    lcal_r = star_catalog['rmag'][close_cat_idx]
    lcal_g = star_catalog['gmag'][close_cat_idx]
    linst_ri = linst_r - linst_i    # Catalogue column order is red -> blue
    linst_gr = linst_g - linst_r
    lcal_ri = lcal_r - lcal_i
    lcal_gr = lcal_g - lcal_r
    
    plot_colour_mag_diagram(params,inst_i, inst_ri, linst_i, linst_ri, 
                            source, blend, RC, 'r', 'i', 'i', tol)
                            
    plot_colour_mag_diagram(params,inst_r, inst_ri, linst_r, linst_ri, 
                            source, blend, RC, 'r', 'i', 'r', tol)
                            
    plot_colour_mag_diagram(params,inst_g, inst_gr, linst_g, linst_gr, 
                            source, blend, RC, 'g', 'r', 'g', tol)
    
    
def plot_colour_mag_diagram(params,mags, colours, local_mags, local_colours, 
                            source, blend, RC, blue_filter, red_filter, 
                            yaxis_filter, tol):
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
                 marker='d',markersize=6, label='Blend')
    
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
    
    plt.legend()
    
    if red_filter == 'i' and blue_filter == 'r' and yaxis_filter == 'i':
        plt.axis([0.0,2.0,18.5,13.5])
    
    if red_filter == 'i' and blue_filter == 'r' and yaxis_filter == 'r':
        plt.axis([0.0,2.0,19.5,13.5])
        
    if red_filter == 'r' and blue_filter == 'g':
        plt.axis([0.0,3.0,21.0,13.5])
        
    plt.savefig(plot_file)

    plt.close(1)
    
    print 'Colour-magnitude diagram output to '+plot_file
    
def plot_colour_colour_diagram(params,star_catalog,catalog_header,target,RC,
                                                 det_idx,cat_idx,close_cat_idx,
                                                 use_isochrones=True):
    """Function to plot a colour-colour diagram, if sufficient data are
    available within the given star catalog"""
    
    tol = 2.0
    
    filters = { 'ip': 'SDSS-i', 'rp': 'SDSS-r', 'gp': 'SDSS-g' }
    
    try:
    
        inst_i = star_catalog['cal_ref_mag_ip'][det_idx]
        inst_r = star_catalog['cal_ref_mag_rp'][det_idx]
        inst_g = star_catalog['cal_ref_mag_gp'][det_idx]
        inst_gr = inst_g - inst_r - RC['Egr']
        inst_ri = inst_r - inst_i - RC['Eri']
        
        linst_i = star_catalog['cal_ref_mag_ip'][close_cat_idx]
        linst_r = star_catalog['cal_ref_mag_rp'][close_cat_idx]
        linst_g = star_catalog['cal_ref_mag_gp'][close_cat_idx]
        lcal_i = star_catalog['imag'][close_cat_idx]
        lcal_r = star_catalog['rmag'][close_cat_idx]
        lcal_g = star_catalog['gmag'][close_cat_idx]
        linst_gr = linst_g - linst_r - RC['Egr']
        linst_ri = linst_r - linst_i - RC['Eri']
        
        fig = plt.figure(1,(10,10))
        
        ax = plt.axes()
        
        ax.scatter(inst_gr, inst_ri, 
                     c='#E1AE13', marker='.', s=1, label=None)
        
        ax.scatter(linst_gr, linst_ri, marker='*', s=4, c='#8c6931',
                 label='Stars < '+str(round(tol,1))+'arcmin of target')
                     
        if len(target) > 0 and 'gr_0' in target.keys() and 'ri_0' in target.keys():
            plt.plot(target['gr_0'], target['ri_0'],'md',markersize=6)
        
        
        use_isochrones=False
        if use_isochrones and 'none' not in str(params['isochrone_file']).lower():
            
            isochrone_data = isochrone_utilities.read_PARSEC_table(params['isochrone_file'])

            fig = isochrone_utilities.overlay_isochrones(fig,isochrone_data,n_steps=6,label_plot=False)
        
        (spectral_type, luminosity_class, gr_colour, ri_colour) = spectral_type_data.get_spectral_class_data()
        
        plot_dwarfs = True
        plot_giants = True
        for i in range(0,len(spectral_type),1):
            
            spt = spectral_type[i]+luminosity_class[i]
            
            if luminosity_class[i] == 'V':
                c = 'g'
            else:
                c = 'k'
                        
            if luminosity_class[i] == 'III' and plot_giants:
                
                plt.plot(gr_colour[i], ri_colour[i], marker='s', color=c, alpha=0.5)

                plt.annotate(spt, (gr_colour[i], 
                               ri_colour[i]-0.1), 
                                 color=c, size=10, 
                                 rotation=-30.0, alpha=1.0)

            if luminosity_class[i] == 'V' and plot_dwarfs:
                
                plt.plot(gr_colour[i], ri_colour[i], marker='v', color=c, alpha=0.5)

                plt.annotate(spt, (gr_colour[i], 
                               ri_colour[i]+0.1), 
                                 color=c, size=10, 
                                 rotation=-30.0, alpha=1.0)

        plt.xlabel('SDSS (g-r) [mag]')
    
        plt.ylabel('SDSS (r-i) [mag]')
        
        plot_file = path.join(params['red_dir'],'colour_colour_diagram.png')
        
        plt.axis([-1.0,2.0,-0.5,1.0])
    
        plt.grid()
        
        plt.legend()
    
        plt.savefig(plot_file)
    
        plt.close(1)
        
        print 'Colour-colour diagram output to '+plot_file
        
    except AttributeError:
            
        print 'Warning: Insufficient data for colour-colour diagram'
        
def calibrate_instrumental_colour_colour_diagram(params,star_catalog,
                                                 catalog_header,target,
                                                 det_idx,cat_idx,close_cat_idx,
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
            
            print 'Calibration colour-colour diagram output to '+plot_file
            
            print 'Measured offsets in photometry:'
            print deltas
            
        except AttributeError:
            
            deltas = {}
            print 'Warning: Insufficient data for colour-colour diagram'
        
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

def localize_red_clump(star_catalog,close_cat_idx):
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
    inst_gr = inst_g - inst_r  
    cal_ri = cal_r - cal_i 
    cal_gr = cal_g - cal_r
    
    print 'Median (r-i), i: ',np.median(inst_ri), np.median(inst_i)
    print 'Median (g-r), g: ',np.median(inst_gr), np.median(inst_g)
    
    ri_min = 0.8 
    ri_max = 1.1 
    i_min = 15.5
    i_max = 16.5
    
    r_min = 16.2
    r_max = 17.5
    
    gr_min = 1.5 
    gr_max = 2.2 
    g_min = 18.5
    g_max = 19.5
    
    print('Selected Red Clump giants between:')
    print('i = '+str(i_min)+' to '+str(i_max))
    print('r = '+str(r_min)+' to '+str(r_max))
    print('(r-i) = '+str(ri_min)+' to '+str(ri_max))
    print('g = '+str(g_min)+' to '+str(g_max))
    print('(g-r) = '+str(gr_min)+' to '+str(gr_max))
    
    idx = select_within_range(inst_i, inst_ri, i_min, i_max, ri_min, ri_max)
    
    RC.ri = np.median(inst_ri[idx])
    RC.sig_ri = np.sqrt( ((inst_ri[idx] - RC.ri)**2).sum() / float(len(idx)) )
    RC.i = np.median(inst_i[idx])
    RC.sig_i = np.sqrt( ((inst_i[idx] - RC.i)**2).sum() / float(len(idx)) )
    
    idx = select_within_range(inst_r, inst_ri, r_min, r_max, ri_min, ri_max)
    
    RC.r = np.median(inst_r[idx])
    RC.sig_r = np.sqrt( ((inst_r[idx] - RC.r)**2).sum() / float(len(idx)) )
    
    idx = select_within_range(inst_g, inst_gr, g_min, g_max, gr_min, gr_max)
    
    RC.gr = np.median(inst_gr[idx])
    RC.sig_gr = np.sqrt( ((inst_gr[idx] - RC.gr)**2).sum() / float(len(idx)) )
    RC.g = np.median(inst_g[idx])
    RC.sig_g = np.sqrt( ((inst_g[idx] - RC.g)**2).sum() / float(len(idx)) )
    
    print('\nCentroid of Red Clump Stars at:')
    print(RC.summary(show_mags=True))
    print(RC.summary(show_colours=True))
    
    RC.transform_to_JohnsonCousins()
    
    print(RC.summary(johnsons=True))
    
    return RC

def measure_RC_offset(params,RC,target):
    """Function to calculate the offset of the Red Clump from its expected 
    values, taken from Bensby et al. (2017), 2017, A&A, 605A, 89 for V, I bands and
    Hawkins et al. (2017) MNRAS, 471, 722 for 2MASS J,H,Ks.
    """
    
    in_use = False
    
    RC = red_clump_utilities.get_essential_parameters(RC=RC)
    
    if in_use:
        RC.transform_2MASS_to_SDSS()
    
    print('\n Red Clump colours and absolute SDSS magnitudes:')
    print('Mg_RC,0 = '+str(RC['M_g_0'])+' +/- '+str(RC['sigMg_0'])+'mag')
    print('Mr_RC,0 = '+str(RC['M_r_0'])+' +/- '+str(RC['sigMr_0'])+'mag')
    print('Mi_RC,0 = '+str(RC['M_i_0'])+' +/- '+str(RC['sigMi_0'])+'mag')
    print('(g-r)_RC,0 = '+str(RC['g-r_0'])+' +/- '+str(RC['siggr_0'])+'mag')
    print('(r-i)_RC,0 = '+str(RC['r-i_0'])+' +/- '+str(RC['sigri_0'])+'mag')
    
    if in_use:
        print('\n Red Clump NIR colours and magnitudes:')
        print('J_RC,0 = '+str(RC['M_J_0'])+' +/- '+str(RC['sig_J_0'])+'mag')
        print('H_RC,0 = '+str(RC['M_H_0'])+' +/- '+str(RC['sig_H_0'])+'mag')
        print('Ks_RC,0 = '+str(RC['M_Ks_0'])+' +/- '+str(RC['sig_Ks_0'])+'mag')
        print('(J-H)_RC,0 = '+str(RC['J-H_0'])+' +/- '+str(RC['sigJH_0'])+'mag')
        print('(H-Ks)_RC,0 = '+str(RC['H-K_0'])+' +/- '+str(RC['sigHK_0'])+'mag')
    
    RC.D = red_clump_utilities.calc_red_clump_distance(params['target_ra'],params['target_dec'])
    RC = red_clump_utilities.calc_apparent_magnitudes(RC)
    
    print('\n Red Clump apparent SDSS magnitudes at distance '+str(RC['D_RC'])+'Kpc')
    print('g_RC,app = '+str(RC['m_g_0'])+' +/- '+str(RC['sigmg_0'])+'mag')
    print('r_RC,app = '+str(RC['m_r_0'])+' +/- '+str(RC['sigmr_0'])+'mag')
    print('i_RC,app = '+str(RC['m_i_0'])+' +/- '+str(RC['sigmi_0'])+'mag')
    
    if in_use:
        RC.transform_to_JohnsonCousins()
        
        print('\n Derived Red Clump instrumental colours and magnitudes:')
        print(RC.summary(johnsons=True))
    
    if in_use:
        RC.A_I = RC.I - RC.I_app
        RC.sig_A_I = RC.sig_I
        RC.EVI = RC.VI - RC.VI_0
        RC.sig_EVI = RC.sig_VI
    
    
    RC.A_g = RC.g - RC.m_g_0
    RC.sig_A_g = np.sqrt(RC.sig_g*RC.sig_g + RC.sig_mg_0*RC.sig_mg_0)
    RC.A_r = RC.r - RC.m_r_0
    RC.sig_A_r = np.sqrt(RC.sig_r*RC.sig_r + RC.sig_mr_0*RC.sig_mr_0)
    RC.A_i = RC.i - RC.m_i_0
    RC.sig_A_i = np.sqrt(RC.sig_i*RC.sig_i + RC.sig_mi_0*RC.sig_mi_0)
    
    RC.Egr = RC.gr - RC.gr_0
    RC.sig_Egr = np.sqrt( (RC.sig_gr*RC.sig_gr) + (RC.sig_gr_0*RC.sig_gr_0) )
    RC.Eri = RC.ri - RC.ri_0
    RC.sig_Eri = np.sqrt( (RC.sig_ri*RC.sig_ri) + (RC.sig_ri_0*RC.sig_ri_0) )

    print('\nExtinction, d(g) = '+str(RC.A_g)+' +/- '+str(RC.sig_A_g)+'mag')
    print('Extinction, d(r) = '+str(RC.A_r)+' +/- '+str(RC.sig_A_r)+'mag')
    print('Extinction, d(i) = '+str(RC.A_i)+' +/- '+str(RC.sig_A_i)+'mag')
    print('Reddening, E(g-r) = '+str(RC.Egr)+' +/- '+str(RC.sig_Egr)+'mag')
    print('Reddening, E(r-i) = '+str(RC.Eri)+' +/- '+str(RC.sig_Eri)+'mag')
    
    return RC
    
def calc_phot_properties(target, source, blend, RC):
    """Function to calculate the de-reddened and extinction-corrected 
    photometric properties of the target
    """
    in_use = False
    
    target.calibrate_phot_properties(RC)
    source.calibrate_phot_properties(RC)
    blend.calibrate_phot_properties(RC)

    print('\nSource star extinction-corrected magnitudes and de-reddened colours:')
    print(source.summary(show_cal=True))
    print(source.summary(show_cal=True,show_colours=True))
    
    print('\nBlend extinction-corrected magnitudes and de-reddened colours:')
    print(blend.summary(show_cal=True))
    print(blend.summary(show_cal=True,show_colours=True))
    
    return target,source,blend
    
if __name__ == '__main__':
    
    perform_colour_analysis()
    