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
import red_clump_utilities

def perform_colour_analysis():
    """Function to plot colour magnitude and colour-colour plots"""
    
    tol = 2.0       # Arcmin
    calib_on_colours = False
    
    params = get_args()
    
    (star_catalog,catalog_header) = read_combined_star_catalog(params)
    
    target = find_target_data(params,star_catalog)

    print target
    
    (det_idx, cat_idx, close_cat_idx) = index_valid_star_entries(star_catalog,
                                                                target,tol,
                                                                valid_cat=True)
    
    deltas = calibrate_instrumental_colour_colour_diagram(params,star_catalog,
                                                 catalog_header,target,
                                                 det_idx,cat_idx,close_cat_idx,
                                                 calib=calib_on_colours)

    RC = localize_red_clump(star_catalog,close_cat_idx)
    
    analyse_colour_mag_diagrams(params,star_catalog,catalog_header,target,
                                                 det_idx,cat_idx,close_cat_idx,
                                                 RC)
    
    plot_colour_colour_diagram(params,star_catalog,catalog_header,target,
                                                 det_idx,cat_idx,close_cat_idx)
                                                 
    measure_RC_offset(params,RC,target)
    
def get_args():
    """Function to gather the necessary commandline arguments"""

    params = {}
    
    if len(argv) == 1:
        
        params['catalog_file'] = raw_input('Please enter the path to combined star catalog file: ')
        params['red_dir'] = raw_input('Please enter the path to the output directory: ')
        params['target_ra'] = raw_input('Please enter the RA of target [sexigesimal]:')
        params['target_dec'] = raw_input('Please enter the Dec of target [sexigesimal]:')
        
    else:

        params['catalog_file'] = argv[1]
        params['red_dir'] = argv[2]
        params['target_ra'] = argv[3]
        params['target_dec'] = argv[4]
    
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

def find_target_data(params,star_catalog):
    """Function to identify the photometry for a given target, if present
    in the star catalogue"""
    
    target = {}
    
    if params['target_ra'] != None:
        
        target_location = SkyCoord([params['target_ra']], [params['target_dec']], unit=(u.hourangle, u.deg))
                
        stars = SkyCoord(star_catalog['RA'], star_catalog['DEC'], unit="deg")
        
        tolerance = 2.0 * u.arcsec
        
        match_data = matching.search_around_sky(target_location, stars, 
                                                seplimit=tolerance)    
                                                
        idx = np.argsort(match_data[2].value)
    
        if len(match_data[0]) > 0:
            target['star_index'] = star_catalog['star_index'][match_data[1][idx[0]]]
            target['ra'] = star_catalog['RA'][match_data[1][idx[0]]]
            target['dec'] = star_catalog['DEC'][match_data[1][idx[0]]]
            target['cal_ref_mag_ip'] = star_catalog['cal_ref_mag_ip'][match_data[1][idx[0]]]
            target['cal_ref_mag_err_ip'] = star_catalog['cal_ref_mag_err_ip'][match_data[1][idx[0]]]
            target['cal_ref_mag_rp'] = star_catalog['cal_ref_mag_rp'][match_data[1][idx[0]]]
            target['cal_ref_mag_err_rp'] = star_catalog['cal_ref_mag_err_rp'][match_data[1][idx[0]]]
            target['ref_mag_ip'] = star_catalog['ref_mag_ip'][match_data[1][idx[0]]]
            target['ref_mag_err_ip'] = star_catalog['ref_mag_err_ip'][match_data[1][idx[0]]]
            target['ref_mag_rp'] = star_catalog['ref_mag_rp'][match_data[1][idx[0]]]
            target['ref_mag_err_rp'] = star_catalog['ref_mag_err_rp'][match_data[1][idx[0]]]
            target['separation'] = match_data[2][idx[0]].to_string(unit=u.arcsec)
            try:
                target['cal_ref_mag_gp'] = star_catalog['cal_ref_mag_gp'][match_data[1][idx[0]]]
                target['cal_ref_mag_err_gp'] = star_catalog['cal_ref_mag_err_gp'][match_data[1][idx[0]]]
                target['ref_mag_gp'] = star_catalog['ref_mag_gp'][match_data[1][idx[0]]]
                target['ref_mag_err_gp'] = star_catalog['ref_mag_err_gp'][match_data[1][idx[0]]]
            except AttributeError:
                pass
            print 'Target identified as '+repr(target)
    
        if len(target) > 0 and target['cal_ref_mag_ip'] > 0.0 and \
                                target['cal_ref_mag_rp'] > 0.0:
                                    
            target['ri'] = target['cal_ref_mag_rp'] - target['cal_ref_mag_ip']
            target['sig_ri'] = np.sqrt( (target['cal_ref_mag_err_rp']*target['cal_ref_mag_err_rp']) + \
                            (target['cal_ref_mag_err_ip']*target['cal_ref_mag_err_ip']) ) 
            
            print('Target i_inst = '+str(target['cal_ref_mag_ip'])+' +/- '+str(target['cal_ref_mag_err_ip']))
            print('Target (r-i)_inst = '+str(target['ri'])+' +/- '+str(target['sig_ri']))
            
            target_phot = jester_phot_transforms.transform_SDSS_to_JohnsonCousins(ri=target['ri'], 
                                                                                  sigri=target['sig_ri'])
            target['V-R'] = target_phot['V-R']
            target['sigVR'] = target_phot['sigVR']
            target['Rc-Ic'] = target_phot['Rc-Ic']
            target['sigRI'] = target_phot['sigRI']
            
            print('\nTarget (V-R)_inst = '+str(target_phot['V-R'])+' +/- '+str(target_phot['sigVR'])+'mag')
            print('Target (Rc-Ic)_inst = '+str(target_phot['Rc-Ic'])+' +/- '+str(target_phot['sigRI'])+'mag')
    
        if len(target) > 0 and target['cal_ref_mag_rp'] > 0.0 and \
                                target['cal_ref_mag_gp'] > 0.0:
                    
            target['gr'] = target['cal_ref_mag_gp'] - target['cal_ref_mag_rp']
            target['sig_gr'] = np.sqrt( (target['cal_ref_mag_err_gp']*target['cal_ref_mag_err_gp']) + \
                                        (target['cal_ref_mag_err_rp']*target['cal_ref_mag_err_rp']) ) 
            
            print('Target g_inst = '+str(target['cal_ref_mag_gp'])+' +/- '+str(target['cal_ref_mag_err_gp']))
            print('Target (g-r)_inst = '+str(target['gr'])+' +/- '+str(target['sig_gr']))
            
            target_phot = jester_phot_transforms.transform_SDSS_to_JohnsonCousins(g=target['cal_ref_mag_gp'],
                                                                                  sigg=target['cal_ref_mag_err_gp'],
                                                                                  gr=target['gr'], 
                                                                                  siggr=target['sig_gr'])
            
            target['V'] = target_phot['V']
            target['sigV'] = target_phot['sigV']
            target['B-V'] = target_phot['B-V']
            target['sigBV'] = target_phot['sigBV']
            
            print('\nTarget V_inst = '+str(target_phot['V'])+' +/- '+str(target_phot['sigV'])+'mag')
            print('Target (B-V)_inst = '+str(target_phot['B-V'])+' +/- '+str(target_phot['sigBV'])+'mag')
    
        if 'V' in target.keys() and 'V-R' in target.keys():
            
            target = jester_phot_transforms.calc_derived_colours_JohnsonCousins(target)
    
            print('\n Derived target instrumental colours and magnitudes:')
            print('R_inst = '+str(target['R'])+' +/- '+str(target['sigR'])+'mag')
            print('I_inst = '+str(target['I'])+' +/- '+str(target['sigI'])+'mag')
            print('(V-I)_inst = '+str(target['V-I'])+' +/- '+str(target['sigVI'])+'mag')
    

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
    
def analyse_colour_mag_diagrams(params,star_catalog,catalog_header,target,
                                                 det_idx,cat_idx,close_cat_idx,
                                                 RC):
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
                            target, RC, 'r', 'i', 'i', tol)
                            
    plot_colour_mag_diagram(params,inst_g, inst_gr, linst_g, linst_gr, 
                            target, RC, 'g', 'r', 'g', tol)
    
    
def plot_colour_mag_diagram(params,mags, colours, local_mags, local_colours, 
                            target, RC, blue_filter, red_filter, 
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
    
    blue_key = 'cal_ref_mag_'+blue_filter+'p'
    red_key = 'cal_ref_mag_'+red_filter+'p'
    col_key = blue_filter+red_filter
    col_err_key = 'sig_'+blue_filter+red_filter
    mag_key = 'cal_ref_mag_'+yaxis_filter+'p'
    mag_err_key = 'cal_ref_mag_err_'+yaxis_filter+'p'
    
    if len(target) > 0 and target[blue_key] > 0.0 and \
            target[red_key] > 0.0:
                
        plt.errorbar(target[col_key], target[mag_key], 
                 yerr = target[mag_err_key],
                 xerr = target[col_err_key], color='m',
                 marker='d',markersize=6, label='Source at baseline')
    
    col_key = blue_filter+red_filter
    col_err_key = 'sig_'+blue_filter+red_filter
    mag_err_key = 'sig_'+yaxis_filter
    plt.errorbar(RC[col_key], RC[yaxis_filter], 
                 yerr=RC[mag_err_key], xerr=RC[col_err_key],
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
    
    if red_filter == 'i' and blue_filter == 'r':
        plt.axis([0.0,2.0,18.0,13.5])
    
    if red_filter == 'r' and blue_filter == 'g':
        plt.axis([0.0,3.0,21.0,13.5])
        
    plt.savefig(plot_file)

    plt.close(1)
    
    print 'Colour-magnitude diagram output to '+plot_file
    
def plot_colour_colour_diagram(params,star_catalog,catalog_header,target,
                                                 det_idx,cat_idx,close_cat_idx):
    """Function to plot a colour-colour diagram, if sufficient data are
    available within the given star catalog"""
    
    tol = 2.0
    
    filters = { 'ip': 'SDSS-i', 'rp': 'SDSS-r', 'gp': 'SDSS-g' }
    
    try:
    
        inst_i = star_catalog['cal_ref_mag_ip'][det_idx]
        inst_r = star_catalog['cal_ref_mag_rp'][det_idx]
        inst_g = star_catalog['cal_ref_mag_gp'][det_idx]
        inst_gr = inst_g - inst_r    # Catalogue column order is red -> blue
        inst_ri = inst_r - inst_i   
        
        if len(target) > 0 and target['cal_ref_mag_ip'] > 0.0 and \
            target['cal_ref_mag_rp'] > 0.0 and target['cal_ref_mag_gp']:
            target_gr = target['cal_ref_mag_gp'] - target['cal_ref_mag_rp']
            target_ri = target['cal_ref_mag_rp'] - target['cal_ref_mag_ip']
    
        fig = plt.figure(1,(10,10))
    
        plt.scatter(inst_gr, inst_ri, 
                     c='#2b8c85', marker='.', s=1)
        
        if len(target) > 0 and target['cal_ref_mag_ip'] > 0.0 and \
            target['cal_ref_mag_rp'] > 0.0 and target['cal_ref_mag_gp']:
            plt.plot(target_gr, target_ri,
                     'md',markersize=6)
        
        (spectral_type, luminosity_class, gr_colour, ri_colour) = spectral_type_data.get_spectral_class_data()
        
        #(dwarf_stars, dwarf_labels) = star_colour_data.get_dwarf_star_colours()
        #(giant_stars, giant_labels) = star_colour_data.get_giant_star_colours()
        
        #plt.plot(dwarf_stars[:,1],dwarf_stars[:,2], 'b.')
        #plt.plot(giant_stars[:,1],giant_stars[:,2], 'r.')
        
        plot_dwarfs = False
        plot_giants = True
        for i in np.arange(len(spectral_type)):
            
            spt = spectral_type[i]+luminosity_class[i]
            
            if luminosity_class[i] == 'V':
                c = 'b'
            else:
                c = 'k'
                        
            if luminosity_class[i] == 'III' and plot_giants:
                
                plt.plot(gr_colour[i], ri_colour[i], c+'.', alpha=0.5)

                plt.annotate(spt, (gr_colour[i], 
                               ri_colour[i]-0.1), 
                                 color=c, size=10, 
                                 rotation=-30.0, alpha=0.5)

            if luminosity_class[i] == 'V' and plot_dwarfs:
                
                plt.plot(gr_colour[i], ri_colour[i], c+'.', alpha=0.5)

                plt.annotate(spt, (gr_colour[i], 
                               ri_colour[i]+0.1), 
                                 color=c, size=10, 
                                 rotation=-30.0, alpha=0.5)

        plt.xlabel('SDSS (g-r) [mag]')
    
        plt.ylabel('SDSS (r-i) [mag]')
        
        plot_file = path.join(params['red_dir'],'colour_colour_diagram.png')
        
        plt.axis([-1.0,2.0,-0.5,1.0])
    
        plt.grid()
        
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
    
    RC = {}
    
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
    
    gr_min = 1.5 
    gr_max = 2.2 
    g_min = 18.5
    g_max = 19.5
    
    print('Selected Red Clump giants between:')
    print('i = '+str(i_min)+' to '+str(i_max))
    print('(r-i) = '+str(ri_min)+' to '+str(ri_max))
    print('g = '+str(g_min)+' to '+str(g_max))
    print('(g-r) = '+str(gr_min)+' to '+str(gr_max))
    
    idx = select_within_range(inst_i, inst_ri, i_min, i_max, ri_min, ri_max)
    
    RC['ri'] = np.median(inst_ri[idx])
    RC['sig_ri'] = np.sqrt( ((inst_ri[idx] - RC['ri'])**2).sum() / float(len(idx)) )
    RC['i'] = np.median(inst_i[idx])
    RC['sig_i'] = np.sqrt( ((inst_i[idx] - RC['i'])**2).sum() / float(len(idx)) )
    
    idx = select_within_range(inst_g, inst_gr, g_min, g_max, gr_min, gr_max)
    
    RC['gr'] = np.median(inst_gr[idx])
    RC['sig_gr'] = np.sqrt( ((inst_gr[idx] - RC['gr'])**2).sum() / float(len(idx)) )
    RC['g'] = np.median(inst_g[idx])
    RC['sig_g'] = np.sqrt( ((inst_g[idx] - RC['g'])**2).sum() / float(len(idx)) )
    
    print('\nCentroid of Red Clump Stars at:')
    print('i_inst,RC= '+str(RC['i'])+' +/- '+str(RC['sig_i']))
    print('(r-i)_inst,RC= '+str(RC['ri'])+' +/- '+str(RC['sig_ri']))
    print('g_inst,RC= '+str(RC['g'])+' +/- '+str(RC['sig_g']))
    print('(g-r)_inst,RC= '+str(RC['gr'])+' +/- '+str(RC['sig_gr']))
    
    RC_phot = jester_phot_transforms.transform_SDSS_to_JohnsonCousins(ri=RC['ri'], 
                                                                      sigri=RC['sig_ri'])
    RC['V-R'] = RC_phot['V-R']
    RC['sigVR'] = RC_phot['sigVR']
    RC['Rc-Ic'] = RC_phot['Rc-Ic']
    RC['sigRI'] = RC_phot['sigRI']

    print('\n(V-R)_inst,RC = '+str(RC_phot['V-R'])+' +/- '+str(RC_phot['sigVR'])+'mag')
    print('(Rc-Ic)_inst,RC = '+str(RC_phot['Rc-Ic'])+' +/- '+str(RC_phot['sigRI'])+'mag')
    
    RC_phot = jester_phot_transforms.transform_SDSS_to_JohnsonCousins(g=RC['g'],
                                                                      sigg=RC['sig_g'],
                                                                      gr=RC['gr'], 
                                                                      siggr=RC['sig_gr'])
    RC['V'] = RC_phot['V']
    RC['sigV'] = RC_phot['sigV']
    RC['B-V'] = RC_phot['B-V']
    RC['sigBV'] = RC_phot['sigBV']

    print('\n(B-V)_inst,RC = '+str(RC_phot['B-V'])+' +/- '+str(RC_phot['sigBV'])+'mag')
    print('V_inst,RC = '+str(RC_phot['V'])+' +/- '+str(RC_phot['sigV'])+'mag')
    
    return RC

def measure_RC_offset(params,RC,target):
    """Function to calculate the offset of the Red Clump from its expected 
    values, taken from Nataf et al. (2013), ApJ, 769, 88"""
    
    RC['M_I_0'] = -0.12
    RC['V-I_0'] = 1.06
    
    RC['distance'] = red_clump_utilities.calc_red_clump_distance(params['target_ra'],params['target_dec'])
    RC['I_app'] = red_clump_utilities.calc_I_apparent(RC['distance'])
    
    RC = jester_phot_transforms.calc_derived_colours_JohnsonCousins(RC)
    
    print('\n Derived Red Clump instrumental colours and magnitudes:')
    print('R_inst,RC = '+str(RC['R'])+' +/- '+str(RC['sigR'])+'mag')
    print('I_inst,RC = '+str(RC['I'])+' +/- '+str(RC['sigI'])+'mag')
    print('(V-I)_inst,RC = '+str(RC['V-I'])+' +/- '+str(RC['sigVI'])+'mag')
    
    A = RC['I'] - RC['I_app']
    sigA = RC['sigI']
    EVI = RC['V-I'] - RC['V-I_0']
    sigEVI = RC['sigVI']
    
    print('\nExtinction, d(I) = '+str(A)+' +/- '+str(sigA))
    print('Reddening, E(V-I) = '+str(EVI)+' +/- '+str(sigEVI))
    
if __name__ == '__main__':
    
    perform_colour_analysis()
    