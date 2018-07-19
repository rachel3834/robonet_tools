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

def plot_colour_diagrams():
    """Function to plot colour magnitude and colour-colour plots"""
    
    tol = 2.0       # Arcmin
    
    params = get_args()
    
    (star_catalog,catalog_header) = read_combined_star_catalog(params)
    
    target = find_target_data(params,star_catalog)

    print target
    
    (det_idx, cat_idx, close_cat_idx) = index_valid_star_entries(star_catalog,
                                                                target,tol,
                                                                valid_cat=True)
                                                            
    deltas = calibrate_instrumental_colour_colour_diagram(params,star_catalog,
                                                 catalog_header,target,
                                                 det_idx,cat_idx,close_cat_idx)
    
    (RC_i, sig_RC_i, RC_ri, sig_RC_ri) = localize_red_clump(star_catalog,close_cat_idx,deltas)
    
    plot_colour_mag_diagram(params,star_catalog,deltas,catalog_header,target,
                                                 det_idx,cat_idx,close_cat_idx,
                                                 RC_i, sig_RC_i, RC_ri, sig_RC_ri)
    
    plot_colour_colour_diagram(params,star_catalog,deltas,catalog_header,target,
                                                 det_idx,cat_idx,close_cat_idx)
    
def get_args():
    """Function to gather the necessary commandline arguments"""

    params = {}
    
    if len(argv) == 1:
        
        params['catalog_file'] = raw_input('Please enter the path to combined star catalog file: ')
        params['red_dir'] = raw_input('Please enter the path to the output directory: ')
        params['target_ra'] = raw_input('Please enter the RA of target to highlight, if any [sexigesimal format or none]:')
        params['target_dec'] = raw_input('Please enter the Dec of target to highlight, if any [sexigesimal format or none]:')
        params['i_offset'] = float(raw_input('Please enter the extinction correction in i [mag]: '))
        params['gr_offset'] = float(raw_input('Please enter the reddening correction in g-r [mag]: '))
        params['ri_offset'] = float(raw_input('Please enter the reddening correction in r-i [mag]: '))
        
    else:

        params['catalog_file'] = argv[1]
        params['red_dir'] = argv[2]
        params['target_ra'] = argv[3]
        params['target_dec'] = argv[4]
        params['i_offset'] = float(argv[5])
        params['gr_offset'] = float(argv[6])
        params['ri_offset'] = float(argv[7])
    
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
            target['ref_mag_ip'] = star_catalog['ref_mag_ip'][match_data[1][idx[0]]]
            target['ref_mag_err_ip'] = star_catalog['ref_mag_err_ip'][match_data[1][idx[0]]]
            target['ref_mag_rp'] = star_catalog['ref_mag_rp'][match_data[1][idx[0]]]
            target['ref_mag_err_rp'] = star_catalog['ref_mag_err_rp'][match_data[1][idx[0]]]
            target['separation'] = match_data[2][idx[0]].to_string(unit=u.arcsec)
            try:
                target['ref_mag_gp'] = star_catalog['ref_mag_gp'][match_data[1][idx[0]]]
            except AttributeError:
                pass
            print 'Target identified as '+repr(target)
    
    return target

def index_valid_star_entries(star_catalog,target,tol,valid_cat=False):
    """Function to return an index of all stars with both full instrumental and
    catalogue entries"""
    
    idx1 = np.where(star_catalog['ref_mag_ip'] > 0.0)[0]
    idx2 = np.where(star_catalog['ref_mag_rp'] > 0.0)[0]
    idx3 = np.where(star_catalog['ref_mag_gp'] > 0.0)[0]
    
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
    
def plot_colour_mag_diagram(params,star_catalog,deltas,catalog_header,target,
                                                 det_idx,cat_idx,close_cat_idx,
                                                 RC_i, sig_RC_i, RC_ri, sig_RC_ri):
    """Function to plot a colour-magnitude diagram"""
    
    tol = 2.0
    
    filters = { 'ip': 'SDSS-i', 'rp': 'SDSS-r', 'gp': 'SDSS-g' }
    
    inst_i = star_catalog['ref_mag_ip'][det_idx]
    inst_r = star_catalog['ref_mag_rp'][det_idx]
    cal_i = star_catalog['imag'][cat_idx]
    cal_r = star_catalog['rmag'][cat_idx]
    inst_ri = inst_r - inst_i    # Catalogue column order is red -> blue
    cal_ri = cal_r - cal_i
    
    if len(target) > 0 and target['ref_mag_ip'] > 0.0 and \
            target['ref_mag_rp'] > 0.0:
        target_ri = target['ref_mag_rp'] - target['ref_mag_ip']
        sig_target_ri = np.sqrt( (target['ref_mag_err_rp']*target['ref_mag_err_rp']) + \
                        (target['ref_mag_err_ip']*target['ref_mag_err_ip']) ) 
        
        print('Target i_inst = '+str(target['ref_mag_ip'])+' +/- '+str(target['ref_mag_err_ip']))
        print('Target (r-i)_inst = '+str(target_ri)+' +/- '+str(sig_target_ri))
        
    fig = plt.figure(1,(10,10))
    
    plt.rcParams.update({'font.size': 18})
        
    plt.scatter(inst_ri+deltas['dri']+params['ri_offset'],
                inst_i+deltas['di'],
                 c='#E1AE13', 
                 marker='.', s=1, label=None)
    
#    plt.scatter(cal_ri, cal_r, c='#8c6931', 
#             marker='.', s=1, 
#             label='Catalogue')
    
    linst_i = star_catalog['ref_mag_ip'][close_cat_idx]
    linst_r = star_catalog['ref_mag_rp'][close_cat_idx]
    lcal_i = star_catalog['imag'][close_cat_idx]
    lcal_r = star_catalog['rmag'][close_cat_idx]
    linst_ri = linst_r - linst_i    # Catalogue column order is red -> blue
    lcal_ri = lcal_r - lcal_i
    
    plt.scatter(linst_ri+deltas['dri']+params['ri_offset'],
                linst_i+deltas['di'],
                 c='#8c6931', 
                 marker='*', s=4, 
                 label='Stars < '+str(round(tol,1))+'arcmin of target')
    
    if len(target) > 0 and target['ref_mag_ip'] > 0.0 and \
            target['ref_mag_rp'] > 0.0:
                
        plt.errorbar(target_ri+deltas['dri']+params['ri_offset'], 
                 target['ref_mag_ip']+deltas['di'], 
                 yerr=target['ref_mag_err_ip'],
                 xerr = sig_target_ri,color='m',
                 marker='d',markersize=6, label='Source at baseline')
    
    plt.errorbar(RC_ri+deltas['dri']+params['ri_offset'], 
                 RC_i+deltas['di'], 
                 yerr=sig_RC_i,
                 xerr = sig_RC_ri,color='g',
                 marker='s',markersize=6, label='Red Clump centroid')
                 
    plt.xlabel('SDSS (r-i) [mag]')

    plt.ylabel('SDSS-i [mag]')
    
    [xmin,xmax,ymin,ymax] = plt.axis()
    
    plt.axis([xmin,xmax,ymax,ymin])
    
    plot_file = path.join(params['red_dir'],'colour_magnitude_diagram.png')

    plt.grid()
    
    plt.legend()
    
    plt.axis([0.0,2.0,18.0,13.5])
    
    plt.savefig(plot_file)

    plt.close(1)
    
    print 'Colour-magnitude diagram output to '+plot_file
    
def plot_colour_colour_diagram(params,star_catalog,deltas,catalog_header,target,
                                                 det_idx,cat_idx,close_cat_idx):
    """Function to plot a colour-colour diagram, if sufficient data are
    available within the given star catalog"""
    
    tol = 2.0
    
    filters = { 'ip': 'SDSS-i', 'rp': 'SDSS-r', 'gp': 'SDSS-g' }
    
    try:
    
        inst_i = star_catalog['ref_mag_ip'][det_idx]
        inst_r = star_catalog['ref_mag_rp'][det_idx]
        inst_g = star_catalog['ref_mag_gp'][det_idx]
        inst_gr = inst_g - inst_r    # Catalogue column order is red -> blue
        inst_ri = inst_r - inst_i   
        
        if len(target) > 0 and target['ref_mag_ip'] > 0.0 and \
            target['ref_mag_rp'] > 0.0 and target['ref_mag_gp']:
            target_gr = target['ref_mag_gp'] - target['ref_mag_rp']
            target_ri = target['ref_mag_rp'] - target['ref_mag_ip']
    
        fig = plt.figure(1,(10,10))
    
        plt.scatter(inst_gr+params['gr_offset'], 
                    inst_ri+params['ri_offset'], 
                     c='#2b8c85', 
                     marker='.', s=1)
        
        if len(target) > 0 and target['ref_mag_ip'] > 0.0 and \
            target['ref_mag_rp'] > 0.0 and target['ref_mag_gp']:
            plt.plot(target_gr+params['gr_offset'], 
                     target_ri+params['ri_offset'],
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
                                                 det_idx,cat_idx,close_cat_idx):
    """Function to plot a colour-colour diagram, if sufficient data are
    available within the given star catalog"""
    
    tol = 2.0
    
    filters = { 'ip': 'SDSS-i', 'rp': 'SDSS-r', 'gp': 'SDSS-g' }
    
    try:
    
        inst_i = star_catalog['ref_mag_ip'][close_cat_idx]
        inst_r = star_catalog['ref_mag_rp'][close_cat_idx]
        inst_g = star_catalog['ref_mag_gp'][close_cat_idx]
        inst_gr = inst_g - inst_r    # Catalogue column order is red -> blue
        inst_ri = inst_r - inst_i   
        
        cal_i = star_catalog['imag'][close_cat_idx]
        cal_r = star_catalog['rmag'][close_cat_idx]
        cal_g = star_catalog['gmag'][close_cat_idx]
        cal_gr = cal_g - cal_r    # Catalogue column order is red -> blue
        cal_ri = cal_r - cal_i    
        
        deltas = {}
        deltas['di'] = np.median(cal_i - inst_i)
        deltas['dr'] = np.median(cal_r - inst_r)
        deltas['dg'] = np.median(cal_g - inst_g)
        deltas['dgr'] = np.median(cal_gr - inst_gr)
        deltas['dri'] = np.median(cal_ri - inst_ri)
        
        if len(target) > 0 and target['ref_mag_ip'] > 0.0 and \
            target['ref_mag_rp'] > 0.0 and target['ref_mag_gp']:
            target_gr = target['ref_mag_gp'] - target['ref_mag_rp']
            target_ri = target['ref_mag_rp'] - target['ref_mag_ip']
    
        fig = plt.figure(1, (10,10))
        
        plt.rcParams.update({'font.size': 18})
        
        plt.plot(inst_gr, inst_ri,'.',
                 color='#2b8c85',markersize=2,alpha=0.4,
                 label='Instrumental')
        plt.plot(cal_gr, cal_ri,'.',
                 color='#8c6931',markersize=2,alpha=0.4,
                 label='Catalogue')
        
        plt.plot(inst_gr+deltas['dgr'], inst_ri+deltas['dri'],'.',
                 color='m',markersize=2,alpha=0.4,
                 label='Calibrated')
                 
#        if len(target) > 0 and target['mag1'] > 0.0 and \
#            target['mag2'] > 0.0 and target['mag2']:
#            plt.plot(target_colour2, target_colour1,'md',markersize=6)
            
        plt.xlabel('SDSS (g-r) [mag]')
    
        plt.ylabel('SDSS (r-i) [mag]')
        
        plot_file = path.join(params['red_dir'],'calib_colour_colour_diagram.png')
        
        plt.grid()
        
        plt.axis([0.0,3.0,-1.0,2.0])
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

def localize_red_clump(star_catalog,close_cat_idx,deltas):
    """Function to calculate the centroid of the Red Clump stars in a 
    colour-magnitude diagram"""
    
    inst_i = star_catalog['ref_mag_ip'][close_cat_idx]
    inst_r = star_catalog['ref_mag_rp'][close_cat_idx]
    cal_i = star_catalog['imag'][close_cat_idx]
    cal_r = star_catalog['rmag'][close_cat_idx]
    inst_ri = inst_r - inst_i    # Catalogue column order is red -> blue
    cal_ri = cal_r - cal_i
    
    print 'Median: ',np.median(inst_ri), np.median(inst_i)
    
    ri_min = 0.8 - deltas['dri']
    ri_max = 1.1 - deltas['dri']
    i_min = 15.5 - deltas['di']
    i_max = 16.5 - deltas['di']
    
    print('Selected Red Clump giants between:')
    print('i = '+str(i_min)+' to '+str(i_max))
    print('(r-i) = '+str(ri_min)+' to '+str(ri_max))
    
    idx1 = np.where(inst_ri >= ri_min)[0]
    idx2 = np.where(inst_ri <= ri_max)[0]
    idx3 = np.where(inst_i >= i_min)[0]
    idx4 = np.where(inst_i <= i_max)[0]
    idx = set(idx1).intersection(set(idx2))
    idx = idx.intersection(set(idx3))
    idx = list(idx.intersection(set(idx4)))
    
    RC_ri = np.median(inst_ri[idx])
    sig_RC_ri = np.sqrt( ((inst_ri[idx] - RC_ri)**2).sum() / float(len(idx)) )
    RC_i = np.median(inst_i[idx])
    sig_RC_i = np.sqrt( ((inst_i[idx] - RC_i)**2).sum() / float(len(idx)) )
    
    print('\nCentroid of Red Clump Stars at:')
    print('i_inst,RC= '+str(RC_i)+' +/- '+str(sig_RC_i))
    print('(r-i)_inst,RC= '+str(RC_ri)+' +/- '+str(sig_RC_ri))
    print('i_0,RC= '+str(RC_i+deltas['di'])+' +/- '+str(sig_RC_i))
    print('(r-i)_0,RC= '+str(RC_ri+deltas['dri'])+' +/- '+str(sig_RC_ri))
    
    return RC_i, sig_RC_i, RC_ri, sig_RC_ri
    
if __name__ == '__main__':
    
    plot_colour_diagrams()
    