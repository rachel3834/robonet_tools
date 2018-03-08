# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 11:03:12 2018

@author: rstreet
"""

import os
import sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import glob
import gather_header_qc_data

def analyze_qc_data():
    """Function to conduct analysis of QC data including manually-evaluated
    bad images"""

    dir_path = gather_header_qc_data.get_args()
    
    file_path = os.path.join(dir_path,'image_qc.dat')
    (image_list,image_data, tqc) = read_image_qc_stats(file_path)

    file_path = os.path.join(dir_path,'bad_images_manual_qc.dat')
    bad_images = read_image_qc_reasons(file_path)

    gather_header_qc_data.plot_stats_for_sites(image_list,image_data,dir_path,'sky_temp.png', 7, 
    'Temperature [$^{\circ}$ C]', 'Measured sky temperature across all sites',
    bad_images=bad_images)
    
    gather_header_qc_data.plot_delta_stats_for_sites(image_list,image_data,dir_path,'delta_sky_mag.png', 9, 8, 
    '$\Delta$ sky brightness [mag]', 'Difference between measured and computed sky brightness',
    bad_images=bad_images)
    
    gather_header_qc_data.plot_stats_for_sites(image_list,image_data,dir_path,'focal_position.png', 13, 
    'Focus offset [mm]', 'Actual focal position for all sites',
    bad_images=bad_images)
    gather_header_qc_data.plot_stats_for_sites(image_list,image_data,dir_path,'focal_zp.png', 14, 
    'Focus offset [mm]', 'Telescope focal zeropoint',
    bad_images=bad_images)
    gather_header_qc_data.plot_stats_for_sites(image_list,image_data,dir_path,'focal_therm_offset.png', 15, 
    'Focus offset [mm]', 'Focal thermal offset',
    bad_images=bad_images)
    
    gather_header_qc_data.plot_stats_for_sites(image_list,image_data,dir_path,'moon_frac.png', 11, 
    'Moon fraction', 'Moon fraction for observations',
    bad_images=bad_images)
    
    gather_header_qc_data.plot_compare_stats_for_sites(image_list,image_data,dir_path,'moon_distance.png', 12, 8, 
    'Moon separation [$^{\circ}$]','Sky brightness [mag]', 
    'Measured sky brightness as a function of lunar separation',
    bad_images=bad_images)
    
    gather_header_qc_data.plot_stats_per_filter(image_list,image_data,dir_path,'mean_sky.png', 17, 
    'Mean sky value [ADU]', 'Average sky background data',
    bad_images=bad_images)
    
    gather_header_qc_data.plot_stats_per_filter(image_list,image_data,dir_path,'mean_sky_sigma.png', 18, 
    'Mean sky sigma [ADU]', 'Average sky background deviation',
    bad_images=bad_images)
    
    
    gather_header_qc_data.plot_stats_for_sites(image_list,image_data,dir_path,'mean_fwhm.png', 16, 
    'Mean FWHM [arcsec]', 'Average FWHM',exclude_no_data=True,
    bad_images=bad_images)
    
    gather_header_qc_data.plot_stats_for_sites(image_list,image_data,dir_path,'mean_ellipticity.png', 19, 
    'Mean ellipticity', 'Average ellipticity',exclude_no_data=True,
    bad_images=bad_images)
    
    plot_qc_pie_charts(dir_path,bad_images)
    
def plot_qc_pie_charts(dir_path,bad_images):
    """Function to plot a pie chart of the reasons for image rejections"""

    fig = plt.figure(2)

    reject_count = {}
    
    for image,tlist in bad_images.items():
        
        for reason in tlist:
            
            if reason not in reject_count.keys():
                
                reject_count[reason] = 0
            
            reject_count[reason] += 1

    reject_reasons = []
    counts = []
    
    for reason, n in reject_count.items():
        
        reject_reasons.append(reason)
        
        counts.append( (float(n)/float(len(bad_images)))*100.0 )

    plt.pie(counts, labels=reject_reasons, autopct='%1.1f%%', startangle=90)
    plt.axis('equal')  # Equal aspect r

    plt.savefig( os.path.join(dir_path,'image_reject_reasons.png') )

    plt.close(2)
    
def read_image_qc_reasons(file_path):
    """Function to read a list of rejected images and the reasons why"""
    
    qc_reasons = {}
    
    if os.path.isfile(file_path):
        
        file_lines = open(file_path,'r').readlines()
        
        for line in file_lines:
            
            if line[0:1] != '#':
                
                entries = line.replace('\n','').split()
                
                qc_reasons[entries[0]] = entries[1:]

    return qc_reasons
    
def read_image_qc_stats(file_path):
    """Function to read in an existing file of image quality control statistics"""
    
    n_metrics = 20
        
    if os.path.isfile(file_path):
        
        file_lines = open(file_path,'r').readlines()
        
        image_list = []
        image_data = []
        qc = []
        
        for line in file_lines:
            
            if line[0:1] != '#':
                
                entries = line.replace('\n','').split()
                
                image_list.append(entries[0:3])
                
                tlist = []
                
                for e in entries[3:]:
                    
                    tlist.append(float(e))
                    
                image_data.append(tlist)
                
                if len(entries) > n_metrics+2:
                    
                    qc.append(entries[n_metrics+2:])
                
        image_list = np.array(image_list)
        
        image_data = np.array(image_data)
        
        qc = np.array(qc)
        
    return image_list, image_data, qc
    
if __name__ == '__main__':
    
    analyze_qc_data()
    