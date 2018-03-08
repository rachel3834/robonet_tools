# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 08:40:13 2018

@author: rstreet
"""

import os
import sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import glob

def process_image_set():
    
    dir_path = get_args()
    
    file_list = glob.glob(os.path.join(dir_path,'*.fits'))
    
    file_list.sort()
    
    keys = [ 'WMSHUMID', 'WMSTEMP', 'WMSPRES', 'WINDSPEE', 'WINDDIR', 
            'WMSMOIST', 'WMSDEWPT', 'WMSCLOUD', 'WMSSKYBR', 'SKYMAG', 
            'AIRMASS', 'MOONFRAC', 'MOONDIST', 'FOCPOSN', 'FOCTELZP', 'FOCTOFF',
            'L1FWHM','L1MEAN', 'L1SIGMA', 'L1ELLIP']
    
    (image_list,image_data) = get_image_data(file_list,keys,dir_path)
    
    plot_stats_for_sites(image_list,image_data,dir_path,'sky_temp.png', 7, 
    'Temperature [$^{\circ}$ C]', 'Measured sky temperature across all sites')
    
    plot_delta_stats_for_sites(image_list,image_data,dir_path,'delta_sky_mag.png', 9, 8, 
    '$\Delta$ sky brightness [mag]', 'Difference between measured and computed sky brightness')
    
    plot_stats_for_sites(image_list,image_data,dir_path,'focal_position.png', 13, 
    'Focus offset [mm]', 'Actual focal position for all sites')
    plot_stats_for_sites(image_list,image_data,dir_path,'focal_zp.png', 14, 
    'Focus offset [mm]', 'Telescope focal zeropoint')
    plot_stats_for_sites(image_list,image_data,dir_path,'focal_therm_offset.png', 15, 
    'Focus offset [mm]', 'Focal thermal offset')
    
    plot_stats_for_sites(image_list,image_data,dir_path,'moon_frac.png', 11, 
    'Moon fraction', 'Moon fraction for observations')
    
    plot_compare_stats_for_sites(image_list,image_data,dir_path,'moon_distance.png', 12, 8, 
    'Moon separation [$^{\circ}$]','Sky brightness [mag]', 
    'Measured sky brightness as a function of lunar separation')
    
    plot_stats_for_sites(image_list,image_data,dir_path,'mean_sky.png', 17, 
    'Mean sky value [ADU]', 'Average sky background data for all sites')
    
    plot_stats_for_sites(image_list,image_data,dir_path,'mean_sky_sigma.png', 18, 
    'Mean sky sigma [ADU]', 'Average sky background deviation for all sites')
    
    plot_stats_for_sites(image_list,image_data,dir_path,'mean_fwhm.png', 16, 
    'Mean FWHM [arcsec]', 'Average FWHM',exclude_no_data=True)
    
    plot_stats_for_sites(image_list,image_data,dir_path,'mean_ellipticity.png', 19, 
    'Mean ellipticity', 'Average ellipticity',exclude_no_data=True)
    
def plot_stats_for_sites(image_list,image_data,dir_path,plot_file,key_index,
                         ylabel,title,exclude_no_data=False,bad_images=None):
    """Function to plot the values of cloud / sky temperture"""
    
    sites = { 'lsc': 'm', 'cpt': 'b', 'coj': 'c' }
    
    fig = plt.figure(1)
    
    for s in sites.keys():
        
        ydata = []
        bad_xdata = []
        bad_ydata = []
        
        for i in range(0,len(image_list),1):
            
            if image_list[i,1] == s:
                
                if exclude_no_data and image_data[i,key_index] != -99.999:
                    
                    ydata.append(image_data[i,key_index])
                    
                    if bad_images != None and image_list[i,0] in bad_images.keys():
                        
                        bad_xdata.append(i)
                        bad_ydata.append(image_data[i,key_index])

                elif not exclude_no_data:
                    
                    ydata.append(image_data[i,key_index])
                    
                    if bad_images != None and image_list[i,0] in bad_images.keys():
                        
                        bad_xdata.append(i)
                        bad_ydata.append(image_data[i,key_index])
                        
        xdata = range(0,len(ydata),1)
        
        plt.plot(xdata,np.array(ydata),sites[s]+'.',label=s)
        
        if bad_images != None:
            
            plt.plot(bad_xdata,np.array(bad_ydata),'rx')
            
    plt.xlabel('Image index')
    
    plt.ylabel(ylabel)
    
    plt.title(title)
    
    plt.legend()
        
    plt.savefig( os.path.join(dir_path,plot_file) )

    plt.close(1)


def plot_delta_stats_for_sites(image_list,image_data,dir_path,plot_file,
                               key_index1,key_index2,ylabel,title,bad_images=None):
    """Function to plot the values of cloud / sky temperture"""
    
    sites = { 'lsc': 'm', 'cpt': 'b', 'coj': 'c' }
    
    fig = plt.figure(1)
    
    for s in sites.keys():
        
        xdata = []
        ydata = []
        bad_xdata = []
        bad_ydata = []
        
        for i in range(0,len(image_list),1):
            
            if image_list[i,1] == s:
                
                xdata.append(i)
                ydata.append(image_data[i,key_index1]-image_data[i,key_index2])
        
                if bad_images != None and image_list[i,0] in bad_images.keys():
                        
                    bad_xdata.append(i)
                    bad_ydata.append(image_data[i,key_index1]-image_data[i,key_index2])
        
        plt.plot(np.array(xdata),np.array(ydata),sites[s]+'.',label=s)
        
        if bad_images != None:
            
            plt.plot(bad_xdata,np.array(bad_ydata),'rx')
            
    plt.xlabel('Image index')
    
    plt.ylabel(ylabel)
    
    plt.title(title)
    
    plt.legend()
        
    plt.savefig( os.path.join(dir_path,plot_file) )

    plt.close(1)
    
def plot_compare_stats_for_sites(image_list,image_data,dir_path,plot_file,
                               key_index1,key_index2,xlabel,ylabel,title,bad_images=None):
    """Function to plot one statistic against another"""
    
    sites = { 'lsc': 'm', 'cpt': 'b', 'coj': 'c' }
    
    fig = plt.figure(1)
    
    for s in sites.keys():
        
        xdata = []
        ydata = []
        bad_xdata = []
        bad_ydata = []
        
        for i in range(0,len(image_list),1):
            
            if image_list[i,1] == s:
                
                xdata.append(image_data[i,key_index1])
                ydata.append(image_data[i,key_index2])
        
                if bad_images != None and image_list[i,0] in bad_images.keys():
                        
                    bad_xdata.append(image_data[i,key_index1])
                    bad_ydata.append(image_data[i,key_index2])
                    
        plt.plot(xdata,ydata,sites[s]+'.',label=s)
        
        if bad_images != None:
            
            plt.plot(bad_xdata,np.array(bad_ydata),'rx')
            
    plt.xlabel(xlabel)
    
    plt.ylabel(ylabel)
    
    plt.title(title)
    
    plt.legend()
        
    plt.savefig( os.path.join(dir_path,plot_file) )

    plt.close(1)
    
def get_image_data(file_list,keys,dir_path):
    """Function to extract the image header keyword parameters relating 
    to weather and other factors influencing data quality from a set
    of images, and output them to file and as a np.array.
    """

    image_list = []
    image_data = []
    
    foutput = open(os.path.join(dir_path,'image_qc.dat'),'w')
    foutput.write('# Image  Site  Filter '+' '.join(keys)+'\n')
    
    for image_file in file_list:
        
        (data,f) = parse_image_header(image_file,keys)
        
        site = os.path.basename(image_file)[0:3]
        
        tstr = os.path.basename(image_file)+' '+site+' '+f+' '
        
        for i in range(0,len(keys),1):
            
            tstr += str(data[i])+' '
        
        foutput.write(tstr+'\n')
        
        if len(data) > 0:
            
            image_list.append([os.path.basename(image_file),site])
            image_data.append(data)
                
    foutput.close()
    
    image_list = np.array(image_list)
    image_data = np.array(image_data)
    
    return image_list,image_data
    
def parse_image_header(image_file,keys):
    """Function to extract the image header keyword parameters relating 
    to weather and other factors influencing data quality.
    
    Assumes LCO-standard header format.    
    """
    
    data = []
    f = None
    
    if os.path.isfile(image_file):
        
        header = fits.getheader(image_file)
        
        f = header['FILTER']
        
        for key in keys:
            
            try:
            
                data.append( float(header[key]) )
            
            except ValueError:
                
                data.append( -99.999 )
                
            except KeyError:
                
                data.append( -99.999 )
                
    return data,f
    
def get_args():
    """Function to obtain or prompt for commandline arguments"""
    
    if len(sys.argv) == 1:
        
        dir_path = raw_input('Please enter the path to the image directory: ')
    
    else:
        
        dir_path = sys.argv[1]
        
    return dir_path
    

if __name__ == '__main__':
    
    process_image_set()
    