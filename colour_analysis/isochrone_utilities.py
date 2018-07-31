# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 09:01:46 2018

@author: rstreet
"""

from os import path
import numpy as np
import matplotlib.pyplot as plt

def analyze_isochrones(gr, ri, isochrone_file, log=None):
    """Function to identify the properties of a star with the given 
    (g-r), (r-i) colour data from PARSEC isochrone models. 
    """
    
    data = read_PARSEC_table(isochrone_file)
    
    (star_data, teff, sig_teff) = find_isochrone_nearest_entries(data, gr, ri, log=log)

    #plot_isochrones(data)
    
    #overlay_isochrones(None,data)
    
    return star_data, teff, sig_teff
    
def find_isochrone_nearest_entries(data, gr, ri, log=None):
    """Function to interpolate over the isochrone data, and return the 
    absolute g, r, i magnitudes closest to the target values."""
    
    tol = 0.005
    
    iso_gr = data[:,4] - data[:,5]
    iso_ri = data[:,5] - data[:,6]
    
    dgr = iso_gr - gr
    dri = iso_ri - ri
    dcol = np.sqrt(dgr*dgr + dri*dri)
    
    idx1 = np.where(dcol >= dcol.min()-tol)[0]
    idx2 = np.where(dcol <= dcol.min()+tol)[0]
    idx = list(set(idx1).intersection(set(idx2)))
    
    sort_idx = dcol[idx].argsort()
    sort_idx = np.array(idx)[sort_idx]
    
    output = '\nAnalysing isochrones\n'
    output += 'Closest entries matching target colours: (g-r)='+str(gr)+\
                                                ' and (r-i)='+str(ri)+'\n'
    output += 'Age[Gyr]  Teff[K]  log(g) (g-r)_iso  (r-i)_iso  (g-r)_t  (r-i)_t  dist  gmag   rmag   imag\n'
    
    for i in sort_idx:
        output += str(round(data[i,0]/1e9,4))+' '+str(round(10**(data[i,2]),1))+\
                    ' '+str(data[i,3])+' '+str(round(iso_gr[i],3))+' '+\
                    str(round(iso_ri[i],3))+' '+str(round(dcol[i],4))+' '+\
                    str(round(data[i,4],3))+' '+str(round(data[i,5],3))+' '+\
                    str(round(data[i,6],3))+'\n'
    
    teff = (10**data[sort_idx,2]).mean()
    teff_min = (10**data[sort_idx,2]).min()
    teff_max = (10**data[sort_idx,2]).max()
    sig_teff = (teff_max - teff_min)/2.0
    
    output += 'Range of Teff = '+str(round(teff_min,1))+\
                                ' to '+\
                                str(round(teff_max,1))+'\n'
    
    output += 'Teff of target = '+str(round(teff,1))+' +/- '+\
                                str(round(sig_teff,1))+'\n'
    
    logg = data[sort_idx,3].mean()
    logg_min = data[sort_idx,3].min()
    logg_max = data[sort_idx,3].max()
    sig_logg = (teff_max - teff_min)/2.0
    
    print data[:,3]
    
    output += 'Range of log(g) = '+str(round(logg_min,1))+\
                                ' to '+\
                                str(round(logg_max,1))+'\n'
    
    output += 'log(g) of target = '+str(round(logg,1))+' +/- '+\
                                str(round(sig_logg,1))+'\n'
               
    if log != None:
        log.info(output)
    else:
        print(output)
    
    return data[sort_idx], teff, sig_teff
    
def read_PARSEC_table(table_file):
    """Function to read the relevant data from an isochrone PARSEC table, 
    extracted from the CMD interface at:
    http://stev.oapd.inaf.it/cgi-bin/cmd
    with SDSS colour data selected.
    
    Returns the following columns in a numpy array:
    2 Age
    5 log(L)
    6 log(Te)
    7 log(g)
    25 gmag
    26 rmag
    27 imag
    """
    
    if path.isfile(table_file) == False:
        
        print('ERROR: Cannot find isochrone table file '+table_file)
        exit()
        
    lines = open(table_file, 'r').readlines()
    
    data = []
    
    for l in lines:
        
        if l[0:1] != '#':
            
            entries = l.replace('\n','').split()
            
            col_data = [ ]
            
            for c in [ 1,4,5,6,24,25,26 ]:
                
                col_data.append(float(entries[c]))
            
            data.append(col_data)
    
    return np.array(data)

def plot_isochrones(data):
    """Function to plot isochrones by themselves"""
    
    fig = plt.figure(1,(10,10))
        
    
    fig = overlay_isochrones(fig,data)
    
    plt.xlabel('(g-r) [mag]')
    plt.ylabel('(r-i) [mag]')
    
    plt.axis([-1.0,2.0,-0.5,1.0])
    
    plt.legend()
    
    plt.savefig('isochrone_plot.png')
    
def overlay_isochrones(fig,data,n_steps=2,label_plot=True):
    """Function to add isochrones to an existing (g-r).vs.(r-i) colour-colour
    diagram"""
    
    ages = np.unique(data[:,0])
        
    r = 10.0/185.0
    b = 10.0/185.0
    ginit = 245.0
    gincr = ginit/(float(len(ages))/float(n_steps))
        
    for i,j in enumerate(range(0,len(ages),n_steps)):
        
        a = ages[j]
        
        idx = np.where(data[:,0] == a)
        
        iso_gr = data[idx,4] - data[idx,5]
        iso_ri = data[idx,5] - data[idx,6]
        
        sort = iso_gr.argsort()[0]
                        
        g = (ginit - (i*gincr))/255.0
        
        if label_plot:
            
            plt.scatter(iso_gr, iso_ri, color=(r,g,b), marker='v', s=4, alpha=0.5, label="{0:.4e}".format(a)+' yr')
            
        else:
            
            plt.scatter(iso_gr[:,sort], iso_ri[:,sort], color=(r,g,b), s=4, marker='v', alpha=0.5)
    
    return fig
    
if __name__ == '__main__':
    
    isochrone_file = '/Users/rstreet/ROMEREA/2018/OGLE-2018-BLG-0022/colour_analysis2/isochrones/solar_metallicity_age_range_isochrones.iso'
    gr = 0.88
    ri = 0.33
    
    (star_data, teff, sig_teff) = analyze_isochrones(gr,ri,isochrone_file)
    