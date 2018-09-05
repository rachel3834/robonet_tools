# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 16:59:19 2018

@author: rstreet
"""

from os import path
from sys import argv
import numpy as np
import json

class Star():
    """Class to describe the photometric properties of a single stellar object"""
    
    def __init__(self,dm=None):
        self.name = None
        self.Mg = None
        self.Mr = None
        self.Mi = None
        self.Mz = None
        self.MJ = None
        self.MH = None
        self.MKs = None
        self.W149 = None
        self.Z087 = None
        self.distance_modulus = dm
        self.ZP = 25.0
        
    def calculate_colours(self):
        
        colours = [('g','r'), ('r','i'), ('g','i'), 
                   ('J','H'), ('H','Ks'), ('J','Ks')]
        
        for band1,band2 in colours:
            
            col = band1+band2
            
            try:
                m1 = getattr(self,'M'+band1)
                m2 = getattr(self,'M'+band2)
                
                setattr(self,col,(m1-m2))
                                
            except AttributeError:
                pass
    
    def calculate_apparent_magnitudes(self):
        
        for band in ['g','r','i','z','J','H','Ks','W149','Z087']:
            
            m = getattr(self,'M'+band)
            m = m + self.distance_modulus
            setattr(self,'m'+band,m)
    
    def summary(self):
        return self.name+' '+str(round(self.Mg,3))+' '+str(round(self.Mr,3))+' '+str(round(self.Mi,3))+' '+str(round(self.Mz,3))+\
                         ' '+str(round(self.MJ,3))+' '+str(round(self.MH,3))+' '+str(round(self.MKs,3))
        
    def summary_app(self):
        return self.name+' '+str(round(self.mg,3))+' '+str(round(self.mr,3))+' '+str(round(self.mi,3))+' '+str(round(self.mz,3))+\
                         ' '+str(round(self.mW149,3))+' '+str(round(self.mZ087,3))
        #                 ' '+str(round(self.mJ,3))+' '+str(round(self.mH,3))+' '+str(round(self.mKs,3))+\
    
    def estimate_wfirst_abs_magnitudes(self):
        
        z = getattr(self,'Mz')
        j = getattr(self,'MJ')
        h = getattr(self,'MH')
        fj = mag_to_flux(j,self.ZP)
        fh = mag_to_flux(h,self.ZP)
        
        fcomb = fj + fj
        mcomb = flux_to_mag(fcomb,self.ZP)
            
        setattr(self,'MW149',mcomb)
        setattr(self,'MZ087',z)
        
class BinaryStar():
    """Class to described the photometric properties of a binary star"""
    
    def __init__(self,dm=None):
        self.star1 = Star()
        self.star2 = Star()
        self.ZP = self.star1.ZP
        self.distance_modulus = dm
        
    def summary(self):
        
        print 'Star      M_g      Mr     Mi     MJ      MH     MKs'
        print self.star1.summary()
        print self.star2.summary()
        print 'DM = ',self.distance_modulus
        
    def read_stellar_properties(self, file_path):
        
        if path.isfile(file_path) == False:
            print('Error: Cannot find input file '+file_path)
            exit()
        
        data = json.load(open(file_path,'r'))
        
        for key,value in data['primary']:
                        
            setattr(self.star1, key, float(value))
        
        
        for key,value in data['secondary']:
                        
            setattr(self.star1, key, float(value))
        
        self.distance_modulus = float(data['distance_modulus'])
        
    def calculate_combined_light(self):
        
        for band in ['g','r','i','J','H','Ks']:
            
            m = getattr(self.star1,'M'+band)
            f1 = mag_to_flux(m,self.star1.ZP)
            setattr(self.star1,'f'+band,f1)
            
            m = getattr(self.star2,'M'+band)
            f2 = mag_to_flux(m,self.star1.ZP)
            setattr(self.star2,'f'+band,f2)
            
            fcomb = f1 + f2
            mcomb = flux_to_mag(fcomb,self.ZP)
            
            setattr(self,'f'+band+'_combined',fcomb)
            setattr(self,'M'+band+'_combined',mcomb)
            
            #print 'M'+band+'_combined = ',getattr(self,'M'+band+'_combined')
            
    def calculate_apparent_magnitudes(self):
        
        self.star1.calculate_apparent_magnitudes()
        self.star2.calculate_apparent_magnitudes()
            
        for band in ['g','r','i','J','H','Ks']:
            
            m = getattr(self,'M'+band+'_combined')
            m = m + self.distance_modulus
            setattr(self,'m'+band+'_combined',m)
    
    def calculate_combined_colours(self):
        
        passbands1 = [ 'g', 'r', 'J', 'H', 'J' ]
        passbands2 = [ 'r', 'i', 'H', 'Ks', 'Ks' ]
        
        for i in range(0,len(passbands1),1):
            
            band1 = passbands1[i]
            band2 = passbands2[i]
            
            m1 = getattr(self,'M'+band1+'_combined')
            m2 = getattr(self,'M'+band2+'_combined')
            
            col = m1 - m2
            
            setattr(self,band1+band2+'_combined',col)
            
    def output_table(self, output_path):
                
        band_order = ['g', 'r', 'i', 'J', 'H', 'Ks']
        colour_order = ['gr', 'ri', 'JH', 'HKs', 'JKs']
        
        t = open(output_path, 'w')
    
        t.write('\\begin{table}[h!]\n')
        t.write('\\centering\n')
        t.write('\\caption{Predicted photometric properties of the binary star.} \label{tab:binaryphot}\n')
        t.write('\\begin{tabular}{lccc}\n')
        t.write('\\hline\n')
        t.write('\\hline\n')
        t.write('Quantity		& K-star (M=0.90$M_{\odot}$) & M-star (M=0.25$M_{\odot}$) & Combined light\\\\\n')
        
        for key in band_order:
            
            m1 = getattr(self.star1,'M'+key)
            m2 = getattr(self.star2,'M'+key)
            mcomb = getattr(self,'M'+key+'_combined')
            
            t.write('$M_{'+key+'}$ [mag]         & '+str(round(m1,3))+'   & '+str(round(m2,3))+'  & '+str(round(mcomb,3))+'\\\\\n')
        
        for key in colour_order:
            
            col1 = getattr(self.star1,key)
            col2 = getattr(self.star2,key)
            colcomb = getattr(self,key+'_combined')
            
            t.write('$('+key+')$ [mag]            & '+str(round(col1,3))+'   & '+str(round(col2,3))+'  & '+str(round(colcomb,3))+'\\\\\n')
        
        
        for key in band_order:
            
            m1 = getattr(self.star1,'m'+key)
            m2 = getattr(self.star2,'m'+key)
            mcomb = getattr(self,'m'+key+'_combined')
            
            t.write('$m_{'+key+'}$ [mag]          & '+str(round(m1,3))+'   & '+str(round(m2,3))+'  & '+str(round(mcomb,3))+'\\\\\n')
        
        t.write('\\hline\n')
        t.write('\\end{tabular}\n')
        t.write('\\end{table}\n')
    
        t.close()

def calc_combined_light():
    """Function to calculate the combined light from the photometric properties 
    of two separate stars"""
    
    if len(argv) > 1:
        file_path = argv[1]
        output_table = argv[2]
    else:
        primary_data = raw_input('Please enter the path to the input file: ')
        output_table = raw_input('Please enter the path to the output table: ')
    
    binary = BinaryStar()
    
    binary.read_stellar_properties(file_path)

    binary.summary()
    
    binary.calculate_combined_light()
    
    binary.calculate_apparent_magnitudes()
    
    binary.calculate_colours()
    
    binary.output_table(output_table)
    
def mag_to_flux(mag, ZP):
    """Function to convert magnitudes into flux units"""
    
    flux = 10**((mag-ZP)/-2.5)
    
    #ferr = mag_err/(2.5*np.log10(np.e)) * flux
    
    return flux

    
def flux_to_mag(flux, ZP):
    """Function to convert the flux of a star from its fitted PSF model 
    and its uncertainty onto the magnitude scale.
    
    :param float flux: Total star flux
    
    Returns:
    
    :param float mag: Measured star magnitude
    """
    
    mag = ZP - 2.5 * np.log10(flux)
    
    return mag
    

if __name__ == '__main__':
    
    calc_combined_light()
    