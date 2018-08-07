# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 16:59:19 2018

@author: rstreet
"""

from os import path
from sys import argv
import numpy as np

class BinaryStar():
    """Class to described the photometric properties of a binary star"""
    
    def __init__(self):
        self.Mg_star1 = None
        self.Mr_star1 = None
        self.Mi_star1 = None
        self.MJ_star1 = None
        self.MH_star1 = None
        self.MKs_star1 = None
        
        self.Mg_star2 = None
        self.Mr_star2 = None
        self.Mi_star2 = None
        self.MJ_star2 = None
        self.MH_star2 = None
        self.MKs_star2 = None
        
        self.distance_modulus = None

        self.ZP = 25.0
    
    def summary(self):
        
        print 'Star      M_g      Mr     Mi     MJ      MH     MKs'
        print 'Star 1 ', self.Mg_star1, self.Mr_star1, self.Mi_star1, self.MJ_star1, self.MH_star1, self.MKs_star1
        print 'Star 2 ', self.Mg_star2, self.Mr_star2, self.Mi_star2, self.MJ_star2, self.MH_star2, self.MKs_star2
        print 'DM = ',self.distance_modulus
        
    def read_stellar_properties(self, file_path):
        
        if path.isfile(file_path) == False:
            print('Error: Cannot find input file '+file_path)
            exit()
            
        lines = open(file_path,'r').readlines()
        
        for l in lines:
            (key, value) = l.replace('\n','').split()
            
            setattr(self, key, float(value))
    
    def calculate_combined_light(self):
        
        for band in ['g','r','i','J','H','Ks']:
            
            m = getattr(self,'M'+band+'_star1')
            f1 = mag_to_flux(m,self.ZP)
            setattr(self,'f'+band+'_star1',f1)
            
            m = getattr(self,'M'+band+'_star2')
            f2 = mag_to_flux(m,self.ZP)
            setattr(self,'f'+band+'_star2',f2)
            
            fcomb = f1 + f2
            mcomb = flux_to_mag(fcomb,self.ZP)
            
            setattr(self,'f'+band+'_combined',fcomb)
            setattr(self,'M'+band+'_combined',mcomb)
            
            print 'M'+band+'_combined = ',getattr(self,'M'+band+'_combined')
            
    def calculate_apparent_magnitudes(self):
        
        for band in ['g','r','i','J','H','Ks']:
            
            m = getattr(self,'M'+band+'_star1')
            m = m + self.distance_modulus
            setattr(self,'m'+band+'_star1',m)
            
            m = getattr(self,'M'+band+'_star2')
            m = m + self.distance_modulus
            setattr(self,'m'+band+'_star2',m)
            
            m = getattr(self,'M'+band+'_combined')
            m = m + self.distance_modulus
            setattr(self,'m'+band+'_combined',m)
    
    def calculate_colours(self):
        
        passbands1 = [ 'g', 'r', 'J', 'H', 'J' ]
        passbands2 = [ 'r', 'i', 'H', 'Ks', 'Ks' ]
        
        for i in range(0,len(passbands1),1):
            
            band1 = passbands1[i]
            band2 = passbands2[i]
            
            for s in ['star1', 'star2', 'combined']:
                m1 = getattr(self,'M'+band1+'_'+s)
                m2 = getattr(self,'M'+band2+'_'+s)
                
                col = m1 - m2
                
                setattr(self,band1+'-'+band2+'_'+s,col)
            
    def output_table(self, output_path):
                
        band_order = ['g', 'r', 'i', 'J', 'H', 'Ks']
        colour_order = ['g-r', 'r-i', 'J-H', 'H-Ks', 'J-Ks']
        
        t = open(output_path, 'w')
    
        t.write('\\begin{table}[h!]\n')
        t.write('\\centering\n')
        t.write('\\caption{Predicted photometric properties of the binary star.} \label{tab:binaryphot}\n')
        t.write('\\begin{tabular}{lccc}\n')
        t.write('\\hline\n')
        t.write('\\hline\n')
        t.write('Quantity		& K-star (M=0.90$M_{\odot}$) & M-star (M=0.25$M_{\odot}$) & Combined light\\\\\n')
        
        for key in band_order:
            
            m1 = getattr(self,'M'+key+'_star1')
            m2 = getattr(self,'M'+key+'_star2')
            mcomb = getattr(self,'M'+key+'_combined')
            
            t.write('$M_{'+key+'}$ [mag]         & '+str(round(m1,3))+'   & '+str(round(m2,3))+'  & '+str(round(mcomb,3))+'\\\\\n')
        
        for key in colour_order:
            
            col1 = getattr(self,key+'_star1')
            col2 = getattr(self,key+'_star2')
            colcomb = getattr(self,key+'_combined')
            
            t.write('$('+key+')$ [mag]            & '+str(round(col1,3))+'   & '+str(round(col2,3))+'  & '+str(round(colcomb,3))+'\\\\\n')
        
        
        for key in band_order:
            
            m1 = getattr(self,'m'+key+'_star1')
            m2 = getattr(self,'m'+key+'_star2')
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
        file_path = raw_input('Please enter the path to the input file: ')
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
    