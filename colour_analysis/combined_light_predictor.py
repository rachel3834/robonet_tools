# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 16:59:19 2018

@author: rstreet
"""

from os import path
from sys import argv
import numpy as np
import json
import jester_phot_transforms

class Star():
    """Class to describe the photometric properties of a single stellar object"""
    
    def __init__(self,dm=None):
        self.name = None
        self.MV = None
        self.Mg = None
        self.Mr = None
        self.Mi = None
        self.Mz = None
        self.MJ = None
        self.MH = None
        self.MKs = None
        self.W149 = None
        self.Z087 = None
        self.mV = None
        self.distance_modulus = dm
        self.ZP = 25.0
        
    def calculate_colours(self):
        
        colours = [('g','r'), ('r','i'), ('g','i'), 
                   ('J','H'), ('H','Ks'), ('J','Ks'),('B','V')]
        
        for band1,band2 in colours:
            
            col = band1+band2
            
            try:
                m1 = getattr(self,'M'+band1)
                m2 = getattr(self,'M'+band2)
                
                setattr(self,col,(m1-m2))
                                
            except AttributeError:
                pass
    
    def calculate_apparent_magnitudes(self):
        
        for band in ['B','V','g','r','i','z','J','H','Ks','W149','Z087']:
            
            try:
                m = getattr(self,'M'+band)
                
                if m != None:
                    m = m + self.distance_modulus
                    setattr(self,'m'+band,m)
            
            except AttributeError:
                pass
        
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
        
    def apply_extinction(self,extinction):
        
        for f in ['V', 'g', 'r', 'i']:
            
            if 'A'+f in extinction.keys() and getattr(self,'m'+f) != None:

                m = getattr(self,'m'+f)
                
                m = m + extinction['A'+f]
                
                setattr(self,'m'+f+'_corr',m)
                
        if 'EBV' in extinction.keys():
            self.BV_corr = self.BV + extinction['EBV']
        
        c1 = ['g', 'r']
        c2 = ['r', 'i']
        
        for i in range(0,2,1):
            
            m1 = getattr(self,'m'+c1[i]+'_corr')
            m2 = getattr(self,'m'+c2[i]+'_corr')
            
            col = m1 - m2
            
            setattr(self, c1[i]+c2[i]+'_corr', col)
            
    def transform_extinc_corr_Johnson_SDSS(self):
        
        results = jester_phot_transforms.transform_JohnsonCousins_to_SDSS(BV=self.BV_corr, 
                                                                          sigBV=0.0,
                                                                          V=self.mV_corr, 
                                                                          sigV=0.0)
        self.mg_corr = results['g']
        self.sig_mg_corr = results['sigg']
        self.mr_corr = results['r']
        self.sig_mr_corr = results['sigr']
        self.gr_corr = results['g-r']
        self.sig_gr_corr = results['siggr']

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
        
        for band in ['B','V','g','r','i','z','J','H','Ks']:
            
            try:
                m1 = getattr(self.star1,'M'+band)
                m2 = getattr(self.star2,'M'+band)
                
                if m1 != None:
                    f1 = mag_to_flux(m1,self.star1.ZP)
                    setattr(self.star1,'f'+band,f1)
                
                if m1 != None:
                    f2 = mag_to_flux(m2,self.star1.ZP)
                    setattr(self.star2,'f'+band,f2)
                
                if m1 != None and m2 != None:
                    fcomb = f1 + f2
                    mcomb = flux_to_mag(fcomb,self.ZP)
                    
                    setattr(self,'f'+band+'_combined',fcomb)
                    setattr(self,'M'+band+'_combined',mcomb)
                
                #print 'M'+band+'_combined = ',getattr(self,'M'+band+'_combined')
            
            except AttributeError:
                pass
            
    def calculate_apparent_magnitudes(self):
        
        self.star1.calculate_apparent_magnitudes()
        self.star2.calculate_apparent_magnitudes()
            
        for band in ['B','V','g','r','i','z','J','H','Ks']:
            
            try:
                m = getattr(self,'M'+band+'_combined')
                
                if m != None:
                    
                    m = m + self.distance_modulus
                    setattr(self,'m'+band+'_combined',m)
                
            except AttributeError:
                pass
            
    def calculate_combined_colours(self):
        
        self.star1.calculate_colours()
        self.star2.calculate_colours()
        
        passbands1 = [ 'B', 'g', 'r', 'J', 'H', 'J' ]
        passbands2 = [ 'V', 'r', 'i', 'H', 'Ks', 'Ks' ]
        
        for i in range(0,len(passbands1),1):
            
            band1 = passbands1[i]
            band2 = passbands2[i]
            
            m1 = getattr(self,'M'+band1+'_combined')
            m2 = getattr(self,'M'+band2+'_combined')
            
            col = m1 - m2
            
            setattr(self,band1+band2+'_combined',col)
    
    def apply_extinction(self,extinction):
        
        for sid in [ 'star1', 'star2' ]:
            
            star = getattr(self,sid)
            
            star.apply_extinction(extinction)
    
            setattr(self,sid,star)
        
        try:
    
            for f in ['V', 'g', 'r', 'i']:
            
                if 'A'+f in extinction.keys():
    
                    m = getattr(self,'m'+f+'_combined')
                    
                    m = m + extinction['A'+f]
                    
                    setattr(self,'m'+f+'_combined_corr',m)
                
            bv = getattr(self,'BV_combined')
            
            bv = bv + extinction['EBV']
            
            setattr(self,'BV_combined_corr',bv)
            
            c1 = ['g', 'r']
            c2 = ['r', 'i']
            
            for i in range(0,2,1):
                
                m1 = getattr(self,'m'+c1[i]+'_combined_corr')
                m2 = getattr(self,'m'+c2[i]+'_combined_corr')
                
                col = m1 - m2
                
                setattr(self, c1[i]+c2[i]+'_combined_corr', col)
        
        except:
            pass
    
    
    def transform_extinc_corr_Johnson_SDSS(self):
        
        results = jester_phot_transforms.transform_JohnsonCousins_to_SDSS(BV=self.BV_combined_corr, 
                                                                          sigBV=0.0,
                                                                          V=self.mV_combined_corr, 
                                                                          sigV=0.0)
        self.mg_combined_corr = results['g']
        self.sig_mg_combined_corr = results['sigg']
        self.mr_combined_corr = results['r']
        self.sig_mr_combined_corr = results['sigr']
        self.gr_combined_corr = results['g-r']
        self.sig_gr_combined_corr = results['siggr']
        
        for sid in [ 'star1', 'star2' ]:
            
            star = getattr(self,sid)
            
            star.transform_extinc_corr_Johnson_SDSS()
            
            setattr(self,sid,star)
            
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
    