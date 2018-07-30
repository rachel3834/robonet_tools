# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 16:28:17 2018

@author: rstreet
"""

import jester_phot_transforms
import bilir_phot_transforms
import numpy as np

class Star:
    """Class describing the photometric parameters of a single stellar object"""
    
    def __init__(self):
        
        self.g = None
        self.sig_g = None
        self.r = None
        self.sig_r = None
        self.i = None
        self.sig_i = None
        
        self.gr = None
        self.sig_gr = None
        self.ri = None
        self.sig_ri = None
        
        self.B = None
        self.sig_B = None
        self.V = None
        self.sig_V = None
        self.R = None
        self.sig_R = None
        self.I = None
        self.sig_I = None
        self.BV = None
        self.sig_BV = None
        self.RI = None
        self.sig_RI = None
        self.VI = None
        self.sig_VI = None
    
    def compute_colours(self,use_inst=True,use_cal=False):
        """Method to calculate the star's colours and colour uncertainties, 
        given measurements in 3 passbands.  
        use_inst=True will apply this calculation to the input magnitudes
        use_cal=True will apply it to the extinction-corrected magnitudes
        """
        
        if use_inst:
            
            (self.gr,self.sig_gr) = self.calc_colour(self.g,self.sig_g,self.r,self.sig_r)
            (self.ri,self.sig_ri) = self.calc_colour(self.r,self.sig_r,self.i,self.sig_i)
        
        if use_cal:
            
            (self.gr_0,self.sig_gr_0) = self.calc_colour(self.g_0,self.sig_g_0,self.r_0,self.sig_r_0)
            (self.ri_0,self.sig_ri_0) = self.calc_colour(self.r_0,self.sig_r_0,self.i_0,self.sig_i_0)
                                    
    def calc_colour(self,col1, sig_col1, col2, sig_col2):
        
        col = col1 - col2
        sig = np.sqrt( (sig_col1*sig_col1)  + \
                        (sig_col2*sig_col2) )
        return col, sig
    
    def summary(self,show_mags=True, show_cal=False, show_colours=False, 
                johnsons=False):
        
        output = ''
        
        if show_mags:
            output = 'g_meas = '+str(round(self.g,3))+' +/- '+str(round(self.sig_g,3))+\
                ' r_meas = '+str(round(self.r,3))+' +/- '+str(round(self.sig_r,3))+\
                ' i_meas = '+str(round(self.i,3))+' +/- '+str(round(self.sig_i,3))
        
        elif show_cal:
            output = 'g_0 = '+str(round(self.g_0,3))+' +/- '+str(round(self.sig_g_0,3))+\
                ' r_0 = '+str(round(self.r_0,3))+' +/- '+str(round(self.sig_r_0,3))+\
                ' i_0 = '+str(round(self.i_0,3))+' +/- '+str(round(self.sig_i_0,3))
            
        elif show_colours:
            output = '(g-r)_meas = '+str(round(self.gr,3))+' +/- '+str(round(self.sig_gr,3))+\
                '\n(r-i)_meas = '+str(round(self.ri,3))+' +/- '+str(round(self.sig_ri,3))
                
        elif show_colours and show_cal:
            output = '(g-r)_0 = '+str(round(self.gr_0,3))+' +/- '+str(round(self.sig_gr_0,3))+\
                '\n(r-i)_0 = '+str(round(self.ri_0,3))+' +/- '+str(round(self.sig_ri_0,3))
        
        elif johnsons:
            
            if self.VR != None and self.RI != None:
                output = '(V-R)_inst = '+str(self.VR)+' +/- '+str(self.sig_VR)+'mag\n'+\
                         '(Rc-Ic)_inst = '+str(self.RI)+' +/- '+str(self.sig_RI)+'mag'
    
            if self.V != None and self.BV != None:
            
                output = 'V_inst = '+str(self.V)+' +/- '+str(self.sig_V)+'mag\n'+\
                         '(B-V)_inst = '+str(self.BV)+' +/- '+str(self.sig_BV)+'mag'
    
            if self.V != None and self.VR != None:
                    
                output = 'Derived target instrumental colours and magnitudes:\n'+\
                         'R_inst = '+str(self.R)+' +/- '+str(self.sig_R)+'mag'+\
                         'I_inst = '+str(self.I)+' +/- '+str(self.sig_I)+'mag'+\
                         '(V-I)_inst = '+str(self.VI)+' +/- '+str(self.sig_VI)+'mag'
            
        return output
    
    def transform_to_JohnsonCousins(self):
        
        if self.ri != None:
            
            target_phot = jester_phot_transforms.transform_SDSS_to_JohnsonCousins(ri=self.ri, 
                                                                                  sigri=self.sig_ri)
            self.VR = target_phot['V-R']
            self.sig_VR = target_phot['sigVR']
            self.RI = target_phot['Rc-Ic']
            self.sig_RI = target_phot['sigRI']
        
        if self.g != None and self.gr != None:
            
            target_phot = jester_phot_transforms.transform_SDSS_to_JohnsonCousins(g=self.g,
                                                                                  sigg=self.sig_g,
                                                                                  gr=self.gr, 
                                                                                  siggr=self.sig_gr)
            
            self.V = target_phot['V']
            self.sig_V = target_phot['sigV']
            self.BV = target_phot['B-V']
            self.sig_BV = target_phot['sigBV']
            
        
        self.R = self.V - self.VR
        self.sig_R = np.sqrt( (self.sig_V*self.sig_V) + (self.sig_VR*self.sig_VR) ) 
        self.I = self.R - self.RI
        self.sig_I = np.sqrt( (self.sig_R*self.sig_R) + (self.sig_RI*self.sig_RI) ) 
        self.VI = self.V - self.I
        self.sig_VI = np.sqrt( (self.sig_V*self.sig_V) + (self.sig_I*self.sig_I) ) 

    def transform_2MASS_to_SDSS(self):
        
        (self.gr_0, self.sig_gr_0, self.ri_0, self.sig_ri_0) = bilir_phot_transforms.transform_2MASS_to_SDSS(JH=self.JH_0, HK=self.HK_0, MH=None)
    
    def calibrate_phot_properties(self, RC, log=None):
        """Function to calculate the de-reddened and extinction-corrected 
        photometric properties of the Star
        """
        
        self.g_0 = self.g - RC.A_g
        self.sig_g_0 = np.sqrt( (self.sig_g*self.sig_g) + (RC.sig_A_g*RC.sig_A_g) )
        self.r_0 = self.r - RC.A_r
        self.sig_r_0 = np.sqrt( (self.sig_r*self.sig_r) + (RC.sig_A_r*RC.sig_A_r) )
        self.i_0 = self.i - RC.A_i
        self.sig_i_0 = np.sqrt( (self.sig_i*self.sig_i) + (RC.sig_A_i*RC.sig_A_i) )
        self.gr_0 = self.gr - RC.Egr
        self.sig_gr_0 = np.sqrt( (self.sig_gr*self.sig_gr) + (RC.sig_Egr*RC.sig_Egr) )
        self.ri_0 = self.ri - RC.Eri
        self.sig_ri_0 = np.sqrt( (self.sig_ri*self.sig_ri) + (RC.sig_Eri*RC.sig_Eri) )

        
        output = '\nSource star extinction-corrected magnitudes and de-reddened colours:'
        output += 'g_S,0 = '+str(self.g_0)+' +/- '+str(self.sig_g_0)+'\n'
        output += 'r_S,0 = '+str(self.r_0)+' +/- '+str(self.sig_r_0)+'\n'
        output += 'i_S,0 = '+str(self.i_0)+' +/- '+str(self.sig_i_0)+'\n'
        output += '(g-r)_S,0 = '+str(self.gr_0)+' +/- '+str(self.sig_gr_0)+'\n'
        output += '(r-i)_S,0 = '+str(self.ri_0)+' +/- '+str(self.sig_ri_0)
        
        if log != None:
            log.info(output)
        else:
            print(output)
            
        phot = jester_phot_transforms.transform_SDSS_to_JohnsonCousins(ri=self.ri_0, 
                                                                       sigri=self.sig_ri_0)
        
        self.VR_0 = phot['V-R']
        self.sig_VR_0 = phot['sigVR']
        self.RI_0 = phot['Rc-Ic']
        self.sig_RI_0 = phot['sigRI']
    
        output = '\n(V-R)_S,0 = '+str(self.VR_0)+' +/- '+str(self.sig_VR_0)+'mag\n'
        output += '(R-I)_S,0 = '+str(self.RI_0)+' +/- '+str(self.sig_RI_0)+'mag\n'
        
        phot = jester_phot_transforms.transform_SDSS_to_JohnsonCousins(g=self.g_0,
                                                                       sigg=self.sig_g_0,
                                                                       gr=self.gr_0, 
                                                                       siggr=self.sig_gr_0)
        self.V_0 = phot['V']
        self.sig_V_0 = phot['sigV']
        self.BV_0 = phot['B-V']
        self.sig_BV_0 = phot['sigBV']
    
        output += '\n(B-V)_S,0 = '+str(self.BV_0)+' +/- '+str(self.sig_BV_0)+'mag\n'
        output += 'V_S,0 = '+str(self.V_0)+' +/- '+str(self.sig_V_0)+'mag'
    
        if log != None:
            log.info(output)
        else:
            print(output)
            