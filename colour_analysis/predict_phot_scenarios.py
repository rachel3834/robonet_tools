# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 15:27:52 2018

@author: rstreet
"""

from os import path
from sys import argv
import combined_light_predictor
import jester_phot_transforms

def predict_photometry(output_path,dm,extinction):
    """Driver code to set up scenarios and generate output"""
    
    kstar = Kdwarf(name='K star', dm=dm)
    kstar.calculate_apparent_magnitudes()
    kstar.apply_extinction(extinction)
    #kstar.transform_extinc_corr_Johnson_SDSS()
    
    mstar = Mdwarf(name='M star', dm=dm)
    mstar.calculate_apparent_magnitudes()
    mstar.apply_extinction(extinction)
    #mstar.transform_extinc_corr_Johnson_SDSS()
    
    wd_hot = WhiteDwarf(100000.0,dm=dm)
    wd_hot.calculate_apparent_magnitudes()
    wd_hot.apply_extinction(extinction)
    #wd_hot.transform_extinc_corr_Johnson_SDSS()
    
    wd_mid = WhiteDwarf(20000.0,dm=dm)
    wd_mid.calculate_apparent_magnitudes()
    wd_mid.apply_extinction(extinction)
    #wd_mid.transform_extinc_corr_Johnson_SDSS()
    
    wd_cool = WhiteDwarf(4000.0,dm=dm)
    wd_cool.calculate_apparent_magnitudes()
    wd_cool.apply_extinction(extinction)
    #wd_cool.transform_extinc_corr_Johnson_SDSS()
    
    ms_binary = scenario1(dm,extinction)
    
    k_wd_binary = scenario2(dm,extinction)
    
    m_wd_binary = scenario3(dm,extinction)
    
    output_latex_table(kstar,mstar,wd_hot,wd_mid,wd_cool, 
                       ms_binary, k_wd_binary, m_wd_binary, output_path)
                       

def Kdwarf(name=None,dm=None):
    """Definition of a single, main sequence K-dwarf star of mass 0.88 Msol"""
    
    star = combined_light_predictor.Star()
    
    star.name = name
    star.Mg = 5.749
    star.Mr = 5.199
    star.Mi = 5.058
    star.MB = 6.368
    star.MV = 5.584
    star.MJ = 4.329
    star.MH = 3.942
    star.MKs = 3.897
    
    star.distance_modulus = dm
    
    star.calculate_colours()
        
    return star

def Mdwarf(name=None,dm=None):
    """Definition of a single, main sequence M-dwarf star of mass 0.25 Msol"""
    
    star = combined_light_predictor.Star()
    
    star.name = name
    star.Mg = 13.650
    star.Mr = 12.081
    star.Mi = 10.714
    star.MB = 14.350
    star.MV = 12.676
    star.MJ = 8.398
    star.MH = 7.910
    star.MKs = 7.676
    
    star.distance_modulus = dm
    
    star.calculate_colours()
    
    return star

def WhiteDwarf(teff, dm=None):
    """Definition of a single White Dwarf star, with parameters determined
    by the t_eff provided
    teff one of { 4000.0, 20000.0, 100000.0 }
    Photometry from Bergeron 1995 PASP 107 1047.
    """

    star = combined_light_predictor.Star()
    
    star.distance_modulus = dm
    
    if teff == 4000.0:
        
        star.name = 'WD (4000K)'
        star.MV = 16.105
        star.MB = 17.081
        star.BV = 0.976
        star.VI = 1.249
        star.RI = 0.620
        star.VKs = 1.309
        star.JH = -0.152
        star.HKs = -0.126
    
    elif teff == 20000.0:
        
        star.name = 'WD (20000K)'
        star.MV = 10.732
        star.MB = 10.688
        star.BV = -0.044
        star.VI = -0.216
        star.RI = -0.119
        star.VKs = -0.705
        star.JH = -0.063
        star.HKs = -0.106
    
    elif teff == 100000.0:
        
        star.name = 'WD (100000K)'
        star.MV = 8.408
        star.MB = 8.079
        star.BV = -0.329
        star.VI = -0.358
        star.RI = -0.210
        star.VKs = -1.093
        star.JH = -0.131
        star.HKs = -0.134
    
    else:
        
        print('ERROR: No photometry available for a white dwarf of teff='+str(teff)+'K')
        
        return star
        
    star.MI = star.MV - star.VI
    star.MKs = star.MV - star.VKs
    star.MH = star.MKs + star.HKs
    star.MJ = star.MH + star.JH
    star.JKs = star.MJ - star.MKs
    
    phot = jester_phot_transforms.transform_JohnsonCousins_to_SDSS(BV=star.BV, sigBV=0.0, 
                                  RI=star.RI, sigRI=0.0,
                                  V=star.MV, sigV=0.0)
    
    star.gr = phot['g-r']
    star.ri = phot['r-i']
    star.Mg = phot['g']
    star.Mr = phot['r']
    star.Mi = star.Mr - star.ri
    
    return star

    
def scenario1(dm,extinction):
    """Define a main sequence binary consisting of a K-dwarf and an M-dwarf"""

    kstar = Kdwarf(name='K star', dm=dm)
    mstar = Mdwarf(name='M star', dm=dm)
    
    binary = combined_light_predictor.BinaryStar()
    binary.distance_modulus = dm
    
    binary.star1 = kstar
    binary.star2 = mstar
    
    binary.calculate_combined_light()
    
    binary.calculate_apparent_magnitudes()
    
    binary.calculate_combined_colours()
    
    binary.apply_extinction(extinction)
    
    return binary

def scenario2(dm,extinction):
    """Define a binary consisting of a mid-range temperature white dwarf 
    plus a main sequence K-star"""
    
    kstar = Kdwarf(name='K star', dm=dm)
    wd = WhiteDwarf(20000.0,dm=dm)
    
    binary = combined_light_predictor.BinaryStar()
    binary.distance_modulus = dm
    
    binary.star1 = kstar
    binary.star2 = wd
    
    binary.calculate_combined_light()
    
    binary.calculate_apparent_magnitudes()
    
    binary.calculate_combined_colours()

    binary.apply_extinction(extinction)
    
    return binary
   

def scenario3(dm,extinction):
    """Define a binary consisting of a mid-range temperature white dwarf 
    plus a main sequence M-star"""
    
    mstar = Mdwarf(name='M star', dm=dm)
    wd = WhiteDwarf(20000.0,dm=dm)
    
    binary = combined_light_predictor.BinaryStar()
    binary.distance_modulus = dm
    
    binary.star1 = wd
    binary.star2 = mstar
    
    binary.calculate_combined_light()
    
    binary.calculate_apparent_magnitudes()
    
    binary.calculate_combined_colours()

    binary.apply_extinction(extinction)
    
    return binary

def roundfloat(value,ndp):
    """Function to convert a floating point value to a string, rounding
    the the given number of decimal places in the process. 
    The resulting string will be end zero-padded to ensure a consistent
    number of characters are returned.
    """
    
    output = str( round(value,ndp) )
    
    while len(output.split('.')[1]) < ndp:
        output = output + '0'
    
    return output
    
def output_latex_table(kstar,mstar,wd_hot,wd_mid,wd_cool, 
                       ms_binary, k_wd_binary, m_wd_binary, output_path):
    """Function to output a LaTeX-format table of the stellar photometry"""
    
    t = open(output_path, 'w')
    
    t.write('\\begin{table}[h!]\n')
    t.write('\\centering\n')
    t.write('\\caption{Predicted photometric properties of different lens system scenarios.  Apparent magnitudes are calculated for the measured lens distance without extinction or reddening, except for the bottom section.  White dwarf is abbreviated to WD, main sequence to MS.} \\label{tab:binaryphot}\n')
    t.write('\\begin{tabular}{lcccccccc}\n')
    t.write('\\hline\n')
    t.write('\\hline\n')
    t.write('Quantity	   & K-dwarf & M-dwarf & WD         & WD          & WD            &MS & WD +   & WD +  \\\\\n')
    t.write('\\[mag\\]		   &		 & 	     & 4000\\,K  & 20,000\\,K & 100,000\\,K  & K+M    & K-dwarf   & M-dwarf  \\\\\n')
    t.write('$M_{B}$       & '+roundfloat(kstar.MB,3)+' & '+roundfloat(mstar.MB,3)+' & '+roundfloat(wd_cool.MB,3)+' & '+roundfloat(wd_mid.MB,3)+' & '+roundfloat(wd_hot.MB,3)+' & '+roundfloat(ms_binary.MB_combined,3)+' & '+roundfloat(k_wd_binary.MB_combined,3)+' & '+roundfloat(m_wd_binary.MB_combined,3)+'\\\\\n')
    t.write('$M_{V}$       & '+roundfloat(kstar.MV,3)+' & '+roundfloat(mstar.MV,3)+' & '+roundfloat(wd_cool.MV,3)+' & '+roundfloat(wd_mid.MV,3)+' & '+roundfloat(wd_hot.MV,3)+' & '+roundfloat(ms_binary.MV_combined,3)+' & '+roundfloat(k_wd_binary.MV_combined,3)+' & '+roundfloat(m_wd_binary.MV_combined,3)+'\\\\\n')
    t.write('$M_{g}$       & '+roundfloat(kstar.Mg,3)+' & '+roundfloat(mstar.Mg,3)+' & '+roundfloat(wd_cool.Mg,3)+' & '+roundfloat(wd_mid.Mg,3)+' & '+roundfloat(wd_hot.Mg,3)+' & '+roundfloat(ms_binary.Mg_combined,3)+' & '+roundfloat(k_wd_binary.Mg_combined,3)+' & '+roundfloat(m_wd_binary.Mg_combined,3)+'\\\\\n')
    t.write('$M_{r}$       & '+roundfloat(kstar.Mr,3)+' & '+roundfloat(mstar.Mr,3)+' & '+roundfloat(wd_cool.Mr,3)+' & '+roundfloat(wd_mid.Mr,3)+' & '+roundfloat(wd_hot.Mr,3)+' & '+roundfloat(ms_binary.Mr_combined,3)+' & '+roundfloat(k_wd_binary.Mr_combined,3)+' & '+roundfloat(m_wd_binary.Mr_combined,3)+'\\\\\n')
    t.write('$M_{i}$       & '+roundfloat(kstar.Mi,3)+' & '+roundfloat(mstar.Mi,3)+' & '+roundfloat(wd_cool.Mi,3)+' & '+roundfloat(wd_mid.Mi,3)+' & '+roundfloat(wd_hot.Mi,3)+' & '+roundfloat(ms_binary.Mi_combined,3)+' & '+roundfloat(k_wd_binary.Mi_combined,3)+' & '+roundfloat(m_wd_binary.Mi_combined,3)+'\\\\\n')
    t.write('$M_{J}$       & '+roundfloat(kstar.MJ,3)+' & '+roundfloat(mstar.MJ,3)+' & '+roundfloat(wd_cool.MJ,3)+' & '+roundfloat(wd_mid.MJ,3)+' & '+roundfloat(wd_hot.MJ,3)+' & '+roundfloat(ms_binary.MJ_combined,3)+' & '+roundfloat(k_wd_binary.MJ_combined,3)+' & '+roundfloat(m_wd_binary.MJ_combined,3)+'\\\\\n')
    t.write('$M_{H}$       & '+roundfloat(kstar.MH,3)+' & '+roundfloat(mstar.MH,3)+' & '+roundfloat(wd_cool.MH,3)+' & '+roundfloat(wd_mid.MH,3)+' & '+roundfloat(wd_hot.MH,3)+' & '+roundfloat(ms_binary.MH_combined,3)+' & '+roundfloat(k_wd_binary.MH_combined,3)+' & '+roundfloat(m_wd_binary.MH_combined,3)+'\\\\\n')
    t.write('$M_{Ks}$      & '+roundfloat(kstar.MKs,3)+' & '+roundfloat(mstar.MKs,3)+' & '+roundfloat(wd_cool.MKs,3)+' & '+roundfloat(wd_mid.MKs,3)+' & '+roundfloat(wd_hot.MKs,3)+' & '+roundfloat(ms_binary.MKs_combined,3)+' & '+roundfloat(k_wd_binary.MKs_combined,3)+' & '+roundfloat(m_wd_binary.MKs_combined,3)+'\\\\\n')
    t.write('$(B-V)$       & '+roundfloat(kstar.BV,3)+' & '+roundfloat(mstar.BV,3)+' & '+roundfloat(wd_cool.BV,3)+' & '+roundfloat(wd_mid.BV,3)+' & '+roundfloat(wd_hot.BV,3)+' & '+roundfloat(ms_binary.BV_combined,3)+' & '+roundfloat(k_wd_binary.BV_combined,3)+' & '+roundfloat(m_wd_binary.BV_combined,3)+'\\\\\n')
    t.write('$(g-r)$       & '+roundfloat(kstar.gr,3)+' & '+roundfloat(mstar.gr,3)+' & '+roundfloat(wd_cool.gr,3)+' & '+roundfloat(wd_mid.gr,3)+' & '+roundfloat(wd_hot.gr,3)+' & '+roundfloat(ms_binary.gr_combined,3)+' & '+roundfloat(k_wd_binary.gr_combined,3)+' & '+roundfloat(m_wd_binary.gr_combined,3)+'\\\\\n')
    t.write('$(r-i)$       & '+roundfloat(kstar.ri,3)+' & '+roundfloat(mstar.ri,3)+' & '+roundfloat(wd_cool.ri,3)+' & '+roundfloat(wd_mid.ri,3)+' & '+roundfloat(wd_hot.ri,3)+' & '+roundfloat(ms_binary.ri_combined,3)+' & '+roundfloat(k_wd_binary.ri_combined,3)+' & '+roundfloat(m_wd_binary.ri_combined,3)+'\\\\\n')
    t.write('$(J-H)$       & '+roundfloat(kstar.JH,3)+' & '+roundfloat(mstar.JH,3)+' & '+roundfloat(wd_cool.JH,3)+' & '+roundfloat(wd_mid.JH,3)+' & '+roundfloat(wd_hot.JH,3)+' & '+roundfloat(ms_binary.JH_combined,3)+' & '+roundfloat(k_wd_binary.JH_combined,3)+' & '+roundfloat(m_wd_binary.JH_combined,3)+'\\\\\n')
    t.write('$(H-Ks)$      & '+roundfloat(kstar.HKs,3)+' & '+roundfloat(mstar.HKs,3)+' & '+roundfloat(wd_cool.HKs,3)+' & '+roundfloat(wd_mid.HKs,3)+' & '+roundfloat(wd_hot.HKs,3)+' & '+roundfloat(ms_binary.HKs_combined,3)+' & '+roundfloat(k_wd_binary.HKs_combined,3)+' & '+roundfloat(m_wd_binary.HKs_combined,3)+'\\\\\n')
    t.write('$(J-Ks)$      & '+roundfloat(kstar.JKs,3)+' & '+roundfloat(mstar.JKs,3)+' & '+roundfloat(wd_cool.JKs,3)+' & '+roundfloat(wd_mid.JKs,3)+' & '+roundfloat(wd_hot.JKs,3)+' & '+roundfloat(ms_binary.JKs_combined,3)+' & '+roundfloat(k_wd_binary.JKs_combined,3)+' & '+roundfloat(m_wd_binary.JKs_combined,3)+'\\\\\n')
    t.write('$m_{B}$       & '+roundfloat(kstar.mB,3)+' & '+roundfloat(mstar.mB,3)+' & '+roundfloat(wd_cool.mB,3)+' & '+roundfloat(wd_mid.mB,3)+' & '+roundfloat(wd_hot.mB,3)+' & '+roundfloat(ms_binary.mB_combined,3)+' & '+roundfloat(k_wd_binary.mB_combined,3)+' & '+roundfloat(m_wd_binary.mB_combined,3)+'\\\\\n')
    t.write('$m_{V}$       & '+roundfloat(kstar.mV,3)+' & '+roundfloat(mstar.mV,3)+' & '+roundfloat(wd_cool.mV,3)+' & '+roundfloat(wd_mid.mV,3)+' & '+roundfloat(wd_hot.mV,3)+' & '+roundfloat(ms_binary.mV_combined,3)+' & '+roundfloat(k_wd_binary.mV_combined,3)+' & '+roundfloat(m_wd_binary.mV_combined,3)+'\\\\\n')
    t.write('$m_{g}$       & '+roundfloat(kstar.mg,3)+' & '+roundfloat(mstar.mg,3)+' & '+roundfloat(wd_cool.mg,3)+' & '+roundfloat(wd_mid.mg,3)+' & '+roundfloat(wd_hot.mg,3)+' & '+roundfloat(ms_binary.mg_combined,3)+' & '+roundfloat(k_wd_binary.mg_combined,3)+' & '+roundfloat(m_wd_binary.mg_combined,3)+'\\\\\n')
    t.write('$m_{r}$       & '+roundfloat(kstar.mr,3)+' & '+roundfloat(mstar.mr,3)+' & '+roundfloat(wd_cool.mr,3)+' & '+roundfloat(wd_mid.mr,3)+' & '+roundfloat(wd_hot.mr,3)+' & '+roundfloat(ms_binary.mr_combined,3)+' & '+roundfloat(k_wd_binary.mr_combined,3)+' & '+roundfloat(m_wd_binary.mr_combined,3)+'\\\\\n')
    t.write('$m_{i}$       & '+roundfloat(kstar.mi,3)+' & '+roundfloat(mstar.mi,3)+' & '+roundfloat(wd_cool.mi,3)+' & '+roundfloat(wd_mid.mi,3)+' & '+roundfloat(wd_hot.mi,3)+' & '+roundfloat(ms_binary.mi_combined,3)+' & '+roundfloat(k_wd_binary.mi_combined,3)+' & '+roundfloat(m_wd_binary.mi_combined,3)+'\\\\\n')
    t.write('$m_{J}$       & '+roundfloat(kstar.mJ,3)+' & '+roundfloat(mstar.mJ,3)+' & '+roundfloat(wd_cool.mJ,3)+' & '+roundfloat(wd_mid.mJ,3)+' & '+roundfloat(wd_hot.mJ,3)+' & '+roundfloat(ms_binary.mJ_combined,3)+' & '+roundfloat(k_wd_binary.mJ_combined,3)+' & '+roundfloat(m_wd_binary.mJ_combined,3)+'\\\\\n')
    t.write('$m_{H}$       & '+roundfloat(kstar.mH,3)+' & '+roundfloat(mstar.mH,3)+' & '+roundfloat(wd_cool.mH,3)+' & '+roundfloat(wd_mid.mH,3)+' & '+roundfloat(wd_hot.mH,3)+' & '+roundfloat(ms_binary.mH_combined,3)+' & '+roundfloat(k_wd_binary.mH_combined,3)+' & '+roundfloat(m_wd_binary.mH_combined,3)+'\\\\\n')
    t.write('$m_{Ks}$      & '+roundfloat(kstar.mKs,3)+' & '+roundfloat(mstar.mKs,3)+' & '+roundfloat(wd_cool.mKs,3)+' & '+roundfloat(wd_mid.mKs,3)+' & '+roundfloat(wd_hot.mKs,3)+' & '+roundfloat(ms_binary.mKs_combined,3)+' & '+roundfloat(k_wd_binary.mKs_combined,3)+' & '+roundfloat(m_wd_binary.mKs_combined,3)+'\\\\\n')
    t.write('\hline\n')
    t.write('$m_{V,corr}$  & '+roundfloat(kstar.mV_corr,3)+' & '+roundfloat(mstar.mV_corr,3)+' & '+roundfloat(wd_cool.mV_corr,3)+' & '+roundfloat(wd_mid.mV_corr,3)+' & '+roundfloat(wd_hot.mV_corr,3)+' & '+roundfloat(ms_binary.mV_combined_corr,3)+' & '+roundfloat(k_wd_binary.mV_combined_corr,3)+' & '+roundfloat(m_wd_binary.mV_combined_corr,3)+'\\\\\n')
    t.write('$(B-V)_{corr}$& '+roundfloat(kstar.BV_corr,3)+' & '+roundfloat(mstar.BV_corr,3)+' & '+roundfloat(wd_cool.BV_corr,3)+' & '+roundfloat(wd_mid.BV_corr,3)+' & '+roundfloat(wd_hot.BV_corr,3)+' & '+roundfloat(ms_binary.BV_combined_corr,3)+' & '+roundfloat(k_wd_binary.BV_combined_corr,3)+' & '+roundfloat(m_wd_binary.BV_combined_corr,3)+'\\\\\n')
    t.write('$m_{g,corr}$  & '+roundfloat(kstar.mg_corr,3)+' & '+roundfloat(mstar.mg_corr,3)+' & '+roundfloat(wd_cool.mg_corr,3)+' & '+roundfloat(wd_mid.mg_corr,3)+' & '+roundfloat(wd_hot.mg_corr,3)+' & '+roundfloat(ms_binary.mg_combined_corr,3)+' & '+roundfloat(k_wd_binary.mg_combined_corr,3)+' & '+roundfloat(m_wd_binary.mg_combined_corr,3)+'\\\\\n')
    t.write('$m_{r,corr}$  & '+roundfloat(kstar.mr_corr,3)+' & '+roundfloat(mstar.mr_corr,3)+' & '+roundfloat(wd_cool.mr_corr,3)+' & '+roundfloat(wd_mid.mr_corr,3)+' & '+roundfloat(wd_hot.mr_corr,3)+' & '+roundfloat(ms_binary.mr_combined_corr,3)+' & '+roundfloat(k_wd_binary.mr_combined_corr,3)+' & '+roundfloat(m_wd_binary.mr_combined_corr,3)+'\\\\\n')
    t.write('$m_{i,corr}$  & '+roundfloat(kstar.mi_corr,3)+' & '+roundfloat(mstar.mi_corr,3)+' & '+roundfloat(wd_cool.mi_corr,3)+' & '+roundfloat(wd_mid.mi_corr,3)+' & '+roundfloat(wd_hot.mi_corr,3)+' & '+roundfloat(ms_binary.mi_combined_corr,3)+' & '+roundfloat(k_wd_binary.mi_combined_corr,3)+' & '+roundfloat(m_wd_binary.mi_combined_corr,3)+'\\\\\n')
    t.write('$(g-r)_{corr}$& '+roundfloat(kstar.gr_corr,3)+' & '+roundfloat(mstar.gr_corr,3)+' & '+roundfloat(wd_cool.gr_corr,3)+' & '+roundfloat(wd_mid.gr_corr,3)+' & '+roundfloat(wd_hot.gr_corr,3)+' & '+roundfloat(ms_binary.gr_combined_corr,3)+' & '+roundfloat(k_wd_binary.gr_combined,3)+' & '+roundfloat(m_wd_binary.gr_combined_corr,3)+'\\\\\n')
    t.write('$(r-i)_{corr}$& '+roundfloat(kstar.ri_corr,3)+' & '+roundfloat(mstar.ri_corr,3)+' & '+roundfloat(wd_cool.ri_corr,3)+' & '+roundfloat(wd_mid.ri_corr,3)+' & '+roundfloat(wd_hot.ri_corr,3)+' & '+roundfloat(ms_binary.ri_combined_corr,3)+' & '+roundfloat(k_wd_binary.ri_combined,3)+' & '+roundfloat(m_wd_binary.ri_combined_corr,3)+'\\\\\n')
    t.write('\hline\n')
    t.write('\end{tabular}\n')
    t.write('\end{table}\n')
    
    t.close()

def option_or_none(option):
    
    try:
        value = float(option)
        
    except ValueError:
        value = None
    
    return value
    
if __name__ == '__main__':
    
    if len(argv) == 1:
        
        output_path = raw_input('Please enter the path to the output file: ')
        dm = float(raw_input('Please enter the distance modulus: '))
        AV = option_or_none(raw_input('Please enter the extinction Av [or None]: '))
        EBV = option_or_none(raw_input('Please enter the reddening, E(B-V) [or None]: '))
        Ag = option_or_none(raw_input('Please enter the extinction A(SDSS-g) [or None]: '))
        Ar = option_or_none(raw_input('Please enter the extinction A(SDSS-r) [or None]: '))
        Ai = option_or_none(raw_input('Please enter the extinction A(SDSS-i) [or None]: '))
        
    else:
        output_path = argv[1]
        dm = float(argv[2])
        AV = option_or_none(argv[3])
        EBV = option_or_none(argv[4])
        Ag = option_or_none(argv[5])
        Ar = option_or_none(argv[6])
        Ai = option_or_none(argv[7])
    
    extinction = { 'AV': AV, 'Ag': Ag, 'Ar': Ar, 'Ai': Ai, 'EBV': EBV }
    
    predict_photometry(output_path,dm,extinction)